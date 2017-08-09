%------------------------------
% Top-level simulation function for the PolarCube mission
%------------------------------
% Franklin Hinckley
% 30 January 2016
%------------------------------
% Inputs:
%   SV (struct): spacecraft properties
%       .I (3 x 3 double matrix): inertia tensor in body 
%           coordinates [kg m^2]
%       .Im (3 x 3 double matrix): mirror inertia tensor in mirror
%           coordinates [kg m^2]
%       .Gs (3 x 3 double matrix): matrix of wheel spin axes
%       .Gt1 (3 x 3 double matrix): matrix of wheel first transverse axes
%       .Gt2 (3 x 3 double matrix): matrix of wheel second transverse axes
%       .Js (3 x 3 double matrix): vector of wheel inertias about the
%           spin axis [kg m^2]
%       .Jt1 (3 x 1 double vector): vector of wheel inertias about the
%           first transverse axis [kg m^2]
%       .Jt2 (3 x 1 double vector): vector of wheel inertias about the
%           second transverse axis [kg m^2]
%       .mu0 (3 x 1 double vector): initial mirror angle [rad]
%       .muDot (function handle): function handle for mirror
%           rotation rate [rad/s]
%       .muDoubDot (function handle): function handle for mirror
%           acceleration rate [rad/s^2]
%       .ms0 (3 x 1 double vector): mirror short axis at initial time
%       .ml0 (3 x 1 double vector): mirror long axis at initial time
%       .mr0 (3 x 1 double vector): mirror rotation axis at initial time
%       .BM (3 x 3 double vector): mirror to body rotation matrix
%           at initial time
%   SIM (struct): simulation configuration
%       .simTime (double): simulation duration [sec]
%       .simFreq (double): simulation rate for dynamics [Hz]
%       .measFreq (double): rate for measurement update [Hz]
%       .controlFreq (double): rate for main control [Hz]
%       .L (3 x 1 array of function handles): defines known external
%           torques about each body axis
%       .dL (3 x 1 array of function handles): defines unknown external
%           torques about each body axis
%       .omega0 (3 x 1 double array): initial SV angular rates [rad/s]
%       .MRP0 (3 x 1 double array): initial attitude in MRPs
%       .MRPRef (3 x 1 array of function handles): 
%           reference trajectory in MRPs
%       .omegaRef (3 x 1 array of function handles):
%           reference trajectory angular rate [rad/s]
%       .omegaDotRef (3 x 1 array of function handles): 
%           reference trajectory angular acceleration [rad/s^2]
%       .RWdyn (double): flag for reaction wheel dynamics
%       .RWdmode (double): flag for raction wheel control
%       .integralControl (double): flag for integral control
%       .mirrorControl (double): flag for mirror model control
% Outputs:
%   SIM (struct): simulation results
%   SV (struct): final SV state
% Dependencies:
%   MRPdiff
%------------------------------

function [SIM,SV] = ADCS_Sim(SIM,SV,code_source)
mex control.c GPS/vector3D.c svd.c
%% Set up simulation
% Generate time series
tvec = 0 : 1/SIM.simFreq : SIM.simTime;

% Determine loop indices
controlInd = SIM.simFreq/SIM.controlFreq;
measInd = SIM.simFreq/SIM.measFreq;
gyroInd = SIM.simFreq/SV.gyro1.sampFreq;
STInd = SIM.simFreq/SV.starTrk.sampFreq;

% Initialize output matrices
len = length(tvec);

MRP = NaN*ones(3,len);
refMRP = NaN*ones(3,len);
omega = NaN*ones(3,len);
h = NaN*ones(3,len);
hSV = NaN*ones(3,len);
hVSCMG = NaN*ones(3,len);
hM = NaN*ones(3,len);
T = NaN*ones(1,len);
Pth = NaN*ones(1,len);
Pan = NaN*ones(1,len);

pointReq = NaN*ones(1,len);
knowReq = NaN*ones(1,len);

N = length(SV.Omega);
Omega = NaN*ones(N,len);
OmegaDot = NaN*ones(N,len);

mu = NaN*ones(len,1);
muDot = NaN*ones(len,1);

currMode = zeros(len,1);

% Set initial conditions
MRP(:,1) = SIM.MRP0;
refMRP(:,1) = zeros(3,1);
omega(:,1) = SIM.omega0;
Omega(:,1) = SV.Omega;

SV.mu = SV.mu0;

SV.us = zeros(N,1);
SV.ug = zeros(N,1);

SIM.omegaDot = [0;0;0];

MRPRef = zeros(3,1);
MRPRefPrev  = zeros(3,1);
omegaRefPrev = 0;

sigInt = [0;0;0];
SIM.z = zeros(3,len);

Lsim = zeros(3,len);

wStar1 = zeros(3,len);
wStar2 = zeros(3,len);
wStarFilt = zeros(3,25);
bTrue = zeros(3,len);

bBStar1 = zeros(3,len);
bBStar2 = zeros(3,len);

PFR = NaN*ones(3,len);

sunAngle = zeros(len,1);
excAngle = zeros(len,1);
txAngle = zeros(len,1);

Lt = zeros(3,1);
SV.Lt = zeros(3,1);

Bdot = zeros(3,len);
SV.Bdot = zeros(3,1);

BbodyS = zeros(3,len);

tSCI = NaN;
z = zeros(3,1);
errOmega0 = zeros(3,1);

mode = 0;
modePrev = 0;
tMode = -1e5;

%% Initialize Control Filter
% Set matrix for previous control values
SV.Lprev = zeros(3,2);

%% Initialize Attitude Determination Filter
% Define standard deviations for rate and bias rate
std_w = mean(SV.gyro1.noise);
std_b = std_w*1e-3;

% Define initial covariance
P_sig = 0.175;
P_w = 0.005;
P = diag([P_sig P_sig P_sig P_w P_w P_w]);

% Define measurement error covariance matrix
dMRP = 100*4.84814e-6; %100 arc-sec
R = dMRP*eye(3);

% Process noise
% Q = [std_w^2*eye(3) zeros(3)       ;...
%      zeros(3)       std_b^2*eye(3)];
% Q = [5e-5*eye(3) zeros(3)       ;...
%      zeros(3)       1e-16*eye(3)]; % from O'Keefe and Schaub
Q = [5e-6*eye(3) zeros(3)       ;...
     zeros(3)       1e-16*eye(3)]; 
% Q = [std_w^2*eye(3) zeros(3)       zeros(3);...
%      zeros(3)       std_b^2*eye(3) zeros(3);...
%      zeros(3)       zeros(3)       std_b^2*eye(3)]*1;

%% Define [B(sigma)] matrix
B = @(mrp)((1-(mrp'*mrp)))*eye(3)+(2*tilde(mrp))+(2*(mrp*mrp'));

%% Define f,g matrices
f = @(mrp,wbar)[0.25*B(mrp)*wbar; zeros(3,1)];
g = @(mrp,etaW,etaB)[-0.25*B(mrp)*etaW; etaB];
% f = @(mrp,wbar)[0.25*B(mrp)*wbar; zeros(6,1)];
% g = @(mrp,etaW,etaB)[-0.25*B(mrp)*etaW; etaB; etaB];

%% Define F, G matrices
F = @(sBar,wBar)[0.5*(sBar*wBar' - wBar*sBar' - tilde(wBar) + sBar'*wBar*eye(3)) -0.25*B(sBar);...
    zeros(3) zeros(3)];
G = @(sBar)[-0.25*B(sBar) zeros(3);...
    zeros(3) eye(3)];
% F = @(sBar,wBar)[0.5*(sBar*wBar' - wBar*sBar' - tilde(wBar) + sBar'*wBar*eye(3)) -0.125*B(sBar) -0.125*B(sBar);...
%     zeros(6,3) zeros(6,6)];
% G = @(sBar)[-0.25*B(sBar) zeros(3,6);...
%     zeros(6,3) eye(6)];

%% Define H matrix
H = [eye(3) zeros(3)];
%H = [eye(3) zeros(3) zeros(3)];

X = zeros(6,1);
X(4:6) = 0.75*SV.gyro1.bias;
xhat = zeros(6,1);
numST = 0;
wBar = 0;
MRP_est = zeros(3,len);
b_est = zeros(3,len);
b_est2 = zeros(3,len);
wBarS = zeros(3,len);
PS = zeros(6,6,len);

%% Main loop
for ii = 2:len
    %% Evaluate external torques
    % Known external torques
    L = [SIM.L{1}(tvec(ii));SIM.L{2}(tvec(ii));SIM.L{3}(tvec(ii))];
    % Unknown external torques
    dL = [SIM.dL{1}(tvec(ii));SIM.dL{2}(tvec(ii));SIM.dL{3}(tvec(ii))];
    
    % Save total external torque
    Lsim(:,ii) = L + dL;
    
    %% Determine current Julian Date
    JD = SIM.JD0 + (ii/86400/SIM.simFreq);
    
    %% Evaluate magnetic field model
    % Determine index into propagated position array for current sim time
    posInd = floor(ii/SIM.simFreq) + 1;
    % Evaluate magnetic field model
    [SIM.Bbody] = BfieldComp(MRP(:,ii-1), SIM.prop.R(1,posInd), ...
        SIM.prop.R(2,posInd), SIM.prop.R(3,posInd), JD);
    BbodyS(:,ii) = SIM.Bbody;
    
    %% Evaluate potential pointing targets
    % Determine the sun vector
    sunVec = sun_vector(JD);

    % Determine COM and science targets
    COMvec = get_target_vector(SIM.latCOM, SIM.lonCOM, JD);
    SCIvec = get_target_vector(SIM.latSCI, SIM.lonSCI, JD);
    
    % Get current attitude as a DCM
    attDCM = mrp2dcm(MRP(:,ii-1));
    
    % Get sun vector in body coordinates
    sunVec_body = attDCM*sunVec;
    
    % Log angle between solar panel normal and sun vector
    SPaxis_body = [-sqrt(2)/2 -sqrt(2)/2 0]';
    sunDot = dot(sunVec_body,SPaxis_body);
    sunAngle(ii) = acos(sunDot);
    
    % Log star camera sun exclusion angle
    SCaxis_body = [0 1 0]';
    excDot = dot(sunVec_body,SCaxis_body);
    excAngle(ii) = acos(excDot);
    
    % Log S-Band transmit antenna pointing error
    TXaxis_body = [1 0 0]';
    GSaxis = COMvec - SIM.prop.R(:,posInd);
    GSaxis = GSaxis/norm(GSaxis);
    GSaxis_body = attDCM*GSaxis;
    txDot = dot(TXaxis_body,GSaxis_body);
    txAngle(ii) = acos(txDot);
    
    %% Determine current mode
    % Check ranges to targets
    COMrange = norm(SIM.prop.R(:,posInd) - COMvec);
    SCIrange = norm(SIM.prop.R(:,posInd) - SCIvec);
    
    % Check detumble time
    if ((ii/SIM.simFreq) > SIM.detTime)
        % Check target ranges
        if COMrange < SIM.COMrngLim
            modeTmp = 2;
        elseif SCIrange < SIM.SCIrngLim
            modeTmp = 3;
        elseif SCIrange < 3000
            modeTmp = 4;
        else
            modeTmp = 1;
        end
    else
        modeTmp = 0;
    end
    
%     % Check for mode switch
%     if mode ~= modePrev
%         tMode = tvec(ii);
%         mode = 4;
%     elseif (tvec(ii) - tMode) < 120
%         mode = 4;
%     else
%         mode = modeTmp;
%     end
    mode = modeTmp;
    
    % Reset previous mode
    modePrev = mode;
    
    %% Set parameters for current mode
    % Set gains, flags
    [SV, SIM] = modeParamSet(SV, SIM, mode);
    
    % Define requirements for current mode 
    switch mode
        case 0 % Detumble
            pointReq(ii) = NaN;
            knowReq(ii) = 2.5*pi/180;
        case 1 % Bus priority (sun-pointing)
            pointReq(ii) = 14*pi/180;
            knowReq(ii) = 2.5*pi/180;
        case 2 % S-Band COM
            pointReq(ii) = 15*pi/180;
            knowReq(ii) = 2.5*pi/180;
        case 3 % Science
            pointReq(ii) = 2.3*pi/180;
            knowReq(ii) = 0.23*pi/180;
        case 4 % Pre-science
            pointReq(ii) = 2.3*pi/180;
            knowReq(ii) = 0.23*pi/180;
    end   
    
    % Set science time
    if mode ~= 3
        tSCI = NaN;
    else
        if isnan(tSCI)
            tSCI = tvec(ii);
        end
    end
    
    % Save mode
    currMode(ii) = mode;
    
    %% Evaluate reference trajectory
    if (ii==26) || (mod(ii, controlInd) == 0)
        [MRPRef, omegaRef, omegaDotRef] = ...
            modeRefSet(modeTmp, sunVec, COMvec, SIM.prop.R(:,posInd));
    end
    
    % Save current reference
    refMRP(:,ii) = MRPRef;
    
    %% Simulate sensor response
    % Gyroscopes
    [wStar1(:,ii),wStar2(:,ii),bTrue1,bTrue2] = ...
        A3G4250_Model(omega(:,ii-1),SV,tvec(ii));
    %wStarAvg(:,ii) = mean([wStar1, wStar2],2);
    
    % Magnetometer
    [bBStar1(:,ii),bBStar2(:,ii)] = LIS3MDL_Model(SIM.Bbody,SV);
    
    % Star tracker
    if mod(ii,STInd) == 0
        sStar = ST_Model(MRP(:,ii-1),SV);
    end
    
    %% Attitude Determination Filter and Control
    if ii > (50 * gyroInd) % 25 gyroscope samples
        if mod(ii, gyroInd) == 0
        % Average previous gyro measurements
        if mode ~=3
%             wStar1avg = mean(wStar1(:,ii-24:ii),2);
%             wStar2avg = mean(wStar2(:,ii-24:ii),2);
            wStar1avg = mean(wStar1(:,ii-9:ii),2);
            wStar2avg = mean(wStar2(:,ii-9:ii),2);
        else
            wStar1avg = mean(wStar1(:,ii-3:ii),2);
            wStar2avg = mean(wStar2(:,ii-3:ii),2);
        end
        %wStar = wStarP(:,ii);
%         wStar1avg = wStar1(:,ii);
%         wStar2avg = wStar1(:,ii);
        
        wStarAvg = mean([wStar1avg wStar2avg],2);
        
        bTrue(:,ii) = mean([bTrue1 bTrue2],2);

        % Low-pass filter
%         a = 1;
%         b = [1/4 1/4 1/4 1/4];
%         wStar1avg = filter(b,a,wStar1(:,ii-24:ii));
%         wStar1avg = wStar1avg(:,end);
%         wStar2avg = filter(b,a,wStar2(:,ii-24:ii));
%         wStar2avg = wStar2avg(:,end);
        
        % Check for mode transition
%         if modeTrans == 1
%             wStarFilt = wStarAvg(:,ii);
%         else
%             [~,lw] = size(wStarFilt);
%             if lw < 25
%                 wStarFilt = [wStarFilt wStarAvg(:,ii)];
%             else
%                 wStarFilt = [wStarFilt(:,2:24) wStarAvg(:,ii)];
%             end
%         end
%         wStar = mean(wStarFilt,2);
        

        %Pdot = Fm*P + P*Fm' + Gm*Q*Gm';

        % Integrate state
        gyroUpdateInterval = (1/SV.gyro1.sampFreq);
        [X, P, wBar] = ADCS_EKF_timeUpdate(X, P, wStarAvg, gyroUpdateInterval);
        
        % Compute measurement options
        if (mod(ii,SIM.simFreq*5) == 0) &&... % every 5 sec
                (ii > 10*60*SIM.simFreq) &&...  % after ten minutes
                (norm(omega(:,ii-1)) < (0.25*pi/180)) % rate < thres              
            sStar = MRP(:,ii-1) + dMRP*randn(3,1);
            
            if numST < 10
                SIM.RWmode = 0;
                SV.OmegaDot = zeros(3,1);
                SV.Lt = zeros(3,1);
                numST = numST + 1;
                R = 10*dMRP*eye(3);
            else
                SIM.RWmode = 1;
                R = dMRP*eye(3);
            end

            % Filter measurement update
            [X, P, y, PFR(:,ii)] = ADCS_EKF_measUpdate(X, P, sStar, R);

%             % Post-fit residuals
%             PFR(:,ii) = y - H*(K*y);
        end 
    
    %% Compute control inputs
    % Compute error MRP
    errMRP = MRPdiff(MRPRef,X(1:3));
    
    % Compute angular rate error
    errOmega = wBar - omegaRef;
    
    % Store initial delta-omega
    if ii == (26*gyroInd)
        errOmega0 = errOmega;
    end
    
    %% Integral term
    % Accumulate z-term
    if (ii==(26*gyroInd)) || (mod(ii, measInd) == 0)
        K = SV.K;
        I_VSCMG = 0;
        for kk = 1 : 3
            tmpS = SV.Js(kk)  * SV.Gs(:,kk)  * SV.Gs(:,kk)';
            tmpT = SV.Jt1(kk) * SV.Gt1(:,kk) * SV.Gt1(:,kk)';
            tmpG = SV.Jt2(kk) * SV.Gt2(:,kk) * SV.Gt2(:,kk)';
            I_VSCMG = I_VSCMG + tmpS + tmpT + tmpG;
        end
        I = SV.I + I_VSCMG;
        sigInt = sigInt + errMRP*(1/SIM.measFreq);
        z = K*sigInt + I*(errOmega-errOmega0);  
    end
    SIM.z(:,ii) = z;
    
%% Control
    % Check if control is enabled
    if (SIM.RWmode == 1) && (ii > 1*SIM.controlFreq)
        % Check if time to update control
        if (ii==26) || (mod(ii, controlInd) == 0)
                   if code_source == 1
                [OmegaDotD, SV] = ADCS_Control(errMRP,z,...
                    omega(:,ii-1),errOmega,omegaRef,omegaDotRef,L,SV,SIM);
                   elseif code_source == 2
                    [OmegaDotD, SV] = ADCS_Control_C(errMRP,z,...
                    omega(:,ii-1),errOmega,omegaRef,omegaDotRef,L,SV,SIM);
                   end
            SV.OmegaDot = OmegaDotD; 
        end
    end
    
    end % end pre-filter loop
    end
    
    % Log filter output
    MRP_est(:,ii) = X(1:3);
    b_est(:,ii) = X(4:6);
    wBarS(:,ii) = wBar;  
    PS(:,:,ii) = P;
    
    Bdot(:,ii) = SV.Bdot;
    Lt(:,ii) = SV.Lt;
               
    %% Integrate Dynamics
    LI = L + dL + SV.Lt;
    [MRP(:,ii), omega(:,ii),SV,SIM] = ADCS_Integrate(MRP(:,ii-1),...
        omega(:,ii-1),LI,SIM,SV,1/SIM.simFreq,tvec(ii)-tSCI);  
    
    % Save torques
    SIM.us(:,ii) = SV.us;
    
    % Check for switching condition
    if norm(MRP(:,ii)) > 1
        MRP(:,ii) = -MRP(:,ii)/(norm(MRP(:,ii))^2);
    end
        
    % Save results
    Omega(:,ii) = SV.Omega;
    OmegaDot(:,ii) = SV.OmegaDot;
    mu(ii) = mod(SV.mu,2*pi);
    muDot(ii) = SV.muDotA;
    
    %% Simulation checks
    % Compute angular momentum, kinetic energy, power
    hVSCMGtmp = 0;
    Ttmp = 0;
    Ptmp = 0;
    for jj = 1 : N        
        ws = SV.Gs(:,jj)'*omega(:,ii);
        wt = SV.Gt1(:,jj)'*omega(:,ii);
        wg = SV.Gt2(:,jj)'*omega(:,ii);
        hVSCMGtmp = hVSCMGtmp + ...
            SV.Js(jj)*(ws + Omega(jj,ii))*SV.Gs(:,jj) +...
            SV.Jt1(jj)*wt*SV.Gt1(:,jj) +...
            SV.Jt2(jj)*(wg)*SV.Gt2(:,jj);
        Ttmp = Ttmp +...
            SV.Js(jj)*((Omega(jj,ii)+ws)^2) +...
            SV.Jt1(jj)*(wt^2) +...
            SV.Jt2(jj)*((+wg)^2);
        Ptmp = Ptmp +...
            SV.Omega(jj)*SV.us(jj);          
    end  
    wMN = SV.muDot(tvec(ii))*SV.BM(:,3) + omega(:,ii);
    hVSCMG(:,ii) = hVSCMGtmp;
    hSV(:,ii) = SV.I*omega(:,ii);
    hM(:,ii) = SV.BM*SV.Im*SV.BM' * (omega(:,ii)+SV.muDot(tvec(ii))*SV.BM(:,3));
    h(:,ii) = hSV(:,ii) + hVSCMG(:,ii) + hM(:,ii);
    T(ii) = 0.5*Ttmp + 0.5*omega(:,ii)'*SV.I*omega(:,ii) + 0.5*wMN'*(SV.BM*SV.Im*SV.BM')*wMN;
    Pth(ii) = omega(:,ii)'*LI + Ptmp;
    Pan(ii) = (T(ii) - T(ii-1))/(1/SIM.simFreq);  
    
    %% Progress bar
    if (ii == 2) || (mod(ii,5*SIM.simFreq) == 0)
        progressbar(ii/len)
    end
end

% Clean up progress bar
progressbar(1)

% Assign outputs
SIM.MRP = MRP;
SIM.refMRP = refMRP;
SIM.omega = omega;
SIM.t = tvec;
SIM.h = h;
SIM.hSV = hSV;
SIM.hVSCMG = hVSCMG;
SIM.hM = hM;
SIM.T = T;
SIM.Pth = Pth;
SIM.Pan = Pan;

SIM.pointReq = pointReq;
SIM.knowReq = knowReq;

SIM.Omega = Omega;
SIM.OmegaDot = OmegaDot;
SIM.mu = mu;
SIM.muDot = muDot;

SIM.Lsim = Lsim;

SIM.MRP_est = MRP_est;
SIM.b_est = b_est;
SIM.wBar = wBarS;
SIM.P = PS;
SIM.PFR = PFR;

SIM.wStar1 = wStar1;
SIM.wStar2 = wStar2;

SIM.sunAngle = sunAngle;
SIM.excAngle = excAngle;
SIM.txAngle = txAngle;

SIM.mode = currMode;

SIM.Bdot = Bdot;
SIM.Bbody = BbodyS;

SIM.bTrue = bTrue;

end

    
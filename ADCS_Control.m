%------------
% Computes the control vector for the PolarCube satellite
%------------
% Franklin Hinckley
% 30 September 2015
%------------
% Inputs:
%   errMRP (3 x 1 double vector): spacecraft attitude error in MRPs
%   z (3 x 1 double vector): integral term [MRP-sec]
%   omega (3 x 1 double vector): spacecraft angular rate [rad/s]
%   errOmega (3 x 1 double vector): spacecraft angular rate error [rad/s]
%   omegaRef (3 x 1 double vector): spacecraft angular rate reference [rad/s]
%   omegaDotRef (3 x 1 double vector): spacecraft angular acceleration 
%       reference [rad/s^2]
%   L (3 x 1 double array): known external torques [Nm]
%   SV (struct): satellite properties
%   SIM (struct): simulation properties
% Outputs:
%   OmegaDot (N x 1 double vector): wheel angular accelerations [rad/s^2]
%   SV (struct): updated SV properties
%------------

function [OmegaDot,SV] = ADCS_Control(errMRP,z,omega,errOmega,omegaRef,...
    omegaDotRef,L,SV,SIM)

% Define control parameters
K = SV.K;
P = SV.P;
Km = SV.Km;

% Check for integral control
if SIM.integralControl == 1
    Ki = SV.Ki;
else
    Ki = diag([0;0;0]);
end

% Get number of VSCMGs
[~,N] = size(SV.Gs);

% Pull out G Matrices
Gs = SV.Gs;
Gt1 = SV.Gt1;
Gt2 = SV.Gt2;
Gt = SV.Gt;

% Pull out J vectors
Js = SV.Js;
Jt1 = SV.Jt1;
Jt2 = SV.Jt2;

% Compute summation of VSCMG inertias
I_VSCMG = 0;
for ii = 1 : N
    tmpS = Js(ii) * Gs(:,ii) * Gs(:,ii)';
    tmpT = Jt1(ii) * Gt1(:,ii) * Gt1(:,ii)';
    tmpG = Jt2(ii) * Gt2(:,ii) * Gt2(:,ii)';
    I_VSCMG = I_VSCMG + tmpS + tmpT + tmpG;
end

% Compute total inertia
I = SV.I + I_VSCMG;

% Compute velocity projections
omegaS=zeros(N,1);
omegaT=zeros(N,1);
omegaG=zeros(N,1);
for ii = 1 : N
    omegaS(ii) = Gs(:,ii)'*omega;
    omegaT(ii) = Gt1(:,ii)'*omega;
    omegaG(ii) = Gt2(:,ii)'*omega;
end

% Mirror angle error
dm = SV.mu - SV.mu0;
ms =  cos(dm)*SV.ms0 + sin(dm)*SV.ml0;
ml = -sin(dm)*SV.ms0 + cos(dm)*SV.ml0;
mr = SV.mr0;
BM = [ms ml mr];

% Mirror rotation model
ImB = BM*SV.Im*BM';
I = I + ImB;
muDot = 1.3*2*pi;
Im = SV.Im; 
omegaMB = muDot*SV.BM(:,3);
% fomegaM = tilde(omegaMB)*BM*Im*BM'*(omega+omegaMB)-...
%         BM*Im*BM'*tilde(omegaMB)*omega+...
%         BM*Im*BM'*muDoubDot*SV.BM(:,3)+...
%         tilde(omega)*BM*Im*BM'*omegaMB;
%fomegaM = tilde(omegaMB)*BM*Im*BM'*omegaMB;

% % Check for mirror model
% if SIM.mirrorControl == 1
%     LrM = fomegaM;
% else
%     LrM = [0;0;0];
% end

% Compute required torque
LrRW = 0;
hs = zeros(3,1);
for ii = 1:N
    LrRW = LrRW + ...
        Js(ii)*(SV.Omega(ii)*omegaG(ii)*Gt1(:,ii) - ...
                SV.Omega(ii)*omegaT(ii)*Gt2(:,ii));
        hs(ii) = Js(ii)*(SV.Omega(ii) + Gs(:,ii)'*omega);
end

LrTmp = -I*(omegaDotRef - tilde(omega)*omegaRef) + K*errMRP + P*errOmega + ...
    P*Ki*z - (tilde(omegaRef) - tilde(Ki*z))*(I*omega + Gs*hs) + L; 

%% Momentum management
% Wheel speed error
delOmega = SV.Omega - SV.OmegaNom;

% Compute magnetic torque for dumping
B_body = SIM.Bbody;
uStar = -Km*diag(Js)*delOmega;
muStar = -pinv(tilde(B_body)*Gt)*Gs*uStar;
tau_mag = -tilde(B_body)*Gt*muStar;
delta_u = inv(Gs)\(tau_mag - Gs*uStar);

SV.Lt = tau_mag;

% Compute total control input
LrT = LrTmp + uStar + delta_u;
%LrT = LrTmp;

% Filter required torque
% if SV.contFilt == 1
%     Lr = mean([SV.Lprev LrT],2);
%     SV.Lprev = [SV.Lprev(:,2) LrT];
% else
%     Lr = LrT;
% end
Lr = LrT;
% Set deadband
if norm(Lr) < SV.deadband
    Lr = zeros(3,1);
end

% Define D matrix
D0 = zeros(3,N);
for ii = 1:N
    D0(:,ii) = Gs(:,ii)*Js(ii);
end
    
% Define Q matrix
Q = D0;

% Defing weighting matrix
W = eye(N);

% Check for no control
if (SIM.RWmode == 0)
    OmegaDot = zeros(N,1);
else
    % Solve for desired control eta-dot (eta = [Omega; gamma])
    OmegaDot = W*Q'*((Q*W*Q')^(-1))*Lr;
end

end

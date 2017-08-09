%--------------------------------
% Test the EKF
%--------------------------------
% Franklin Hinckley
% 16 December 2015
%--------------------------------

%% Clean up
clearvars
close all
clc

addpath subroutines

%% Define simulation properties
% Define satellite principal inertias
I1 = 4;
I2 = 4;
I3 = 3;

SV = satProps;

% Define initial states
omega0 = [-0.2 0.2 0.1]'*(pi/180); %[rad/s]
MRP0 = [0.3 0.1 -0.5]';

% Define gyro bias
%wb = [-1 2 -3]'*(pi/180)/3600*100; % [rad/sec]
wb = SV.gyro1.bias;

% Define standard deviations for rate and bias rate
%std_w = sqrt(10)*10^(-7);
std_w = mean(SV.gyro1.noise);
%std_b = sqrt(10)*10^(-10);
std_b = std_w*1e-3;

% Define initial filter estimates
MRP_est0 = [0 0 0]';
wb_est0 = [0 0 0]';

% Define initial covariance
P_sig = 0.175;
P_w = 0.005;
%P_w = 0.25*(pi/180);
P = diag([P_sig P_sig P_sig P_w P_w P_w]);

% Define measurement noise
dMRP = 100*4.84814e-6; %100 arc-sec

% Define measurement error covariance matrix
%R = 0.0004*eye(3)*5;
R = dMRP*1*eye(3);

lmb_w = 0.001;

% Process noise
Q = [std_w^2*eye(3) zeros(3);zeros(3) std_b^2*eye(3)]*5;

%Q = 3e-9*eye(6);

%% Define [B(sigma)] matrix
B = @(mrp)((1-(mrp'*mrp)))*eye(3)+(2*tilde(mrp))+(2*(mrp*mrp'));

%% Define f,g matrices
f = @(mrp,b,wStar)[0.25*B(mrp)*(wStar-b); zeros(3,1)];
g = @(mrp,etaW,etaB)[-0.25*B(mrp)*etaW; etaB];

%% Define F, G matrices
F = @(sBar,wBar)[0.5*(sBar*wBar' - wBar*sBar' - tilde(wBar) + sBar'*wBar*eye(3)) -0.25*B(sBar);...
    zeros(3) zeros(3)];
G = @(sBar)[-0.25*B(sBar) zeros(3);...
    zeros(3) eye(3)];

%% Define H matrix
H = [eye(3) zeros(3)];

%% Run EKF
tstepD = 0.01;
tstep = 0.1;
tSim = 0:tstepD:(60*60);

omega = zeros(3,length(tSim));
omega(:,1) = omega0;
wBarS = zeros(3,length(tSim));
MRP = zeros(3,length(tSim));
MRP(:,1) = MRP0;
MRP_est = zeros(3,length(tSim));
b_est = zeros(3,length(tSim));
Ps = zeros(6,6,length(tSim));
sy = zeros(3,length(tSim));
X = [MRP_est0; wb_est0];

wStarP = zeros(3,length(tSim));

wBar = 0;

for ii = 2:length(tSim) 
    % Propagate dynamics
    w1=omega(1,ii-1);
    w2=omega(2,ii-1);
    w3=omega(3,ii-1);
    omegaDot=[(1/I1)*(-(I3-I2)*w2*w3);...
              (1/I2)*(-(I1-I3)*w3*w1);...
              (1/I3)*(-(I2-I1)*w1*w2)];
           
    % Determine new angular velocity
    omega(:,ii) = omega(:,ii-1) + omegaDot*tstepD;
    
    % Propagate true attitude
    mrpDot = body2mrp(omega(:,ii),MRP(:,ii-1));
    %Compute new MRPs
    MRP(:,ii) = MRP(:,ii-1) + mrpDot*tstepD;
    if norm(MRP(:,ii)) > 1
        MRP(:,ii) = -MRP(:,ii)/(norm(MRP(:,ii))^2);
    end
        
    %wStar = omega(:,ii) + wb + std_w*randn(3,1);
    [wStar1,wStar2] = A3G4250_Model(omega(:,ii),SV);
    wStarP(:,ii) = mean([wStar1, wStar2],2);
    
    %X = X + f(X(1:3),bbar,wStar)*tstep;
        
    %if mod(ii,10) == 0
    if ii > 25
        % Average previous gyro measurements
        wStar = mean(wStarP(:,ii-24:ii),2);
        %wStar = wStarP(:,ii);
        
        % Compute state estimates
        bbar = X(4:6);
        wBar = wStar - bbar;
        sBar = X(1:3);
    
        % Propagate state estimate
        eta_w = std_w*randn(3,1);
        eta_b = std_b*randn(3,1);

        fm = f(sBar,bbar,wStar);
        gm = g(X(1:3),eta_w,eta_b);
        Fm = F(sBar,wBar);
        Gm = G(sBar);

        xDot = fm + gm;
        %mrpDot = body2mrp(wBar,sBar);
        %xDot = [mrpDot; zeros(3,1)] + gm;
        A = expm(-[-Fm Gm*Q*Gm'; zeros(6) Fm']*tstepD);
        Phi = A(7:12,7:12)';
        Qk = Phi*A(1:6,7:12);
        %Pdot = Fm*P + P*Fm' + Gm*lmb_w*Gm' + Q;
        Pdot = Fm*P + P*Fm' + Gm*Q*Gm';

        % Integrate state

        X = X + xDot*tstepD;
        P = P + Pdot*tstepD;
        %X = Phi*X;
        %P = Phi*P*Phi' + Qk;

        % Check for MRP switch
        sig = X(1:3);
        s = norm(sig);
        if s > 1
            X = [-sig/(s^2); bbar];
            Sm = 2*s^(-4)*(sig*sig') - s^(-2)*eye(3);
            P = [Sm*P(1:3,1:3)*Sm' Sm*P(1:3,4:6);...
                   P(4:6,1:3)*Sm'     P(4:6,4:6)];
        end

        % Compute measurement options
        if mod(ii,500) == 0 % Every 5 sec
            sStar = MRP(:,ii) + dMRP*randn(3,1);
            y = sStar - X(1:3);
            %y = MRPmult(sStar,-X(1:3));
            sS = norm(sStar);
            if sS > (1/3)
                sStarS = -sStar/(sS^2);
                yS = sStarS - X(1:3);
                %yS = MRPmult(sStarS,-X(1:3));
                if norm(yS) < norm(y)
                    y = yS;
                end
            end
            sy(:,ii) = y;

            % Compute Kalman gain
            K = P*H'*((H*P*H' + R)^(-1));

            % Update state 
            X = X + K*y;
            bbar = X(4:6);
            % Update covariance (Joseph Formulation)
            P = (eye(6) - K*H)*P*((eye(6) - K*H)') + K*R*K';
        end
        %P = (eye(6) - K*H)*P;
        Ps(:,:,ii) = P;

        % Check again for switch
    %     sig = X(1:3);
    %     s = norm(sig);
    %     if s > 1
    %         X = [-sig/(s^2); bbar];
    %         Sm = 2*s^(-4)*(sig*sig') - s^(-2)*eye(3);
    %         P = [Sm*P(1:3,1:3)*Sm' Sm*P(1:3,4:6);...
    %                P(4:6,1:3)*Sm'     P(4:6,4:6)];
    %     end

        
    end
    % Save
    MRP_est(:,ii) = X(1:3);
    b_est(:,ii) = X(4:6);
    wBarS(:,ii) = wBar;
        
end

%% Plot results
tSim = tSim/60; %convert to minutes
figure
plot(tSim,MRP')
title('True Attitude')
xlabel('Simulation Time [min]')
ylabel('MRPs')

% figure
% plot(tSim,omega')

figure
plot(tSim,MRP_est')
title('Estimated Attitude')
xlabel('Simulation Time [min]')
ylabel('MRPs')

figure
plot(tSim,b_est')
title('Gyroscope Bias')
xlabel('Simulation Time [min]')
ylabel('Bias [rad/s]')

figure
plot(tSim,MRP'-MRP_est')
title('Attitude Error')
xlabel('Simulation Time [min]')
ylabel('MRPs')

figure
semilogy(tSim,abs(MRP'-MRP_est'))
title('Attitude Error')
xlabel('Simulation Time [min]')
ylabel('MRPs')

figure
plot(tSim,sy')
title('Measurement Residuals')
xlabel('Simulation Time [min]')
ylabel('Residual')

figure
subplot(3,1,1)
hold on
plot(tSim,b_est(1,:))
plot(tSim,ones(size(tSim))*wb(1),'--k')
hold off
%ylim([-0.0006 0])

subplot(3,1,2)
hold on
plot(tSim,b_est(2,:))
plot(tSim,ones(size(tSim))*wb(2),'--k')
hold off
%ylim([0.0006 0.0012])

subplot(3,1,3)
hold on
plot(tSim,b_est(3,:))
plot(tSim,ones(size(tSim))*wb(3),'--k')
hold off

phiErr = zeros(size(tSim));
for ii = 1:length(tSim)
    mrpErr = MRPdiff(MRP_est(:,ii),MRP(:,ii));
    Cerr = mrp2dcm(mrpErr);
    [phiErr(ii),~] = dcm2prv(Cerr,'ntol',1);
end

figure
semilogy(tSim,phiErr)
hold on
semilogy([tSim(1) tSim(end)],[0.14*pi/180 0.14*pi/180],'-k')
title('\Phi')
xlabel('Simulation Time [s]')
ylabel('Attitude Error [rad]')

figure
subplot(3,1,1)
title('3\sigma Covariance Envelopes')
hold on
plot(tSim,sy(1,:),'.b')
plot(tSim,reshape(3*Ps(1,1,:),1,length(tSim)),'--k')
plot(tSim,reshape(-3*Ps(1,1,:),1,length(tSim)),'--k')
hold off

subplot(3,1,2)
hold on
plot(tSim,sy(2,:),'.b')
plot(tSim,reshape(3*Ps(2,2,:),1,length(tSim)),'--k')
plot(tSim,reshape(-3*Ps(2,2,:),1,length(tSim)),'--k')
hold off

subplot(3,1,3)
hold on
plot(tSim,sy(3,:),'.b')
plot(tSim,reshape(3*Ps(3,3,:),1,length(tSim)),'--k')
plot(tSim,reshape(-3*Ps(3,3,:),1,length(tSim)),'--k')
hold off

figure
hold on
plot(tSim,wBarS)
%plot(
hold off
title('Angular Velocity')

figure
subplot(3,1,1)
title('Angular Velocity')
hold on
plot(tSim,wBarS(1,:))
plot(tSim,omega(1,:),'--k')
hold off

subplot(3,1,2)
hold on
plot(tSim,wBarS(2,:))
plot(tSim,omega(2,:),'--k')
hold off

subplot(3,1,3)
hold on
plot(tSim,wBarS(3,:))
plot(tSim,omega(3,:),'--k')
hold off


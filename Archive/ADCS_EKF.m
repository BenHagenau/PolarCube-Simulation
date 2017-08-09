%------------------------------
% EKF for PolarCube ADCS
% Returns the updae
%------------------------------
% Franklin Hinckley
% 24 January 2017
%------------------------------
%
%------------------------------

function [X, P] = ADCS_EKF(X, P, wStar1, wStar2, sStar)

% Define [B(sigma)] matrix
B = @(mrp)((1-(mrp'*mrp)))*eye(3)+(2*tilde(mrp))+(2*(mrp*mrp'));

% Define f,g matrices
f = @(mrp,b,wStar)[0.25*B(mrp)*(wStar-b); zeros(3,1)];
g = @(mrp,etaW,etaB)[-0.25*B(mrp)*etaW; etaB];

% Define F, G matrices
F = @(sBar,wBar)[0.5*(sBar*wBar' - wBar*sBar' - tilde(wBar) + sBar'*wBar*eye(3)) -0.25*B(sBar);...
    zeros(3) zeros(3)];
G = @(sBar)[-0.25*B(sBar) zeros(3);...
    zeros(3) eye(3)];

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
Pdot = Fm*P + P*Fm' + Gm*Q*Gm';

% Integrate state
X = X + xDot*(1/SV.gyro1.sampFreq);
P = P + Pdot*(1/SV.gyro1.sampFreq);

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
if mod(ii,SIM.simFreq*SV.starTrk.sampFreq) == 0 
    y = sStar - X(1:3);
    sS = norm(sStar);
    if sS > (1/3)
        sStarS = -sStar/(sS^2);
        yS = sStarS - X(1:3);
        if norm(yS) < norm(y)
            y = yS;
        end
    end

    % Compute Kalman gain
    K = P*H'*((H*P*H' + R)^(-1));

    % Update state 
    X = X + K*y;
    % Update covariance (Joseph Formulation)
    P = (eye(6) - K*H)*P*((eye(6) - K*H)') + K*R*K';
end 

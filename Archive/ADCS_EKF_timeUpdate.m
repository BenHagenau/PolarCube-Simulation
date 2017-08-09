%------------------------------
% EKF for PolarCube ADCS
% Evaluates the time update for the EKF
%------------------------------
% Franklin Hinckley
% 24 January 2017
%------------------------------
%
%------------------------------

function [X, P, w_bar] = ADCS_EKF_timeUpdate(X, P, wStar, gyroUpdateInterval)

Q = [5e-6*eye(3) zeros(3)       ;...
     zeros(3)       1e-16*eye(3)]; 

% Define [B(sigma)] matrix
B = @(mrp)((1-(mrp'*mrp)))*eye(3)+(2*tilde(mrp))+(2*(mrp*mrp'));

f = @(mrp,wbar)[0.25*B(mrp)*wbar; zeros(3,1)];

% Define F, G matrices
F = @(sBar,wBar)[0.5*(sBar*wBar' - wBar*sBar' - tilde(wBar) + sBar'*wBar*eye(3)) -0.25*B(sBar);...
    zeros(3) zeros(3)];
G = @(sBar)[-0.25*B(sBar) zeros(3);...
    zeros(3) eye(3)];

% Pull out MRP and bias estimates
b_bar = X(4:6);
s_bar = X(1:3);

% Compute estimated body rate
w_bar = wStar - b_bar;

% Propagate state estimate
fm = f(s_bar,w_bar);
Fm = F(s_bar,w_bar);
Gm = G(s_bar);

xDot = fm;
Pdot = Fm*P + P*Fm' + Gm*Q*Gm';

% Integrate state
X = X + xDot*gyroUpdateInterval;
P = P + Pdot*gyroUpdateInterval;

% Check for MRP set switch
MRP = X(1:3);
s = norm(X(1:3));
if s > 1
    X = [-MRP/(s^2); b_bar];
    Sm = 2*s^(-4)*(MRP*MRP') - s^(-2)*eye(3);
    P = [Sm*P(1:3,1:3)*Sm' Sm*P(1:3,4:6);...
           P(4:6,1:3)*Sm'     P(4:6,4:6)];
end

end


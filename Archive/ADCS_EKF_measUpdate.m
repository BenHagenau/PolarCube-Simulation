%------------------------------
% EKF for PolarCube ADCS
% Evaluates the measurement update for the EKF
%------------------------------
% Franklin Hinckley
% 24 January 2017
%------------------------------
%
%------------------------------

function [X, P, y, PFR] = ADCS_EKF_measUpdate(X, P, sStar, R)

% Define measurement mapping
H = [eye(3) zeros(3)];

% Pre-fit residual
y = sStar - X(1:3);

% Check if magnitude over threshold, check other MRP set
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

% Post-fit residuals
PFR = y - H*(K*y);

end 

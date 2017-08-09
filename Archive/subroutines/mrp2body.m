%------------------------------
% Converts MRP rates to body angular rates
%------------------------------
% Franklin Hinckley
% 9 October 2015
%------------------------------
% Inputs:
%   mrp (3x1 double vector): Modified Rodrigues Parameters
%   mrpDot (3x1 double vector): Modified Rodrigues Parameter Rates
% Outputs:
%   omega (3x1 double vector): angular velocity [rad/s]
% Dependencies:
%   tilde
%------------------------------

function [omega] = mrp2body(mrp,mrpDot)

% Compute norm of the MRP
sig = norm(mrp);

% Define kinematic differential equation
B = (1-sig^2)*eye(3)+(2*tilde(mrp))+(2*(mrp*mrp'));

% Compute omega
omega = (4/((1+sig^2)^2))*B'*mrpDot;

end
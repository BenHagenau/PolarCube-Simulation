%------------------------------
% Converts a body angular velocity to Modified Rodrigues 
%   Parameter Rates
%------------------------------
% Franklin Hinckley
% 1 October 2015
%------------------------------
% Inputs:
%   mrp (3x1 double vector): Modified Rodrigues Parameters
%   omega (3x1 double vector): angular velocity [rad/s]
% Outputs:
%   mrpDot (3x1 double vector): Modified Rodrigues Parameter Rates
% Dependencies:
%   tildeMat
%------------------------------

function [mrpDot]=body2mrp(omega,mrp)

% Compute norm of the MRP
sig = norm(mrp);

% Define kinematic differential equation
mrpDot = 0.25*((1-sig^2)*eye(3)+(2*tilde(mrp))+(2*(mrp*mrp')))*omega;

end
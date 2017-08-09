%------------------------------
% Converts a given set of Euler Parameters to Modified Rodrigues
% Parameters
%------------------------------
% Franklin Hinckley
% 28 August 2015
%------------------------------
% Inputs:
%   beta (4x1 double vector): Euler Parameters
% Outputs:
%   sigma (3x1 double vector): Modified Rodrigues Parameters
%------------------------------

function [sigma] = ep2mrp(beta)
% Check direction of rotation and use correct MRP set
if beta(1) > 0
    % Main set
    sigma = beta(2:4)/(1+beta(1));
else
    % Shadow set
    sigma = -beta(2:4)/(1-beta(1));
end
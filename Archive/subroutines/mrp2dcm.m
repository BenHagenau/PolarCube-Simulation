%------------------------------
% Converts a set of MRPs to the corresponding DCM
%------------------------------
% Franklin Hinckley
% 28 August 2015
%------------------------------
% Inputs:
%   MRP (3x1 double vector): Modified Rodrigues Parameters
% Outputs:
%   DCM (3x3 double matrix): Direction Cosine Matrix
%------------------------------

function [DCM] = mrp2dcm(MRP)
%Compute the norm of the mrp
s=norm(MRP);

%Compute the DCM
DCM = eye(3) + ((8*tilde(MRP)^2 - 4*(1-s^2)*tilde(MRP)) / ((1 + s^2)^2));

end
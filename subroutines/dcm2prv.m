%------------------------------
% Converts a given DCM to the corresponding Principal Rotation Vector
%------------------------------
% Franklin Hinckley
% 31 August 2015
%------------------------------
% Inputs:
%   C (3x3 double matrix): Direction Cosine Matrix
%   flag (string): flag for checking if angle is zero, 'ntol' disables
%       check, 'tol' sets check and tol variable is used
%   tol (double): tolerance for checking zero rotation angle
% Outputs:
%   phi (double): Principal Rotation Angle [rad]
%   ehat (3x1 double vector): Unit Principal Rotation Axis
%------------------------------

function [phi,ehat]=dcm2prv(C,flag,tol)

% Determine rotation angle
phi = acos(0.5*(trace(C)-1));

% Check for zero rotation (to within tolerance)
if strcmp(flag,'tol') && phi < tol
    error('Principal Rotation Angle is zero to within tolerance')
else
%Determine rotation axis
    ehat = zeros(3,1);
    ehat(1) = (1/(2*sin(phi)))*(C(2,3)-C(3,2));
    ehat(2) = (1/(2*sin(phi)))*(C(3,1)-C(1,3));
    ehat(3) = (1/(2*sin(phi)))*(C(1,2)-C(2,1));
end

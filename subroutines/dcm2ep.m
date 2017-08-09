%------------------------------
% Converts a given DCM to the corresponding Euler Parameter
% Uses Sheppard's Method, will return the short rotation
%------------------------------
% Franklin Hinckley
% 28 August 2015
%------------------------------
% Inputs:
%   C (3 x 3 double matrix): Direction Cosine Matrix
% Outputs:
%   ep (4 x 1 double matrix): Euler Parameters
%       ep(1) is the scalar term
%------------------------------

function [ep]=dcm2ep(C)
% Compute squares
b0sq = 0.25*(1+trace(C));
b1sq = 0.25*(1+2*C(1,1)-trace(C));
b2sq = 0.25*(1+2*C(2,2)-trace(C));
b3sq = 0.25*(1+2*C(3,3)-trace(C));
betaSq = [b0sq b1sq b2sq b3sq];

% Determine largest and take square root
[~,ind] = max(betaSq);
bi = sqrt(betaSq(ind));

% Determine parameters
switch ind
    case 1
        b0 = bi;
        b1 = (C(2,3)-C(3,2))/(4*b0);
        b2 = (C(3,1)-C(1,3))/(4*b0);
        b3 = (C(1,2)-C(2,1))/(4*b0);
    case 2
        b1 = bi;
        b0 = (C(2,3)-C(3,2))/(4*b1);
        b2 = (C(1,2)+C(2,1))/(4*b1);
        b3 = (C(3,1)+C(1,3))/(4*b1);
    case 3
        b2 = bi;
        b0 = (C(3,1)-C(1,3))/(4*b2);
        b1 = (C(1,2)+C(2,1))/(4*b2);
        b3 = (C(2,3)+C(3,2))/(4*b2);
    case 4
        b3 = bi;
        b0 = (C(1,2)-C(2,1))/(4*b3);
        b1 = (C(3,1)+C(1,3))/(4*b3);
        b2 = (C(2,3)+C(3,2))/(4*b3);
end

%Form output
ep = [b0;b1;b2;b3];

% Check sign of the scalar term
if ep(1) < 0
    ep = -ep;
end

end

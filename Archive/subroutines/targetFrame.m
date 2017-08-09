%--------------------------------
% Determines the reference DCM (input to target)
%   Uses the Triad method
%--------------------------------
% Franklin Hinckley
% 13 January 2016
%--------------------------------
% Inputs:
%   priTarg (3 x 1 double vector): primary pointing target
%   secTarg (3 x 1 double vector): secondary pointing target
% Outputs:
%   FT (3 x 3 double matrix): DCM for input frame relative to target frame 
%--------------------------------

function [FT] = targetFrame(priTarg, secTarg)

% Set first vector to primary target
t1 = priTarg/norm(priTarg);

% Determine second vector as cross-product of target vectors
t2 = cross(priTarg,secTarg)/norm(cross(priTarg,secTarg));

% Complete right-handed frame
t3 = cross(t1,t2);

% Form DCM
FT = [t1 t2 t3];

end
function RotationMatrix = rotMatrix(axis,angle)
% function : rotation.m
% Creates a rotation matrix for a rotation of a spacified "angle" about a 
% specified "axis". 
% 
% Inputs: 
%   axis [int] = axis the rotation is being performed [x=1,y=2,z=3]
%   angle  = angle for the rotation in +/- degrees.
% Outputs:
%   RotationMatrix = rotation matrix
%
% ------------------------------------------------------------------------

% Change the angle from degrees to radians
angle = deg2rad(angle);

% Determine about which axis the rotation is performed
switch axis
    case 1
        RotationMatrix = [1 0 0; ...
                        0 cos(angle) sin(angle); ...
                        0 -sin(angle) cos(angle)];
    case 2
        RotationMatrix = [cos(angle) 0 -sin(angle);...
                          0 1 0;...
                          sin(angle) 0 cos(angle)];
        
    case 3
        RotationMatrix = [cos(angle) sin(angle) 0; ...
                        -sin(angle) cos(angle) 0; ...
                        0 0 1];
    otherwise
        error('Please input a correct value for the axis.')
end

end
%------------------------------
% Compute the angle between a vehicle axis and a ground target
%------------------------------
% Franklin Hinckley
% 11 March 2017
%------------------------------
% 
%------------------------------

function [ang] = groundPointAngleCheck(vehiclePos,DCM,targetPos)

% Compute vector from vehicle location to ground location
relVec = targetPos - vehiclePos;
relUnitVec = relVec/nor(relVec);

% Extract vehicle axes
v_x = DCM(1,:);
v_y = DCM(2,:);
v_z = DCM(3,:);

% Compute angles
ang(1) = acos(dot(relUnitVec,v_x));
ang(2) = acos(dot(relUnitVec,-v_x));

ang(3) = acos(dot(relUnitVec,v_y));
ang(4) = acos(dot(relUnitVec,-v_y));

ang(5) = acos(dot(relUnitVec,v_z));
ang(6) = acos(dot(relUnitVec,-v_z));
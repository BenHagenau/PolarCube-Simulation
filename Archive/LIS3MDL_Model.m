%--------------------------------
% Simulate the response of the LIS3MDL MEMS magnetometer
%--------------------------------
% Franklin Hinckley
% 6 January 2016
%--------------------------------
% Inputs:
%   Bb (3 x 1 double vector): true magnetic field components in
%       body coordinates [mG]
%   SV (struct): structure of satellite properties
%       .mag1.bias (3 x 1 double vector): field direction bias of the 
%           magnetometer [rad/s]
%       .mag1.noise (3 x 1 double vector): standard deviation of the 
%           direction measurement noise 
%       .mag2.bias (3 x 1 double vector): field direction bias of the 
%           magnetometer [rad/s]
%       .mag2.noise (3 x 1 double vector): standard deviation of the 
%           direction measurement noise 
% Outputs:
%   bSim1 (3 x 1 double vector): simulated magnetic field components
%       in body coordinates from the first magnetometer [mG]
%   bSim2 (3 x 1 double vector): simulated magnetic field components
%       in body coordinates from the first magnetometer [mG]
%--------------------------------

function [bSim1, bSim2] = LIS3MDL_Model(Bb,SV)

% Simulate measurements
bSim1 = Bb;% + SV.mag1.bias + SV.mag1.noise.*randn(3,1);
bSim2 = Bb;% + SV.mag2.bias + SV.mag2.noise.*randn(3,1);

end

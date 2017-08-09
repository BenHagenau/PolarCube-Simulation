%--------------------------------
% Simulate the response of the A3G4250D MEMS gyroscope
%--------------------------------
% Franklin Hinckley
% 6 January 2016
%--------------------------------
% Inputs:
%   omega (3 x 1 double vector): true angular velocity [rad/s]
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

function [bSim1, bSim2, bTrue1, bTrue2] = A3G4250_Model(omega,SV,t)

% Compute true bias
bTrue1 = SV.gyro1.bias;% + SV.gyro1.drift*t;
bTrue2 = SV.gyro2.bias;% + SV.gyro2.drift*t;

% Simulate measurements
bSim1 = omega + SV.gyro1.noise.*randn(3,1) + bTrue1;
bSim2 = omega + SV.gyro2.noise.*randn(3,1) + bTrue2;

% % Load data file
% global gyroData
% 
% % Get indices
% % ind1 = ceil(93767*rand(1));
% % ind2 = ceil(93767*rand(1));
% ind1 = ceil(93757*rand(1));
% ind2 = ceil(93757*rand(1));
% 
% % Simulation measurements from the data file
% % bSim1 = omega + [gyroData(ind1,1:2)';gyroData(ind1,1)];
% % bSim2 = omega + [gyroData(ind2,1:2)';gyroData(ind2,1)];
% data1 = gyroData(ind1:ind1+9,:);
% data2 = gyroData(ind2:ind2+9,:);
% bSim1 = omega + mean(data1)';
% bSim2 = omega + mean(data2)';

end

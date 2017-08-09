%------------------------------
% Generates a series of ECI magnetic field data points
%------------------------------
% Franklin Hinckley
% 14 April 2016
%------------------------------
% 
%------------------------------


%% Clean up workspace
clearvars
close all
clc

%% Set up simulation
% Mission launch date
year = 2016; month = 01; day = 18;
hour = 4; minute = 45; second = 0;

% Propagator start time
SIM.prop.tStart = [year month day hour minute second];

% Propagator time step [sec] (this is effectively the GPS update rate)
SIM.prop.tStep = 1;

% Propagator duration [sec]
SIM.prop.tProp = 86400;  

% Location of TLE file
SIM.prop.TLEfname = 'polarcube.txt';

% Initial Julian Date
SIM.JD0 = juliandate(year, month,day, hour, minute, second);  

%% Run propagator
SIM = orbitPropagator(SIM);

%% Get magnetic field coefficients
global gh
[gh,~] = GetIGRF11_Coefficients(1);

%% Main loop
numPoints = floor(SIM.prop.tProp/SIM.prop.tStep);
B_ECI = zeros(3,numPoints);
for ii = 1:numPoints
    
    %% Determine current Julian Date
    JD = SIM.JD0 + (ii/86400*SIM.prop.tStep);
    
    %% Evaluate magnetic field model
    B_ECI(:,ii) = BfieldECI(SIM.prop.R(1,ii), ...
        SIM.prop.R(2,ii), SIM.prop.R(3,ii), JD);

end

    
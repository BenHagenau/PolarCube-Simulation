%------------------------------
% Orbit propagator main function for ADCS simulation
%------------------------------
% Franklin A. Hinckley
% Initial Version: 15 December 2015
% Latest Revision: 1.0
%   15 December 2015
%------------------------------

function [SIM] = orbitPropagator(SIM)

%% Pull out propagator settings struct
prop = SIM.prop;

%% Get TLE
%Parse TLE
TLE = readTLE(prop.TLEfname);

disp('Orbit Data Loaded')
fprintf('\n')

%% Orbit propagation
% Convert simulation start time to decimal day
tStart = dateEpoch(prop.tStart);

% Convert simulation time step to decimal day
tStep = prop.tStep/86400; 

% Set simulation stop time
tStop = tStart + (prop.tProp/86400);

% Run propagator
[R,V] = sgp4prop(TLE, tStart, tStep, tStop);
disp('Propagation Complete')
fprintf('\n')

% Assign output
SIM.prop.R = R;
SIM.prop.V = V;

end
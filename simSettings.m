%------------------------------
% Define simulation settings for PolarCube ADCS simulation
%------------------------------
% Franklin A. Hinckley
% Initial Version: 15 December 2015
% Latest Revision: 1.0
%   15 December 2015
%------------------------------

function [SIM] = simSettings()

%% Orbit propagator settings
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

%% Julian Date
% Initial Julian Date
SIM.JD0 = juliandate(year, month,day, hour, minute, second);  

%% Pointing targets
% COM mode
SIM.COMrngLim = 3000;         % Range limit to start COM operations, km
%SIM.COMrngLim = 0;
SIM.latCOM = 40.007648;         % Boulder GND station latitude (COSGC)
SIM.lonCOM = -105.261732;       % Boulder GND station longitude (COSGC)

% Science mode
SIM.SCIrngLim = 2000;     % Range limit to start SCI operations, km
SIM.latSCI = 90;            % North pole latitude
SIM.lonSCI = 0;              % North pole longitude

%% Mode controls
% Time for detumble
SIM.detTime = 15*60; % Fifteen minutes

%% External torques
% Define known external torques
Lx = @(t)0;
Ly = @(t)0;
Lz = @(t)0;

SIM.L = {Lx;Ly;Lz};

% Define unknown external torques
dLx = @(t)0;
dLy = @(t)0;
dLz = @(t)0;

SIM.dL = {dLx;dLy;dLz};

%% Simulation setup
% Simulation time [sec]
SIM.simTime = 60*60;

% Simulation frequency [Hz], this is the update rate for the dynamics
SIM.simFreq = 100;

% Control update frequency [Hz]
%SIM.controlFreq = 4;
SIM.controlFreq = 16;

%
SIM.measFreq = 25;

%% Initial conditions
% Define initial conditions
SIM.MRP0 = [0.5;0.25;-0.25];
%SIM.MRP0 = [0;0;0];
SIM.omega0 = [10;20;-10]*(pi/180)/4;
%SIM.omega0 = [0;0;0];
%SIM.omega0 = -[-0.00057703;0.0009007;0.012223];

%% Simulation Flags
% Enable dynamics
SIM.RWdyn  = 1;

% Enable control 
SIM.RWmode  = 1;
SIM.integralControl = 0;
SIM.mirrorControl = 0;

%% Plot control (1 enables the plot)
% Mode-specific plots
SIM.plots.detumble = 1; % MRPs, rate
SIM.plots.bus      = 1; % MRPs vs ref, rate, log-scale attitude/knowledge error
SIM.plots.science  = 1; % MRPs vs ref, rate, log-scale attitude/knowledge error
SIM.plots.COM      = 1; % MRPs vs ref, rate, log-scale attitude/knowledge error

% Full simulation plots
SIM.plots.MRP      = 1; % True vs reference MRPs
SIM.plots.MRPest   = 1; % Estimated vs true MRPs
SIM.plots.detErr   = 1; % Log-scale attitude determination error
SIM.plots.omega    = 1; % True angular rate
SIM.plots.omegaEst = 1; % Estimated angular rate
SIM.plots.omegaErr = 1; % Error in velocity estimate
SIM.plots.gyroBias = 1; % Estimated gyroscope bias
SIM.plots.H        = 0; % System angular momentum
SIM.plots.Hrel     = 0; % Relative variation in system angular momentum
SIM.plots.HSV      = 0; % Spacecraft bus angular momentum
SIM.plots.HRW      = 0; % Reaction wheel angular momentum
SIM.plots.HM       = 0; % Mirror angular momentum
SIM.plots.T        = 0; % System kinetic energy
SIM.plots.P        = 0; % Power computed by power eqn and finite diff of T
SIM.plots.us       = 1; % Control torques
SIM.plots.Omega    = 1; % Wheel speeds
SIM.plots.OmegaDot = 0; % Wheel accelerations
SIM.plots.mu       = 0; % Mirror angular position
SIM.plots.muDot    = 0; % Mirror angular rate
SIM.plots.Phi      = 0; % Pointing error
SIM.plots.pointErr = 1; % Log-scale pointing error
SIM.plots.cov      = 1; % Attitude errors with 3-sigma covariance envelope

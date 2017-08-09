%------------------------------
% Control Function for PolarCube ADCS Simulation
%------------------------------
% Franklin A. Hinckley
% Initial Version: 15 December 2015
% Latest Revision: 1.0
%   15 December 2015
%------------------------------

function init

%% Clean up workspace and include subfolders
clearvars
close all
clc

addpath subroutines
addpath propagator
code_source = 2; % 1 for MATLAB, 2 for C code
%% Print header
% fprintf(['\n',...
%     'PolarCube ADCS Mission Simulation\n\n', ...
%     'An open source GNSS SDR software project initiated by:\n\n', ...
%     '              Danish GPS Center/Aalborg University\n\n', ...
%     'The code was improved by GNSS Laboratory/University of Colorado.\n\n',...
%     'The software receiver softGNSS comes with ABSOLUTELY NO WARRANTY;\n',...
%     'for details please read license details in the file license.txt. This\n',...
%     'is free software, and  you  are  welcome  to  redistribute  it under\n',...
%     'the terms described in the license.\n\n']);
% fprintf('                   -------------------------------\n\n');

%% Get satellite properties
SV = satProps();

%% Get simulation settings
SIM = simSettings();

%% Get magnetic field coefficients
global gh
[gh,~] = GetIGRF11_Coefficients(1);

%% Determine satellite position/velocity
%prop = input('Enter "1" to initiate GNSS processing or "0" to exit : ');
SIM = orbitPropagator(SIM);

%% Run simulation
[simResults,SV] = ADCS_Sim(SIM,SV,code_source);
disp('Simulation Complete')
fprintf('\n')

%% Check simulation results
[checks] = ADCS_Checks(simResults,SV,(1/SIM.simFreq));
simResults.checks = checks;

%% Save results
save('simulationResults', 'simResults', 'SV', 'SIM');
disp('Results Saved')
fprintf('\n')

%% Plot results
ADCS_Plot(simResults,SIM,SV)

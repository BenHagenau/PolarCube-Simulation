%--------------------------------
% Checks the current satellite state and returns the correct mode
%--------------------------------
% Franklin Hinckley
% 23 October 2015
%--------------------------------
% Inputs:
%
% Outputs:
%
%--------------------------------

function [mode] = modeCheck(R)

% COM mode when in range of facility
INrange_COM = 3000;         % Range limit to start COM operations, km
latCOM = 40.007648;         % Boulder GND station latitude (COSGC)
lonCOM = -105.261732;       % Boulder GND station longitude (COSGC)

% SCIENCE mode, where the mirror is pointed nadir
INrange_SCI = 2000;     % Range limit to start SCI operations, km
latSCI = 90;            % North pole latitude
lonSCI= 0;              % North pole longitude
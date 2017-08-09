%--------------------------------
% Sets the parameters for the current mode
%--------------------------------
% Franklin Hinckley
% 12 January 2016
%--------------------------------
% Inputs:
%   SV (struct): satellite properties
%   SIM (struct): simulation properties
% Outputs:
%   SV (struct): satellite properties
%   SIM (struct): simulation properties
%   mode (string): defines current mode
%--------------------------------

function [SV, SIM] = modeSet(SV, SIM, mode)

%% Check mode
switch mode
    
%% Detumble
    case 'det'
    % Reference MRP (not actually used in detumble)    
    SIM.MRPRef = {@(t)0;@(t)0;@(t)0};

    % Reference Rate
    SIM.omegaRef = @(t,mrp)[0;0;0];

    % Gains
    SV.K = 0;
    SV.P = diag([0.0075; 0.0075; 0.0075]);
    SV.Ki = diag([0.0025;0.0025;0.0025]);

    % Flags
    SIM.RWdyn  = 1;
    SIM.mirrorDyn = 1;

    SIM.RWmode  = 1;
    SIM.integralControl = 0;

%% Bus Priority
    case 'bus'
    % Reference MRP
    sX = @(t)0.25*cos(0.01*t);
    sY = @(t)0.1*sin(0.01*t);
    sZ = @(t)0.1*sin(0.01*t);
    SIM.MRPRef = {sX;sY;sZ};

    % Reference Rate
    SIM.omegaRef = @(t,mrp)mrp2body(mrp,...
        [-0.0025*sin(0.01*t) 0.001*cos(0.01*t) 0.001*cos(0.01*t)]');

    % Gains
    SV.K = 0.005;
    SV.P = diag([0.0075; 0.0075; 0.0075]);
    SV.Ki = diag([0.0025;0.0025;0.0025]);

    % Flags
    SIM.RWdyn  = 1;
    SIM.mirrorDyn = 1;

    SIM.RWmode  = 1;
    SIM.integralControl = 0;

%% S-Band COM
    case 'com'
    % Reference MRP
    sX = @(t)0.25*cos(0.01*t);
    sY = @(t)0.1*sin(0.01*t);
    sZ = @(t)0.1*sin(0.01*t);
    SIM.MRPRef = {sX;sY;sZ};

    % Reference Rate
    SIM.omegaRef = @(t,mrp)mrp2body(mrp,...
        [-0.0025*sin(0.01*t) 0.001*cos(0.01*t) 0.001*cos(0.01*t)]');

    % Gains
    SV.K = 0.005;
    SV.P = diag([0.0075; 0.0075; 0.0075]);
    SV.Ki = diag([0.0025;0.0025;0.0025]);

    % Flags
    SIM.RWdyn  = 1;
    SIM.mirrorDyn = 1;

    SIM.RWmode  = 1;
    SIM.integralControl = 0;

%% Science
    case 'sci'
    % Reference MRP
    sX = @(t)0.25*cos(0.01*t);
    sY = @(t)0.1*sin(0.01*t);
    sZ = @(t)0.1*sin(0.01*t);
    SIM.MRPRef = {sX;sY;sZ};

    % Reference Rate
    SIM.omegaRef = @(t,mrp)mrp2body(mrp,...
        [-0.0025*sin(0.01*t) 0.001*cos(0.01*t) 0.001*cos(0.01*t)]');

    % Gains
    SV.K = 0.005;
    SV.P = diag([0.0075; 0.0075; 0.0075]);
    SV.Ki = diag([0.0025;0.0025;0.0025]);

    % Flags
    SIM.RWdyn  = 1;
    SIM.mirrorDyn = 1;

    SIM.RWmode  = 1;
    SIM.integralControl = 1;

end 
%------------------------------
% Define satellite properties for PolarCube ADCS Simulation
%------------------------------
% Franklin A. Hinckley
% Initial Version: 15 December 2015
% Latest Revision: 1.0
%   15 December 2015
%------------------------------

function [SV] = satProps()

%% Satellite properties
% Inertia [kg*m^2]
SV.I = [0.05220503  -0.00026736   0.000475987;...
       -0.00026736   0.052387143 -0.000342469;...
        0.000475987 -0.000342469  0.004468478];
         
%% Reaction wheel parameters
% Inertias [kg*m^2]
SV.Js = [1.1249e-5; 1.1249e-5; 1.1249e-5];
SV.Jt1 = [5.6735e-6; 5.6735e-6; 5.6735e-6];
SV.Jt2 = [5.6735e-6; 5.6735e-6; 5.6735e-6];

% Wheel frame axes
SV.Gs = [-1 0 0;...
          0 1 0;...
          0 0 1];
SV.Gt1 = [0 0 1;...
         -1 0 0;...
          0 1 0];
SV.Gt2 = [0 1 0;...
          0 0 1;...
         -1 0 0];
      
% Nominal wheel speed [rad/s]      
SV.OmegaNom = [250; 250; 250]/60/(2*pi);

% Torque limit [Nm]
SV.RWtorqueLim = 0.27/1000;

%% Torque rod parameters
% Axes
SV.Gt = [-1 0 0;...
          0 1 0;...
          0 0 1];
      
%% Mirror parameters
% Mirror inertia [kg*m^2]
SV.Im = [305.19 -0.53 0.18;...
         -0.53 306.11 96.05;...
         0.61 94.68 332.18]/1000/100/100/3;
     
% Initial mirror rotation angle
SV.mu0 = 0;
SV.mu = 0;
SV.contFilt = 0;
SV.Lprev = zeros(3,2);
% Initial mirror frame axes
SV.ms0 = [1/sqrt(2) 1/sqrt(2) 0]';
SV.ml0 = [-sind(45)*sind(45) sind(45)*cosd(45) cosd(45)]';
SV.mr0 = [sind(45)*cosd(45) -sind(45)*sind(45) cosd(45)]';
SV.BM = [SV.ms0 SV.ml0 SV.mr0];

% Mirror rate and acceleration
% SV.muDot = @(t)1.3*2*pi + (0.05 * 1.3*2*pi)*sin(2*pi*0.1*t);
% SV.muDoubDot = @(t)2*pi*0.1*(0.05 * 1.3*2*pi)*cos(2*pi*0.1*t);
%SV.muDot = @(t)1.3*2*pi;
SV.muDot = @(t)0;
SV.muDoubDot = @(t)0;

%% Control gains
% SV.K = 0.005;
% SV.P = diag([0.0075; 0.0075; 0.0075]);
% SV.Ki = diag([0.0025;0.0025;0.0025]);
SV.K = 0.005;
SV.P = diag([0.0075; 0.0075; 0.0075]);
SV.Ki = zeros(3,3);%diag([0.0025;0.0025;0.0025]);
SV.Km = 0.0025;
SV.Kb = 2.5e-3;

%% Define initial conditions for reaction wheels
SV.Omega = [1;1;1]*250/60/(2*pi); %[rad/s] (250 rpm)
SV.OmegaDot = [0;0;0];

%% Gyroscope parameters
% Sampling frequency [Hz]
SV.gyro1.sampFreq = 100;
SV.gyro2.sampFreq = 100; % not currently used in sim, assumes same as gyro1

% Rate measurement bias [rad/s]
SV.gyro1.bias = [-0.0137; 0.2418; -0.0202]*(pi/180);
SV.gyro2.bias = [-0.2418; 0.0202; 0.0137]*(pi/180);

% Rate measurement noise [rad/s]
SV.gyro1.noise = [0.1328; 0.1203; 0.1385]*(pi/180);
SV.gyro2.noise = [0.1328; 0.1203; 0.1385]*(pi/180);

% Rate measurement bias drift [rad/s/s]
SV.gyro1.drift = -[0.0001/60; 0.0002/60; 0.0002/60]*(pi/180);
SV.gyro2.drift = -[0.0001/60; 0.0002/60; 0.0002/60]*(pi/180);

% Load gyro data file
% global gyroData
% gyroDataRaw = load('gyroLog10Hz.txt');
% Wrange = 245;    
% Crange = 2^15;  
% C2W = Wrange/Crange; 
% gyroData = gyroDataRaw*C2W;
% SV.gyro1.bias = [mean(gyroData(:,1:2))';mean(gyroData(:,1))];
% SV.gyro2.bias = [mean(gyroData(:,1:2))';mean(gyroData(:,1))];

%% Magnetometer parameters  (need updated to actual values)
% Sampling frequency [Hz]
SV.mag1.sampFreq = 100;
SV.mag2.sampFreq = 100; % not currently used in sim, assumes same as mag1

% Field measurement bias [nT]
SV.mag1.bias = [1000; 1000; 1000];%*100;
SV.mag2.bias = [1000; 1000; 1000];%*100;

% Field measurement noise [nT]
SV.mag1.noise = [3.2; 3.2; 3.2];%*100;
SV.mag2.noise = [3.2; 3.2; 3.2];%*100;

%% Star tracker parameters
% Update rate [Hz]
SV.starTrk.sampFreq = 0.2;

% Measurement noise [rad]
SV.starTrk.noise = 100*4.84814e-6; %100 arc-sec

%% Control Filter Parameters
% Loop bandwidth and damping ratio
SV.gyroFilt.BW = 0.1;  % [Hz]
SV.gyroFilt.DR = 0.7; % sqrt(2) for optimally flat response

SV.deadband = 0.01/1000;


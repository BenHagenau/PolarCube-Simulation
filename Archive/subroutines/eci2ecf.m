function [ECF] = eci2ecf(X, Y, Z, JD)

%This function takes the input of ECI Coordinates (X,Y,Z) as well as time
%elapsed since Jan 1, 2006 00h 00m 00s to calculate the ECG cartesian
%coordinates (x,y,z)
t = (JD - 2453736.5) * 86400;

%greenwich hour angle
g=100.43878;    %degrees
%degrees-->radians conversion
grad=g*pi/180;  %radians

%Earth's rotation rate
w=7.2921158553e-5;   %rad/s

%greenwich hour angle at time t
gt=grad+w*t;

%transformation matrix
xform=[cos(gt) sin(gt) 0; -sin(gt) cos(gt) 0; 0 0 1];

[ECF]=xform*[X; Y ;Z];

%------------------------------
% Determines the body-frame magnetic field at the current time
%------------------------------
% Franklin A. Hinckley
% Initial Version: 15 December 2015
% Latest Revision: 1.0
%   15 December 2015
%------------------------------
% Inputs:
%   MRP (3 x 1 double array): current MRP for body relative to inertial
%   X (double): X-position in ECI [km]
%   Y (double): Y-position in ECI [km]
%   Z (double): Z-position in ECI [km]
%   JD (double): current Julian date [days]
% Outputs:
%   Bbody (3 x 1 double array): magnetic field components in body
%       coordinates [nT]
%------------------------------

function [Bbody] = BfieldComp(MRP, X, Y, Z, JD)

% Determine the magnetic field WRT to earth centered inertial frame
r = norm([X Y Z]);
alt =  r - 6378.1363; % Altitude: radius to spacecraft minus radius of earth [km]

% Earth centered frame:
[ECF] = eci2ecf(X, Y, Z, JD);
x = ECF(1);
y = ECF(2);
z = ECF(3);
nlat = asin(z/r); % Latitude
elong = atan2(y,x); % Longitude
year = JD/365.25 - 4800; % Time 

% Determine magnetic field in earth centered frame
[B] = igrf11syn(year,alt,nlat*180/pi,elong*180/pi);

% Convert to earth centered inertial frame
B_eci(1,:) = sin(elong) * B(2) - sin(nlat) * cos(elong) * B(1) + cos(nlat) * cos(elong) * B(3);
B_eci(2,:) = cos(elong) * B(2) - sin(nlat) * sin(elong) * B(1) + cos(nlat) * sin(elong) * B(3);
B_eci(3,:) = cos(nlat) * B(1) + sin(nlat) * B(3);

% Convert magnetic field in ECI to body frame
mag_sigma = sqrt(MRP'*MRP);
sigma_tilde = tilde(MRP);
C = eye(3)+(8*sigma_tilde^2-4*(1-mag_sigma^2)*sigma_tilde)/((1+mag_sigma^2)^2);
Bbody = C * B_eci;


function [B_ECI] = BfieldECI(X, Y, Z, JD)

% Get spacecraft altitude
r = norm([X Y Z]);
alt =  r - 6378.1363; 

% Compute Earth-centered frame (approximation as ignores polar motion)
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
B_ECI(1,:) = sin(elong) * B(2) - sin(nlat) * cos(elong) * B(1) + cos(nlat) * cos(elong) * B(3);
B_ECI(2,:) = cos(elong) * B(2) - sin(nlat) * sin(elong) * B(1) + cos(nlat) * sin(elong) * B(3);
B_ECI(3,:) = cos(nlat) * B(1) + sin(nlat) * B(3);
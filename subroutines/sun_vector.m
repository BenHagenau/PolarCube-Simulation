function sun_vect = sun_vector(JD)
%% Ephemeris calculation
D = JD - 2451545;                           
g = 357.529 + 0.98560028 * D;                   % deg
q = 280.459 + 0.98564736 * D;                   % deg
L = q + 1.915 * sind(g) + 0.020 * sind(2 * g);  % deg
% R = 1.00014 - 0.01671 * cosd(g) - 0.00014 * cosd(2 * g);
e = 23.439 - 3.6e-7 * D;                        % deg

RA = atan2(cosd(e) * sind(L),cosd(L));          % rad
DEC = asin(sind(e) * sind(L));                  % rad

X_sun = sin(pi / 2 - DEC) * cos(RA);
Y_sun = sin(pi / 2 - DEC) * sin(RA);
Z_sun = cos(pi / 2 - DEC);
sun_vect = [X_sun; Y_sun; Z_sun];
end
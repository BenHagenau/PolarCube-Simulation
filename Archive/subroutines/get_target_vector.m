%------------------------------
% Calculates the vector from the center of the earth
%------------------------------

%------------------------------

%------------------------------

function target_vector = get_target_vector(lat,long,JD)
%This function calculates the target vector from earth center to a target
%specified by lat, long, and Julian date.

%Inputs: 
%   Latitude (Degrees)
%   Longitude (Degrees)
%   Julian Date

%Output:
%   Target vector: containing the x, y, z coordinates in earth centered
%   frame

%%
lat = lat * pi / 180;               % Convert from deg. to rad.
long = long * pi /180;              % Convert from deg. to rad.
time = (JD - 245376.5) * 86400;     % Convert from JD to sec
earth_rotation = 7.2921158553e-5;   % [rad/s]

% Compute the rotation of the earth at JD
rotation = 100.43878 * pi / 180 + earth_rotation * time;

%Convert to earth centered frame:
X_ECF = 6378 * sin(pi / 2 - lat) * cos(long);
Y_ECF = 6378 * sin(pi / 2 - lat) * sin(long);
Z_ECF = 6378 * cos(pi / 2 - lat);
target_vector = [cos(rotation) sin(rotation) 0;
                 -sin(rotation) cos(rotation) 0;
                                           0 0 1]' * [X_ECF; 
                                                     Y_ECF; 
                                                     Z_ECF];
%------------------------------
% Parses out a TLE from a text file
%------------------------------
% Franklin Hinckley
% 15 January 2016
%------------------------------
% Inputs: 
%   TLEInFile (string): file name of TLE, file should be a text file 
% Outputs:
%   TLE (struct): parsed TLE, each entry in the array is a 
%       struct of the fields for a single TLE
% Usage: 
%   TLE is assumed to be in the 3LE format and have all fields 
%       whitespace-delimited (the standard space-track.org 3LE format)
%   The TLE is assumed to be complete
%------------------------------

function [TLE] = readTLE(TLEInFile)
% Open input file
inf = fopen(TLEInFile,'rt');

% Main body
%Initialize counter
i = 1;    

%% Parse fields
% Read lines
try
    line0 = fgetl(inf);
    line1 = fgetl(inf);
    line2 = fgetl(inf);
catch
    error('Full TLE Not Defined')
end

% Parse line 0
TLEtmp.satName = line0;

% Break line 1 
line1entries = textscan(line1(3:end),'%s');
line1entries = line1entries{:};

% Assign line 1 fields
TLEtmp.satNumber      = line1entries{1};
TLEtmp.class          = TLEtmp.satNumber(end);
TLEtmp.satNumber      = TLEtmp.satNumber(1:end-1);
TLEtmp.intDesig       = line1entries{2};
TLEtmp.epoch          = line1entries{3};
TLEtmp.meanMotionDev1 = line1entries{4};
TLEtmp.meanMotionDev2 = line1entries{5};
TLEtmp.BSTAR          = line1entries{6};
TLEtmp.elSet          = line1entries{8};
TLEtmp.checksum1      = TLEtmp.elSet(end);
TLEtmp.elSet          = TLEtmp.elSet(1:end-1);

% Break line 2 fields
line2entries = textscan(line2(9:end),'%s');
line2entries = line2entries{:};

% Assign line 2 fields
TLEtmp.incl        = line2entries{1};
TLEtmp.RAAN        = line2entries{2};
TLEtmp.ecc         = line2entries{3};
TLEtmp.argPerigee  = line2entries{4};
TLEtmp.meanAnomaly = line2entries{5};
if length(line2entries) == 7
    TLEtmp.meanMotion  = line2entries{6};
    TLEtmp.revNumber   = line2entries{7};
else
    tmp = line2entries{6};
    TLEtmp.meanMotion  = tmp(1:11);
    TLEtmp.revNumber   = tmp(12:end);
end
    
TLEtmp.checksum2   = TLEtmp.revNumber(end);
TLEtmp.revNumber   = TLEtmp.revNumber(1:end-1);

%% Assign string fields
TLE.satName  = TLEtmp.satName(find(isstrprop(TLEtmp.satName,'alpha'),1):end);
TLE.class    = TLEtmp.class;
TLE.intDesig = TLEtmp.intDesig;

%% Fix scale factors and convert to decimal
%Set format long so str<->num conversions retain full accuracy
format long

% First time derivative of mean motion
TLE.meanMotionDev1 = str2double(TLEtmp.meanMotionDev1)*2;

% Second time deriviative of mean motion
sign = 1;
if strcmp(TLEtmp.meanMotionDev2(6),'-')
    sign = -1;
end
TLE.meanMotionDev2 = str2double(['0.' TLEtmp.meanMotionDev2(1:5)])*6*(10^(str2double(TLEtmp.meanMotionDev2(7))*sign));

% BSTAR drag term
sign = 1;
if strcmp(TLEtmp.BSTAR(6),'-')
    sign = -1;
end
TLE(i).BSTAR = str2double(['0.' TLEtmp.BSTAR(1:5)])*6*(10^(str2double(TLEtmp.BSTAR(7))*sign));

%% Convert remaining values to decimal 
TLE.satNumber   = str2double(TLEtmp.satNumber);
TLE.epoch       = str2double(TLEtmp.epoch);
TLE.checksum1   = str2double(TLEtmp.checksum1);
TLE.elSet       = str2double(TLEtmp.elSet);

TLE.satNumber   = str2double(TLEtmp.satNumber);
TLE.incl        = str2double(TLEtmp.incl);
TLE.RAAN        = str2double(TLEtmp.RAAN);
TLE.ecc         = str2double(['0.' TLEtmp.ecc]);   %fix assumed decimal place
TLE.argPerigee  = str2double(TLEtmp.argPerigee);
TLE.meanAnomaly = str2double(TLEtmp.meanAnomaly);
TLE.meanMotion  = str2double(TLEtmp.meanMotion)*(2*pi)/86400;    %convert units to rad/s
TLE.revNumber   = str2double(TLEtmp.revNumber);
TLE.checksum2   = str2double(TLEtmp.checksum2);

%% Add semimajor axis to struct
mu = 398600.4418;
TLE.semiMaj = (mu/(TLE(i).meanMotion^2))^(1/3);

%Set format short for typical command window use
format short

%% Verify fields
fields = fieldnames(TLE);
for ii = 1:length(fields)
    if isnan(TLE.(fields{ii}))
        error('TLE Field Contains NaN')
    end
end

end 
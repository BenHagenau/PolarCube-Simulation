%------------------------------
% Smooths the attitude data
%------------------------------
% Franklin Hinckley
% 15 January 2016
%------------------------------
% Inputs: 
%   TLEInFile (string): file name of TLE, file should be a text file 
% Outputs:
%   TLE (array of structs): parsed TLE, each entry in the array is a 
%       struct of the fields for a single TLE
% Usage: 
%   TLEs are assumed to be in the 3LE format and have all fields 
%       whitespace-delimited (the standard space-track.org 3LE format)
%   Each TLE is assumed to be complete
%------------------------------

function [TLE] = TLEMultRead(TLEInFile)
% Open input file
inf = fopen(TLEInFile,'rt');

% Main body
%Initialize counter
i = 1;    

while ~feof(inf)
    %% Parse fields
    %Parse the TLE into lines (3 line format, must have correct spaces)
    line1 = fgetl(inf);
    line2 = fgetl(inf);
    line3 = fgetl(inf);

    %Parse line 1
    TLEtmp.satName = line1;

    %Parse line 2
    line2entries = textscan(line2(3:end),'%s');
    line2entries = line2entries{:};

    TLEtmp.satNumber      = line2entries{1};
    TLEtmp.class          = TLEtmp.satNumber(end);
    TLEtmp.satNumber      = TLEtmp.satNumber(1:end-1);
    TLEtmp.intDesig       = line2entries{2};
    TLEtmp.epoch          = line2entries{3};
    TLEtmp.meanMotionDev1 = line2entries{4};
    TLEtmp.meanMotionDev2 = line2entries{5};
    TLEtmp.BSTAR          = line2entries{6};
    TLEtmp.elSet          = line2entries{8};
    TLEtmp.checksum1      = TLEtmp.elSet(end);
    TLEtmp.elSet          = TLEtmp.elSet(1:end-1);

    %Parse line 3
    line3entries = textscan(line3(9:end),'%s');
    line3entries = line3entries{:};

    TLEtmp.incl        = line3entries{1};
    TLEtmp.RAAN        = line3entries{2};
    TLEtmp.ecc         = line3entries{3};
    TLEtmp.argPerigee  = line3entries{4};
    TLEtmp.meanAnomaly = line3entries{5};
    TLEtmp.meanMotion  = line3entries{6};
    TLEtmp.revNumber   = line3entries{7};
    TLEtmp.checksum2   = TLEtmp.revNumber(end);
    TLEtmp.revNumber   = TLEtmp.revNumber(1:end-1);
    
    %% Assign string fields
    TLE(i).satName  = TLEtmp.satName(find(isstrprop(TLEtmp.satName,'alpha'),1):end);
    TLE(i).class    = TLEtmp.class;
    TLE(i).intDesig = TLEtmp.intDesig;

    %% Fix scale factors and convert to decimal
    %Set format long so str<->num conversions retain full accuracy
    format long
    %First time derivative of mean motion
    TLE(i).meanMotionDev1=str2double(TLEtmp.meanMotionDev1)*2;
    %Second time deriviative of mean motion
    sign=1;
    if strcmp(TLEtmp.meanMotionDev2(6),'-')
        sign=-1;
    end
    TLE(i).meanMotionDev2=str2double(['0.' TLEtmp.meanMotionDev2(1:5)])*6*(10^(str2double(TLEtmp.meanMotionDev2(7))*sign));
    %BSTAR drag term
    sign=1;
    if strcmp(TLEtmp.BSTAR(6),'-')
        sign=-1;
    end
    TLE(i).BSTAR=str2double(['0.' TLEtmp.BSTAR(1:5)])*6*(10^(str2double(TLEtmp.BSTAR(7))*sign));

    %% Convert remaining values to decimal 
    TLE(i).satNumber=str2double(TLEtmp.satNumber);
    TLE(i).epoch=str2double(TLEtmp.epoch);
    TLE(i).checksum1=str2double(TLEtmp.checksum1);
    TLE(i).elSet=str2double(TLEtmp.elSet);

    TLE(i).satNumber=str2double(TLEtmp.satNumber);
    TLE(i).incl=str2double(TLEtmp.incl);
    TLE(i).RAAN=str2double(TLEtmp.RAAN);
    TLE(i).ecc=str2double(['0.' TLEtmp.ecc]);   %fix assumed decimal place
    TLE(i).argPerigee=str2double(TLEtmp.argPerigee);
    TLE(i).meanAnomaly=str2double(TLEtmp.meanAnomaly);
    TLE(i).meanMotion=str2double(TLEtmp.meanMotion)*(2*pi)/86400;    %convert units to rad/s
    TLE(i).revNumber=str2double(TLEtmp.revNumber);
    TLE(i).checksum2=str2double(TLEtmp.checksum2);
    
    %% Add semimajor axis to struct
    mu=398600.4418;
    TLE(i).semiMaj=(mu/(TLE(i).meanMotion^2))^(1/3);
    
    %Set format short for typical command window use
    format short
    %Increment counter
    i=i+1;
    
end %end main loop

end %end function
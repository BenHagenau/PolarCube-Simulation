%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Converts a MATLAB clock time to an epoch time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Franklin Hinckley
% Latest Revision: 18 December 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   clock (MATLAB date vector)
% Outputs:
%   epoch (double): epoch time [TLE epoch]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [epoch]=dateEpoch(clock)
%Parse clock
YYYY=clock(1);
MM=clock(2);
DD=clock(3);
H=clock(4);
M=clock(5);
S=clock(6);

%Form decimal day
sec=(3600*H)+(60*M)+S;
decDay=sec/86400;

%Form day number
switch MM
    case 1
        day=DD;
    case 2
        day=DD+31;
    case 3
        day=DD+59;
    case 4
        day=DD+90;
    case 5
        day=DD+120;
    case 6
        day=DD+151;
    case 7
        day=DD+181;
    case 8
        day=DD+212;
    case 9
        day=DD+243;
    case 10
        day=DD+273;
    case 11
        day=DD+304;
    case 12
        day=DD+334;
end

%Form year
year=(YYYY-2000)*1000;

%Form epoch time
epoch=year+day+decDay;

end
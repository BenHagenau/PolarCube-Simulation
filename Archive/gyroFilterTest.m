%------------------------------
% Gyroscope filtering testing
%------------------------------
% Franklin Hinckley
% 11 February 2016
%------------------------------

%% Clean up
close all
clc

%% Plot time-domain results
figure
plot(simResults.wStar1)

figure
plot(simResults.t,simResults.wStar2)

%% Design filter
% [b,a] = butter(2,0.1257);
% 
% %% Filter data
% wStar1filt = filter(b,a,simResults.wStar1(50000:100000));
% wStar2filt = filter(b,a,simResults.wStar2);
% 
% %% Plot filtered time-domain results
% figure
% plot(wStar1filt)
% 
% figure
% plot(simResults.t,wStar2filt)



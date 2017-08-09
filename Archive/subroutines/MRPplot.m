%--------------------------------
% Plots MRPs as continuous series without steps at the switches
% Saves plot to requested folder
%--------------------------------
% Franklin Hinckley
% 23 October 2015
%--------------------------------
% Inputs:
%   MRP (3 x n double matrix): simulation MRPs
%   MRPref (3 x n double matrix): reference MRPs
%   t (n x 1 double vector): time series
%   plotFol (string): folder for saved plot
%   plotName (string): file name for saved plot
%   psflag (double): plot save flag, setting this to 1 saves the plot
%--------------------------------

function MRPplot(MRP,MRPRef,t,titleStr,plotFol,plotName,psflag)

%% Find switches and insert NaNs for simulation MRPs
% Find switches
sw = [];
for i = 2:length(t)
    s1 = MRP(1,i)*MRP(1,i-1);
    s2 = MRP(2,i)*MRP(2,i-1);
    s3 = MRP(3,i)*MRP(3,i-1);
    if (s1<0) || (s2<0) || (s3<0)
        sw = [sw i];
    end
end

% Insert NaN
if ~isempty(sw)
    %First block
    MRP1=MRP(1,1:sw(1)-1);
    MRP2=MRP(2,1:sw(1)-1);
    MRP3=MRP(3,1:sw(1)-1);
    tvec=t(1:sw(1)-1);
    %Middle blocks
    for i=2:length(sw)
        MRP1=[MRP1 NaN MRP(1,sw(i-1):sw(i)-1)];
        MRP2=[MRP2 NaN MRP(2,sw(i-1):sw(i)-1)];
        MRP3=[MRP3 NaN MRP(3,sw(i-1):sw(i)-1)];
        tvec=[tvec NaN t(sw(i-1):sw(i)-1)];
    end
    %End block
    MRP1=[MRP1 NaN MRP(1,sw(end):end)];
    MRP2=[MRP2 NaN MRP(2,sw(end):end)];
    MRP3=[MRP3 NaN MRP(3,sw(end):end)];
    tvec=[tvec NaN t(sw(end):end)];
else
    tvec=t;
    MRP1=MRP(1,:);
    MRP2=MRP(2,:);
    MRP3=MRP(3,:);
end

%% Find switches and insert NaNs for reference
% Find switches
sw=[];
for i=2:length(t)
    s1=MRPRef(1,i)*MRPRef(1,i-1);
    s2=MRPRef(2,i)*MRPRef(2,i-1);
    s3=MRPRef(3,i)*MRPRef(3,i-1);
    if (s1<0) || (s2<0) || (s3<0)
        sw=[sw i];
    end
end

% Insert NaN
if ~isempty(sw)
    %First block
    MRP1Ref=MRPRef(1,1:sw(1)-1);
    MRP2Ref=MRPRef(2,1:sw(1)-1);
    MRP3Ref=MRPRef(3,1:sw(1)-1);
    tvecRef=t(1:sw(1)-1);
    %Middle blocks
    for i=2:length(sw)
        MRP1Ref=[MRP1Ref NaN MRPRef(1,sw(i-1):sw(i)-1)];
        MRP2Ref=[MRP2Ref NaN MRPRef(2,sw(i-1):sw(i)-1)];
        MRP3Ref=[MRP3Ref NaN MRPRef(3,sw(i-1):sw(i)-1)];
        tvecRef=[tvecRef NaN t(sw(i-1):sw(i)-1)];
    end
    %End block
    MRP1Ref=[MRP1Ref NaN MRPRef(1,sw(end):end)];
    MRP2Ref=[MRP2Ref NaN MRPRef(2,sw(end):end)];
    MRP3Ref=[MRP3Ref NaN MRPRef(3,sw(end):end)];
    tvecRef=[tvecRef NaN t(sw(end):end)];
else
    tvecRef=t;
    MRP1Ref=MRPRef(1,:);
    MRP2Ref=MRPRef(2,:);
    MRP3Ref=MRPRef(3,:);
end

% Clip MRPref
tvecRef = tvecRef(1:end-1);
MRP1Ref = MRP1Ref(1:end-1);
MRP2Ref = MRP2Ref(1:end-1);
MRP3Ref = MRP3Ref(1:end-1);

% Plots
figure
subplot(3,1,1)
hold on
plot(tvec,MRP1,'-b')
plot(tvecRef,MRP1Ref,'--r')
hold off
title(titleStr)
xlim([t(1) t(end)])
ylabel('\sigma_1')
%ylim([-1 1])
%legend('\sigma','\sigma_r','Location','Northeast')

subplot(3,1,2)
hold on
plot(tvec,MRP2,'-b')
plot(tvecRef,MRP2Ref,'--r')
hold off
xlim([t(1) t(end)])
ylabel('\sigma_2')
%ylim([-1 1])

subplot(3,1,3)
hold on
plot(tvec,MRP3,'-b')
plot(tvecRef,MRP3Ref,'--r')
hold off
xlim([t(1) t(end)])
ylabel('\sigma_3')
%ylim([-1 1])
legend('\sigma','\sigma_r','Location','Northeast')

xlabel('Time [s]')

if psflag
    print([plotFol plotName],'-depsc')
end

end
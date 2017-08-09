%--------------------------------
% Plots results from the PolarCube ADCS simulation
%--------------------------------
% Franklin Hinckley
% 23 October 2015
%--------------------------------
% Inputs:
%   simResults (struct): 
%   SIM (struct):
% Outputs:
%   
%--------------------------------

function ADCS_Plot(simResults,SIM,SV)

% Pull out plots struct
plots = SIM.plots;

% Pull out checks
H0 = simResults.checks.H0;
T0 = simResults.checks.T0;

% Define number of points to skip at start of data series
psk = 25;

% Time limits
tmin = simResults.t(psk);
tmax = simResults.t(end);

% True and Reference MRPs
if plots.MRP == 1
    MRPplot(simResults.MRP,simResults.refMRP,simResults.t,'','','',0)
    drawnow
end


% % True and Estimated MRPs
% if plots.MRPest == 1
%     MRPplot(simResults.MRP_est,simResults.MRP,simResults.t,'','','',0)
%     drawnow
% end

mrpErr = zeros(3,length(simResults.t));
estErr = zeros(size(simResults.t));
for ii = 1:length(simResults.t)
    mrpErr(:,ii) = MRPdiff(simResults.MRP_est(:,ii),simResults.MRP(:,ii));
    Cerr = mrp2dcm(mrpErr(:,ii));
    [estErr(ii),~] = dcm2prv(Cerr,'ntol',1);
    saveas(gcf, 'PLOTS/mrps', 'jpg');
end

% Attitude determination error
if plots.detErr == 1
    figure('Name','Determination Error')
    semilogy(simResults.t(psk:end),estErr(psk:end))
    hold on
    semilogy(simResults.t(psk:end),simResults.knowReq(psk:end),'.r',...
        'Linewidth',2)
    title('Real-Time Knowledge Error')
    xlabel('Simulation Time [s]')
    ylabel('Knowledge Error [rad]')
    xlim([0 max(simResults.t)])
    drawnow
end
saveas(gcf, 'PLOTS/att_det_error', 'jpg');

% Compute 3-sigma covariance envelope for attitude errors
covEnv = zeros(6,length(simResults.t));
for ii = 1:length(simResults.t)
    for jj = 1:6
        covEnv(jj,ii) = 3*sqrt(simResults.P(jj,jj,ii));
    end
end

% Compute gyroscope bias error
bErr = zeros(3,length(simResults.t));
for ii = 1:length(simResults.t)
    bErr(:,ii) = simResults.b_est(:,ii) - simResults.bTrue(:,ii);
end

%State errors with 3-sigma covariance envelope
if plots.cov
    % Attitude errors
    figure('Name','Attitude Err/Cov')
    subplot(3,1,1)
    hold on
    plot(simResults.t(psk:end),mrpErr(1,psk:end),'.b')
    plot(simResults.t(psk:end),covEnv(1,psk:end),'--k')
    plot(simResults.t(psk:end),-covEnv(1,psk:end),'--k')
    hold off
    ylabel('\sigma_1')
    xlim([0 max(simResults.t)])
    ylim([-0.1 0.1])
    subplot(3,1,2)
    hold on
    plot(simResults.t(psk:end),mrpErr(2,psk:end),'.b')
    plot(simResults.t(psk:end),covEnv(2,psk:end),'--k')
    plot(simResults.t(psk:end),-covEnv(2,psk:end),'--k')
    hold off
    ylabel('\sigma_2')
    xlim([0 max(simResults.t)])
    ylim([-0.1 0.1])
    subplot(3,1,3)
    hold on
    plot(simResults.t(psk:end),mrpErr(3,psk:end),'.b')
    plot(simResults.t(psk:end),covEnv(3,psk:end),'--k')
    plot(simResults.t(psk:end),-covEnv(3,psk:end),'--k')
    hold off
    ylabel('\sigma_3')
    xlim([0 max(simResults.t)])
    ylim([-0.1 0.1])
    xlabel('Simulation Time [s]')
    saveas(gcf, 'att_err_conv', 'jpg')
    % Gyroscope bias errors
    figure('Name','Gyro Err/Cov')
    subplot(3,1,1)
    hold on
    plot(simResults.t(psk:end),bErr(1,psk:end),'.b')
    plot(simResults.t(psk:end),covEnv(4,psk:end),'--k')
    plot(simResults.t(psk:end),-covEnv(4,psk:end),'--k')
    hold off
    ylabel('b_1')
    xlim([0 max(simResults.t)])
    ylim([-1e-3 1e-3])
    subplot(3,1,2)
    hold on
    plot(simResults.t(psk:end),bErr(2,psk:end),'.b')
    plot(simResults.t(psk:end),covEnv(5,psk:end),'--k')
    plot(simResults.t(psk:end),-covEnv(5,psk:end),'--k')
    hold off
    ylabel('b_2')
    xlim([0 max(simResults.t)])
    ylim([-1e-3 1e-3])
    subplot(3,1,3)
    hold on
    plot(simResults.t(psk:end),bErr(3,psk:end),'.b')
    plot(simResults.t(psk:end),covEnv(6,psk:end),'--k')
    plot(simResults.t(psk:end),-covEnv(6,psk:end),'--k')
    hold off
    ylabel('b_3')
    xlim([0 max(simResults.t)])
    ylim([-1e-3 1e-3])
    xlabel('Simulation Time [s]')
    saveas(gcf, 'PLOTS/gyro_bias_err', 'jpg')
end  

% True angular velocity
if plots.omega == 1
    figure('Name','Angular Rates')
    plot(simResults.t(psk:end),flipud(simResults.omega(:,psk:end)))
    %title('\omega')
    xlabel('Simulation Time [s]')
    ylabel('Angular Velocity [rad/s]')
    %ylim([-0.015 0.025])
    legend('\omega_3','\omega_2','\omega_1');
    xlim([0 max(simResults.t)])
    drawnow
    saveas(gcf, 'PLOTS/true_ang_vel', 'jpg')
end

% Estimated angular velocity
if plots.omegaEst == 1
    figure('Name','Estimated Angular Rates')
    plot(simResults.t(psk:end),flipud(simResults.wBar(:,psk:end)))
    xlabel('Simulation Time [s]')
    ylabel('Estimated Angular Velocity [rad/s]')
    legend('\omega_3','\omega_2','\omega_1');
    xlim([0 max(simResults.t)])
    drawnow
    saveas(gcf, 'PLOTS/est_ang_vel', 'jpg')
end

if plots.omegaErr == 1
    velErr = zeros(3,length(simResults.t));
    for ii = 1:length(simResults.t)
        velErr(:,ii) = simResults.omega(:,ii) - simResults.wBar(:,ii);
    end
    figure('Name','Angular Rate Error')
    plot(simResults.t(psk:end),flipud(velErr(:,psk:end)))
    xlabel('Simulation Time [s]')
    ylabel('Angular Velocity Error [rad/s]')
    legend('\delta\omega_3','\delta\omega_2','\delta\omega_1');
    xlim([0 max(simResults.t)])
    drawnow
    saveas(gcf, 'PLOTS/ang_err_rate', 'jpg')
end

% Gyroscope biases
if plots.gyroBias == 1 
    % Plots
    figure('Name','Gyro Bias')
    subplot(3,1,1)
    hold on
    plot(simResults.t(psk:end),flipud(simResults.b_est(1,psk:end)))
    plot(simResults.t(psk:end),simResults.bTrue(1,psk:end),'r')
    hold off
    xlim([0 max(simResults.t)])
    title('Estimated Gyro 1 Bias [rad/s]')
    subplot(3,1,2)
    hold on
    plot(simResults.t(psk:end),flipud(simResults.b_est(2,psk:end)))
    plot(simResults.t(psk:end),simResults.bTrue(2,psk:end),'r')
    hold off
    xlim([0 max(simResults.t)])
    subplot(3,1,3)
    hold on
    plot(simResults.t(psk:end),flipud(simResults.b_est(3,psk:end)))
    plot(simResults.t(psk:end),simResults.bTrue(3,psk:end),'r')
    hold off  
    xlabel('Simulation Time [s]')
    xlim([0 max(simResults.t)])
    drawnow
    saveas(gcf, 'PLOTS/gyro_bias', 'jpg')
end

% Post-fit residuals
figure('Name','Post-Fit Residuals')
subplot(3,1,1)
hold on
plot(simResults.t(psk:end),simResults.PFR(1,psk:end),'.')
plot([simResults.t(psk) simResults.t(end)],...
    [3*SV.starTrk.noise 3*SV.starTrk.noise],'--r')
plot([simResults.t(psk) simResults.t(end)],...
    [-3*SV.starTrk.noise -3*SV.starTrk.noise],'--r')
hold off
ylabel('\sigma_1 [ ]')
ylim([-3.5*SV.starTrk.noise 3.5*SV.starTrk.noise])
xlim([0 max(simResults.t)])

subplot(3,1,2)
hold on
plot(simResults.t(psk:end),simResults.PFR(2,psk:end),'.')
plot([simResults.t(psk) simResults.t(end)],...
    [3*SV.starTrk.noise 3*SV.starTrk.noise],'--r')
plot([simResults.t(psk) simResults.t(end)],...
    [-3*SV.starTrk.noise -3*SV.starTrk.noise],'--r')
hold off
ylabel('\sigma_2 [ ]')
ylim([-3.5*SV.starTrk.noise 3.5*SV.starTrk.noise])
xlim([0 max(simResults.t)])

subplot(3,1,3)
hold on
plot(simResults.t(psk:end),simResults.PFR(3,psk:end),'.')
plot([simResults.t(psk) simResults.t(end)],...
    [3*SV.starTrk.noise 3*SV.starTrk.noise],'--r')
plot([simResults.t(psk) simResults.t(end)],...
    [-3*SV.starTrk.noise -3*SV.starTrk.noise],'--r')
hold off
ylabel('\sigma_3 [ ]')
ylim([-3.5*SV.starTrk.noise 3.5*SV.starTrk.noise])
xlim([0 max(simResults.t)])
xlabel('Simulation Time [s]')
saveas(gcf, 'PLOTS/post_fit_resid', 'jpg')

% Compute angular momentum magnitude 
hnorm = zeros(length(simResults.t),1);
hSVnorm = zeros(length(simResults.t),1);
hVSCMGnorm = zeros(length(simResults.t),1);
hMnorm = zeros(length(simResults.t),1);
for ii = 1:length(simResults.t)
    hnorm(ii) = norm(simResults.h(:,ii));
    hSVnorm(ii) = norm(simResults.hSV(:,ii));
    hVSCMGnorm(ii) = norm(simResults.hVSCMG(:,ii));
    hMnorm(ii) = norm(simResults.hM(:,ii));
end

% Variation in system angular momentum
if plots.Hrel == 1
    figure
    plot(simResults.t(psk:end),(hnorm(psk:end)-norm(H0))/norm(H0))
    %title('Angular Momentum Relative Error')
    xlabel('Simulation Time [s]')
    ylabel('Angular Momentum Norm [kg-m/s]')
    drawnow
    saveas(gcf, 'PLOTS/ang_mom_rel_error', 'jpg')
end

% System angular momentum
if plots.H == 1
    figure
    plot(simResults.t(psk:end),hnorm(psk:end))
    title('Angular Momentum')
    xlabel('Simulation Time [s]')
    ylabel('Angular Momentum Norm [kg-m/s]')
    drawnow
    saveas(gcf, 'PLOTS/sys_ang_mom', 'jpg')
end

% Satellite bus angular momentum
if plots.HSV == 1
    figure
    plot(simResults.t(psk:end),hSVnorm(psk:end))
    title('hSV')
    drawnow
    saveas(gcf, 'PLOTS/bus_ang_mom', 'jpg')
end

% Reaction wheel angular momentum
if plots.HRW == 1
    figure
    plot(simResults.t(psk:end),hVSCMGnorm(psk:end))
    title('hVSCMG')
    drawnow
    saveas(gcf, 'PLOTS/wheel_ang_mom', 'jpg')
end

% % Mirror angular momentum
% if plots.HM == 1
%     figure
%     plot(simResults.t(psk:end),hMnorm(psk:end))
%     title('hM')
%     drawnow
% end

% % Variation in system kinetic energy
% if plots.T == 1
%     figure
%     plot(simResults.t(psk:end),(simResults.T(psk:end)-T0)/T0)
%     %title('Kinetic Energy Variation')
%     xlabel('Simulation Time [s]')
%     ylabel('Kinetic Energy [J]')
%     drawnow
% end

% Power
if plots.P == 1
    figure
    hold on
    plot(simResults.t(psk:end),simResults.Pth(psk:end))
    plot(simResults.t(psk:300:end),simResults.Pan(psk:300:end),'.','MarkerSize',10)
    hold off
    legend('Analytical','Numerical')
    %title('Power')
    xlabel('Simulation Time [s]')
    ylabel('Power [W]')
    drawnow
    saveas(gcf, 'PLOTS/power', 'jpg')
end

% Reaction wheel torques
if plots.us == 1
    figure('Name','RW Control Torque')
    %set(gca, 'ColorOrder', [0 0 1; 1 0 0; 0 0 0], 'NextPlot', 'replacechildren');
    plot(simResults.t(psk:end),flipud(simResults.us(:,psk:end)*1000))
    %title('u_s')
    ylabel('Control Torque [mNm]')
    ylim([-0.3 0.3])
    xlabel('Simulation Time [s]')
    legend('u_3','u_2','u_1');
    xlim([0 max(simResults.t)])
    drawnow
    saveas(gcf, 'PLOTS/wheel_torques', 'jpg')
end

% Wheel speeds
if plots.Omega == 1
    figure('Name','Wheel Speeds')
    plot(simResults.t(psk:end),flipud(simResults.Omega(:,psk:end))* 9.5493)
    %title('\Omega')
    xlabel('Simulation Time [s]')
    ylabel('Wheel Speed [rpm]')
    legend('\Omega_3','\Omega_2','\Omega_1');
    xlim([simResults.t(psk) simResults.t(end)])
    drawnow
    saveas(gcf, 'PLOTS/wheel_speeds', 'jpg')
end

% Wheel accelerations
if plots.OmegaDot == 1
    figure
    plot(simResults.t(psk:end),simResults.OmegaDot(:,psk:end))
    title('OmegaDot')
    xlabel('Simulation Time [s]')
    ylabel('Wheel Acceleration [rad/s^2]')
    drawnow
    saveas(gcf, 'PLOTS/wheel_accelerations', 'jpg')
end

% % Mirror angular position
% if plots.mu == 1
%     figure
%     plot(simResults.t(psk:end),simResults.mu(psk:end))
%     title('\mu')
%     xlabel('Simulation Time [s]')
%     ylabel('Mirror Angle [rad]')
%     drawnow
% end

% % Mirror rotation rate
% if plots.muDot == 1
%     figure
%     plot(simResults.t(psk:end),simResults.muDot(psk:end))
%     title('\mu dot')
%     xlabel('Simulation Time [s]')
%     ylabel('Mirror Rate [rad/s]')
%     drawnow
% end

% Spectral content of control
% Y = fft(simResults.us(length(simResults.us)/2:end));
% L = length(Y);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% Fs = SIM.simFreq;
% f = Fs*(0:(L/2))/L;
% figure
% plot(f,P1)
% %title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P(f)|')

% Compute error (using Principal Rotation Angle)
phiErr = zeros(size(simResults.t));
for ii = 1:length(simResults.t)
    mrpErr = MRPdiff(simResults.MRP(:,ii),simResults.refMRP(:,ii));
    Cerr = mrp2dcm(mrpErr);
    [phiErr(ii),~] = dcm2prv(Cerr,'ntol',1);
end

% Principal rotation angle error 
if plots.Phi == 1
    figure
    plot(simResults.t(psk:end),phiErr(psk:end))
    title('\Phi')
    xlabel('Simulation Time [s]')
    ylabel('Attitude Error [rad]')
    drawnow
end

% Log-scale pointing error
if plots.pointErr == 1
    figure
    semilogy(simResults.t(psk:end),phiErr(psk:end))
    hold on
    semilogy(simResults.t(psk:end),simResults.pointReq(psk:end),'.r',...
        'Linewidth',2)
    %title('\Phi')
    title('Pointing Error')
    xlabel('Simulation Time [s]')
    ylabel('Attitude Error [rad]')
    xlim([0 max(simResults.t)])
    drawnow
end
 
%% Detumble mode
if plots.detumble
    detInd = SIM.detTime*SIM.simFreq;
    MRPplot(simResults.MRP(:,1:detInd),simResults.refMRP(:,1:detInd),...
        simResults.t(1:detInd),'MRPs: Detumble','','',0)
    drawnow
    saveas(gcf, 'PLOTS/mrps_DETUMBLE', 'jpg')
    figure('Name','Detumble Rates')
    hold on
    plot(simResults.t(psk:detInd)/60,flipud(simResults.omega(:,psk:detInd)))
    plot([simResults.t(psk) simResults.t(detInd)]/60,[0.2*pi/180 0.2*pi/180],...
        'k','linewidth',2)
    plot([simResults.t(psk) simResults.t(detInd)]/60,[-0.2*pi/180 -0.2*pi/180],...
        'k','linewidth',2)
    hold off
    xlabel('Simulation Time [min]')
    ylabel('Angular Velocity [rad/s]')
    ylim([-0.03 0.03])
    legend('\omega_3','\omega_2','\omega_1');
    xlim([0 simResults.t(detInd+10)/60])
    drawnow
    saveas(gcf, 'PLOTS/ang_rates_DETUMBLE', 'jpg')
end   

%% Bus mode
%Find bus indices
busStartInd = find(simResults.mode == 1,1,'first');
for ii = busStartInd:length(simResults.t)-1
    if (simResults.mode(ii) == 1) && (simResults.mode(ii+1) ~= 1)
        busStopInd = ii;
        break
    end
end

if isempty(busStopInd) && ~isempty(busStartInd)
    busStopInd = length(simResults.t);
end

if plots.bus == 1
    % Plot MRPs
    MRPplot(simResults.MRP(:,busStartInd:busStopInd),...
        simResults.refMRP(:,busStartInd:busStopInd),...
        simResults.t(busStartInd:busStopInd),'MRPs: Sun-Soak','','',0)
    drawnow
    saveas(gcf, 'PLOTS/mrps_SCI', 'jpg')
    
    % Plot angular velocity
    figure
    plot(simResults.t(busStartInd:busStopInd),...
        flipud(simResults.omega(:,busStartInd:busStopInd)))
    title('Angular Velocity: Sun-Soak')
    xlabel('Simulation Time [sec]')
    ylabel('Angular Velocity [rad/s]')
    legend('\omega_3','\omega_2','\omega_1');
    xlim([simResults.t(busStartInd) simResults.t(busStopInd)])
    drawnow
    saveas(gcf, 'PLOTS/ang_vel_SCI', 'jpg')
    
    % Plot log-scale knowledge error
    figure
    semilogy(simResults.t(busStartInd:busStopInd),...
        estErr(busStartInd:busStopInd))
    hold on
    semilogy(simResults.t(busStartInd:busStopInd),...
        simResults.knowReq(busStartInd:busStopInd),'.k',...
        'Linewidth',2)
    title('Knowledge Error: Sun-Soak')
    xlabel('Simulation Time [s]')
    ylabel('Knowledge Error [rad]')
    xlim([simResults.t(busStartInd) simResults.t(busStopInd)])
    drawnow
    saveas(gcf, 'PLOTS/log_scale_know_err_BUS', 'jpg')
    
    % Plot log-scale pointing error
    figure('Name','Pointing Error: Sun-Soak')
    semilogy(simResults.t(busStartInd:busStopInd),...
        phiErr(busStartInd:busStopInd))
    hold on
    semilogy(simResults.t(busStartInd:busStopInd),...
        simResults.pointReq(busStartInd:busStopInd),'.k',...
        'Linewidth',2)
    title('Pointing Error: Sun-Soak')
    xlabel('Simulation Time [s]')
    ylabel('Pointing Error [rad]')
    xlim([simResults.t(busStartInd) simResults.t(busStopInd)])
    drawnow 
    saveas(gcf, 'PLOTS/log_scale_point_err_BUS', 'jpg')
end

% Science mode
%Find science indices
sciStartInd = find(simResults.mode == 3,1,'first');
sciStopInd = [];
for ii = sciStartInd:length(simResults.t)
    if simResults.mode(ii) == 3
        sciStopInd = ii;
    end
end

if plots.science == 1 && ~isempty(sciStartInd) && ~isempty(sciStopInd)
    Plot MRPs
    MRPplot(simResults.MRP(:,sciStartInd:sciStopInd),...
        simResults.refMRP(:,sciStartInd:sciStopInd),...
        simResults.t(sciStartInd:sciStopInd),'MRPs: Science','','',0)
    drawnow
    saveas(gcf, 'PLOTS/mrps_SCI', 'jpg')
    
    Plot angular velocity
    figure('Name','Angular Rates: Science')
    plot(simResults.t(sciStartInd:sciStopInd),...
        flipud(simResults.omega(:,sciStartInd:sciStopInd)))
    title('Angular Velocity: Science')
    xlabel('Simulation Time [sec]')
    ylabel('Angular Velocity [rad/s]')
    legend('\omega_3','\omega_2','\omega_1');
    xlim([simResults.t(sciStartInd) simResults.t(sciStopInd)])
    drawnow
    saveas(gcf, 'PLOTS/ang_vel_SCI', 'jpg')
    
    Plot log-scale knowledge error
    figure('Name','Determination Error: Science')
    semilogy(simResults.t(sciStartInd:sciStopInd),...
        estErr(sciStartInd:sciStopInd))
    hold on
    semilogy(simResults.t(sciStartInd:sciStopInd),...
        simResults.knowReq(sciStartInd:sciStopInd),'.k',...
        'Linewidth',2)
    title('Knowledge Error: Science')
    xlabel('Simulation Time [s]')
    ylabel('Knowledge Error [rad]')
    xlim([simResults.t(sciStartInd) simResults.t(sciStopInd)])
    drawnow
    saveas(gcf, 'PLOTS/log_scale_know_err_SCI', 'jpg')
    
    Plot log-scale pointing error
    figure('Name','Pointing Error: Science')
    semilogy(simResults.t(sciStartInd:sciStopInd),...
        phiErr(sciStartInd:sciStopInd))
    hold on
    semilogy(simResults.t(sciStartInd:sciStopInd),...
        simResults.pointReq(sciStartInd:sciStopInd),'.k',...
        'Linewidth',2)
    title('Pointing Error: Science')
    xlabel('Simulation Time [s]')
    ylabel('Pointing Error [rad]')
    xlim([simResults.t(sciStartInd) simResults.t(sciStopInd)])
    drawnow 
    saveas(gcf, 'PLOTS/log_scale_point_err_SCI', 'jpg')
end

%% COM mode
% Find COM indices
comStartInd = find(simResults.mode == 2,1,'first');
comStopInd = [];
for ii = comStartInd : length(simResults.t)
    if simResults.mode(ii) == 2
        comStopInd = ii;
    end
end

if plots.COM == 1 && ~isempty(comStartInd) && ~isempty(comStopInd)
    % Plot MRPs
    MRPplot(simResults.MRP(:,comStartInd:comStopInd),...
        simResults.refMRP(:,comStartInd:comStopInd),...
        simResults.t(comStartInd:comStopInd),'MRPs: COM','','',0)
    drawnow
    
    % Plot angular velocity
    figure
    plot(simResults.t(comStartInd:comStopInd),...
        flipud(simResults.omega(:,comStartInd:comStopInd)))
    title('Angular Velocity: COM')
    xlabel('Simulation Time [sec]')
    ylabel('Angular Velocity [rad/s]')
    legend('\omega_3','\omega_2','\omega_1');
    xlim([simResults.t(comStartInd) simResults.t(comStopInd)])
    drawnow
    saveas(gcf, 'PLOTS/angular_velocity_COM', 'jpg')
    
    % Plot log-scale knowledge error
    figure
    semilogy(simResults.t(comStartInd:comStopInd),...
        estErr(comStartInd:comStopInd))
    hold on
    semilogy(simResults.t(comStartInd:comStopInd),...
        simResults.knowReq(comStartInd:comStopInd),'.k',...
        'Linewidth',2)
    title('Knowledge Error: COM')
    xlabel('Simulation Time [s]')
    ylabel('Knowledge Error [rad]')
    xlim([simResults.t(comStartInd) simResults.t(comStopInd)])
    drawnow
    saveas(gcf, 'PLOTS/log_scale_know_err', 'jpg')
    
    % Plot log-scale pointing error
    figure('Name','Pointing Error: COM')
    semilogy(simResults.t(comStartInd:comStopInd),...
        phiErr(comStartInd:comStopInd))
    hold on
    semilogy(simResults.t(comStartInd:comStopInd),...
        simResults.pointReq(comStartInd:comStopInd),'.k',...
        'Linewidth',2)
    title('Pointing Error: COM')
    xlabel('Simulation Time [s]')
    ylabel('Pointing Error [rad]')
    xlim([simResults.t(comStartInd) simResults.t(comStopInd)])
    drawnow 
    saveas(gcf, 'PLOTS/log_scale_point_err_COM', 'jpg')
end

%% Sun plots
% Sun exclusion angle for star tracker
figure('Name','Sun Exclusion Angle')
plot(simResults.t(psk:end),simResults.excAngle(psk:end)*180/pi)
xlabel('Simulation Time [s]')
ylabel('Exclusion Angle [deg]')
saveas(gcf, 'PLOTS/sun_exclus_angle', 'jpg')

% Angle between solar panel normal and sun
figure('Name','Solar Panel Incidence Angle')
plot(simResults.t(psk:end),simResults.sunAngle(psk:end)*180/pi)
xlabel('Simulation Time [s]')
ylabel('Angle Off \perp [deg]')
ylim([-0.1 10])
saveas(gcf, 'PLOTS/sol_panel_inc_angle', 'jpg')

%% Other mode plots
%Angle between S-Band antenna axis and ground station
figure('Name','S-Band Tracking Error')
plot(simResults.t(psk:end),simResults.txAngle(psk:end)*180/pi)
xlabel('Simulation Time [s]')
ylabel('Tracking Error [deg]')
ylim([-0.1 10])
xlim([tmin tmax])
saveas(gcf, 'PLOTS/s_band_track_err', 'jpg')


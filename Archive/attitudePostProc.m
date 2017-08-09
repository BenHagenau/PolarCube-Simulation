%------------------------------
% Smooths the attitude data
%------------------------------
% Franklin Hinckley
% 15 January 2016
%------------------------------
% Inputs:
%   
% Outputs:
%   
%------------------------------

function attitudePostProc(simResults)

%% Pull out simulation results
MRP = simResults.MRP_est;
wBar = simResults.wBar;
P = simResults.P(1:3,1:3,:);

%% Smoother
% Initialize output
[~,ll] = size(MRP);
MRP_S = zeros(size(MRP));
P_S = zeros(size(P));

MRP_S(:,end) = MRP(:,end);
P_S(:,:,end) = P(:,:,end);

% Smoother loop
for kk = ll-1 : -1 : 1
    % Compute smoother gain
    F = 0.5*(MRP(:,kk)*wBar(:,kk)' - wBar(:,kk)*MRP(:,kk)' - ...
        tilde(wBar(:,kk)) + MRP(:,kk)'*wBar(:,kk)*eye(3));
    C = P(:,:,kk)*F'*P_S(:,:,kk+1);
    
    % Smooth state and covariance
    MRP_S(:,kk) = MRP(:,kk) + C*(MRP_S(:,kk+1) - MRP(:,kk+1));
    P_S(:,:,kk) = P(:,:,kk) + C*(P_S(:,:,kk+1)-P(:,:,kk))*C';
end

%% Plot results
psk = 25;
figure
plot(simResults.t,MRP_S)

phiErr = zeros(size(simResults.t));
for ii = 1:length(simResults.t)
    mrpErr = MRPdiff(MRP_S(:,ii),simResults.MRP(:,ii));
    Cerr = mrp2dcm(mrpErr);
    [phiErr(ii),~] = dcm2prv(Cerr,'ntol',1);
end

figure
semilogy(simResults.t(psk:end),phiErr(psk:end))
hold on
semilogy([simResults.t(psk) simResults.t(end)],[0.13*pi/180 0.13*pi/180],...
        'Linewidth',2)
%title('\Phi')
title('Knowledge Error')
xlabel('Simulation Time [s]')
ylabel('Knowledge Error [rad]')
xlim([0 max(simResults.t)])
drawnow

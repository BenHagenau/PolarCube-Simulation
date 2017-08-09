function [OmegaDot,SV]=ADCS_Control_C(errMRP,z,omega,errOmega,omegaRef,omegaDotRef, L, SV, SIM)
                       
RW_torq = [];
MAG_torq = [];

K_mex = SV.K;
P_mex = diag(SV.P)';
Km_mex = SV.Km;
Ki_mex = diag(SV.Ki)';

%Launching c-code - parameters of s/c and simulation are hard-coded
[RW_torq,MAG_torq] = control(errMRP,z,omega,errOmega,omegaRef,omegaDotRef,...
    L, SV.Omega, SIM.Bbody, K_mex, P_mex, Km_mex, Ki_mex, SV.deadband, RW_torq, MAG_torq);

[~,N] = size(SV.Gs);
Lr = RW_torq';
SV.Lt = MAG_torq';

% Define D matrix
D0 = zeros(3,N);
for ii = 1:N
    D0(:,ii) = SV.Gs(:,ii)*SV.Js(ii);
end
    
% Define Q matrix
Q = D0;

% Defing weighting matrix
W = eye(N);

% Check for no control
if (SIM.RWmode == 0)
    OmegaDot = zeros(N,1);
else
    % Solve for desired control eta-dot (eta = [Omega; gamma])
    OmegaDot = W*Q'*((Q*W*Q')^(-1))*Lr;
end
end

%------------------------------
% B-Dot Controller for detumble
%------------------------------
% Franklin A. Hinckley
% Initial Version: 9 March 2016
% Latest Revision: 1.0
%   9 March 2016
%------------------------------
% Inputs:
%   
% Outputs:
%   
%------------------------------

function [OmegaDot, SV] = Bdot_Controller(Bt1,Bt0,omega,MRP,SV,SIM)

%% Parse satellite struct
% Control gains
Kb = SV.Kb;

% Pull out spin axis matrix
Gs = SV.Gs;

% Pull out spin axis inertia
Js = SV.Js;

%% Evaluate control law
% Determine Bdot
tstep = 1/SIM.simFreq;
Bdot = (Bt1 - Bt0)/tstep;
Bdot = Bdot * 1e-5;
%Bdot = cross(Bdot,Bt1)/norm(Bt1);

%SV.Bdot = Bdot;

Bdot = mrp2body(Bt1,Bdot*1e9);

%Bdot = tilde(omega)*Bt1*1e-5;

% Compute required torque
Lr = -Kb * Bdot;
Lr = cross(Lr,Bt1)/norm(Bt1);

%% Compute reaction wheel control command
% Define D matrix
D0 = zeros(3,3);
for ii = 1:3
    D0(:,ii) = Gs(:,ii)*Js(ii);
end
    
% Define Q matrix
Q = D0;

% Defing weighting matrix
W = eye(3);

% Check for no control
if (SIM.RWmode == 0)
    OmegaDot = zeros(3,1);
else
    % Solve for desired control eta-dot (eta = [Omega])
    OmegaDot = W*Q'*((Q*W*Q')^(-1))*Lr;
end
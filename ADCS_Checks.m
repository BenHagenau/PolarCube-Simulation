%--------------------------------
% Simulation value checks (angular momentum, kinetic energy, power)
%--------------------------------
% Franklin Hinckley
% 23 October 2015
%--------------------------------

function [checks] = ADCS_Checks(simResults,SV,tstep)

% Pull out values
t = simResults.t;
omega = simResults.omega;
Omega = simResults.Omega;
L = simResults.Lsim;
us = simResults.us;

% Initialize return values
h = zeros(length(t),1);
T = zeros(length(t),1);
Pan = zeros(length(t),1);
Pnum = zeros(length(t),1);

% Compute initial values
H0 = SV.I*omega(:,1);
T0 = 0.5*omega(:,1)'*SV.I*omega(:,1);
for jj = 1:3       
    ws = SV.Gs(:,jj)'*omega(:,1);
    H0 = H0 + SV.Js(jj)*(ws + Omega(jj,1))*SV.Gs(:,jj);
    T0 = T0 + 0.5*(SV.Js(jj)*((Omega(jj,1)+ws)^2));
end
checks.H0 = norm(H0);
checks.T0 = norm(T0);

% Loop over all times
for ii = 2:length(t)
    % Compute angular momentum
    hRW = 0;
    for jj = 1 : 3        
        ws = SV.Gs(:,jj)'*omega(:,ii);
        hRW = hRW + SV.Js(jj)*(ws + Omega(jj,ii))*SV.Gs(:,jj);  
    end  
    h(ii) = norm(hRW + SV.I*omega(:,ii));

    % Compute kinetic energy
    Ttmp = 0;
    for jj = 1 : 3        
        ws = SV.Gs(:,jj)'*omega(:,ii);
        Ttmp = Ttmp + SV.Js(jj)*((Omega(jj,ii)+ws)^2);        
    end     
    T(ii) = 0.5*Ttmp + 0.5*omega(:,ii)'*SV.I*omega(:,ii);
end

% Compute power
for ii = 2:length(t)
    Ptmp = 0;
    for jj = 1 : 3        
        Ptmp = Ptmp + Omega(jj,ii)*us(jj,ii);          
    end  
    Pan(ii) = omega(:,ii)'*L(:,ii) + Ptmp;
    Pnum(ii) = (T(ii) - T(ii-1))/tstep;  
end

% Assign outputs
checks.hnorm = h;
checks.T     = T;
checks.P_an  = Pan;
checks.P_num = Pnum;


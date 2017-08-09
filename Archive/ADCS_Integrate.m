%------------------------------
% Integrates the PolarCube Equations of Motion
%------------------------------
% Franklin Hinckley
% 30 October 2015
%------------------------------
% Inputs:
%   MRP (3 x 1 double vector): spacecraft attitude  in MRPs
%   omega (3 x 1 double vector): spacecraft angular rate [rad/s]
%   L (3 x 1 double array): total external torques [Nm]
%   SIM (struct): simulation properties
%   SV (struct): satellite properties
%   tstep (double): integration time step [sec]
%   t (double): current simulation time [sec]
% Outputs:
%   MRP (3 x 1 double vector): new spacecraft attitude  in MRPs
%   omega (3 x 1 double vector): new spacecraft angular rate [rad/s]
%   SIM (struct): updated simulation properties
%   SV (struct): updated satellite properties
%------------------------------

function [MRP,omega,SV,SIM] = ADCS_Integrate(MRP,omega,L,SIM,SV,tstep,t)

% Pull out G Matrices
Gs = SV.Gs;
Gt1 = SV.Gt1;
Gt2 = SV.Gt2;

% Pull out J vectors
Js = SV.Js;
Jt1 = SV.Jt1;
Jt2 = SV.Jt2;

% Pull out mirror parameters
Im = SV.Im;
BM = SV.BM;
muDot = SV.muDot(t);
SV.muDotA = muDot;
muDoubDot = SV.muDoubDot(t);

% Compute summation of inertias
JggS = 0;
JggT = 0;
JggG = 0;
for ii = 1:3
    JggS = JggS + Js(ii)*Gs(:,ii)*Gs(:,ii)';
    JggT = JggT + Jt1(ii)*Gt1(:,ii)*Gt1(:,ii)';
    JggG = JggG + Jt2(ii)*Gt2(:,ii)*Gt2(:,ii)';
end
I_VSCMG = JggS + JggT + JggG;

% Compute total inertia
I = SV.I + I_VSCMG + BM*Im*BM';

% Pull out rates
Omega = SV.Omega;
OmegaDot = SV.OmegaDot;

% Compute velocity projections
ws=zeros(3,1);
wt=zeros(3,1);
wg=zeros(3,1);
for ii = 1:3
    ws(ii) = Gs(:,ii)'*omega;
    wt(ii) = Gt1(:,ii)'*omega;
    wg(ii) = Gt2(:,ii)'*omega;
end

% Compute motor torques
us = zeros(3,1);
for ii = 1:3
    us(ii) = Js(ii)*(OmegaDot(ii)+Gs(:,ii)'*SIM.omegaDot);
end

% Check torque limit
for ii = 1:3
    if abs(us(ii)) > SV.RWtorqueLim;
        us(ii) = sign(us(ii))*SV.RWtorqueLim;
    end
end

SV.us = us;

JsM = diag(Js);
% Compute system mass matrix
M = [eye(3)     zeros(3,3) zeros(3,3) ;...
     zeros(3,3) I          JsM*Gs     ;...
     zeros(3,3) JsM*Gs'    JsM*eye(3)];
 
% Define differential equations
f_sigma    = body2mrp(omega,MRP);
fomegaW = 0;
for ii = 1:3
    fomegaW = fomegaW + Js(ii)*tilde(omega)*Omega(ii)*Gs(:,ii);
end

omegaMB = muDot*SV.BM(:,3);
fomegaM = tilde(omegaMB)*BM*Im*BM'*(omega+omegaMB)-...
        BM*Im*BM'*tilde(omegaMB)*omega+...
        BM*Im*BM'*muDoubDot*SV.BM(:,3)+...
        tilde(omega)*BM*Im*BM'*omegaMB;

f_omega = -tilde(omega)*I*omega + L - fomegaW - fomegaM;

f_Omega    = zeros(3,1);
for ii = 1:3
    % Check RW case
    if SIM.RWdyn == 1
        f_Omega(ii) = us(ii);
    else
        f_Omega(ii) = 0;
    end
end

% Solve differential equations
Xdot = M\[f_sigma' f_omega' f_Omega']';

% Zero out torques by mode
if SIM.RWdyn == 0
    Xdot(7:9) = zeros(3,1);
end

% Integrate
MRP      = MRP      + Xdot(1:3) *tstep;
omega    = omega    + Xdot(4:6) *tstep;
SV.Omega = SV.Omega + Xdot(7:9) *tstep;
SV.mu    = SV.mu    + muDot     *tstep;

% Compute new mirror-frame rotation matrix
dm = SV.mu - SV.mu0;
ms =  cos(dm)*SV.ms0 + sin(dm)*SV.ml0;
ml = -sin(dm)*SV.ms0 + cos(dm)*SV.ml0;
mr = SV.mr0;
SV.BM = [ms ml mr];

end

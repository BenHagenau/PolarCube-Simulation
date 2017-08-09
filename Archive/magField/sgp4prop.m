%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements the SGP4 propagator as defined in Spacetrack Report #3
% Felix R. Hoots, Ronald L. Roehrich, TS Kelso 
% 31 December 1988
% https://celestrak.com/NORAD/documentation/spacetrk.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Franklin Hinckley
% Latest Revision: 18 December 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: 
%   TLE (struct): contains the orbital information read from the TLE
%   ts (double): start time for propagation [TLE epoch]
%   te (double): end time for propagation [TLE epoch]
%   dt (double): time step [decimal day]
% Outputs:
%   r (double): matrix of ECI positions [km], each column is a single
%       position vector
%   v (double): matrix of ECI velocities [km/s], each column is a 
%       single velocity vector     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rf,vf]=sgp4prop(TLE,ts,dt,te)
%% Define constants
RE=6378.135;    %earth radius [km]
k2=5.413080e-4*RE^2;
k4=.62098875e-6*RE^4;
QOMS=1.88027916e-9/(RE^4);
S=1.01222928*RE;
A30=-(-.253881e-5)*RE^3; %
TD=1440.0;  %time units per day
mu=398600.4415; %gravitational parameter of earth
ke=sqrt(mu);

%% Convert degrees to radians
dr=pi/180;
TLE.incl=TLE.incl*dr;
TLE.RAAN=TLE.RAAN*dr;
TLE.argPer=TLE.argPerigee*dr;
TLE.meanAnomaly=TLE.meanAnomaly*dr;

%% Recover original mean motion and semimajor axis
a1=(ke/TLE.meanMotion)^(2/3);
delta1=(3/2)*(k2/(a1^2))*((3*cos(TLE.incl)^2-1)/((1-TLE.ecc^2)^(3/2)));
a0=a1*(1-(delta1/3)-(delta1^2)-((134/81)*delta1^3));
delta0=(3/2)*(k2/(a0^2))*((3*cos(TLE.incl)^2-1)/((1-TLE.ecc^2)^(3/2)));
n0pp=TLE.meanMotion/(1+delta0);    %true mean motion at epoch
a0pp=a0/(1-delta0);    %true semimajor axis at epoch

%% Check perigee height and correct values if appropriate
per=(a0pp*(1-TLE.ecc))-RE;
if per<98
    ss=(20/RE)+RE;
    QOMS=((QOMS^(1/4))+S-ss)^4;
    S=ss;
elseif per<156
    ss=(a0pp*(1-TLE.ecc))-s+RE;
    QOMS=((QOMS^(1/4))+S-ss)^4;
    S=ss;
end

%% Calculate constants
theta=cos(TLE.incl);
zeta=1/(a0pp-S);
beta0=sqrt(1-TLE.ecc^2);
nu=a0pp*TLE.ecc*zeta;
C2=QOMS*zeta^4*n0pp*((1-nu^2)^(-7/2))*(a0pp*(1+((3/2)*nu^2)+(4*TLE.ecc*nu)+...
    (TLE.ecc*nu^3))+((3/2)*((k2*zeta)/(1-nu^2))*((-1/2)+((3/2)*theta^2))*...
    (8+(24*nu^2)+(3*nu^4))));
C1=TLE.BSTAR*C2;
C3=(QOMS*zeta^5*A30*n0pp*RE*sin(TLE.incl))/(k2*TLE.ecc);
C4=(2*n0pp*QOMS*zeta^4*a0pp*beta0^2*((1-nu^2)^(-7/2)))*(((2*nu*(1+(TLE.ecc*nu)))+...
    ((1/2)*TLE.ecc)+((1/2)*nu^3))-(((2*k2*zeta)/(a0pp*(1-nu^2)))*((3*(1-3*theta^2))*...
    (1+((3/2)*nu^2)-(2*TLE.ecc*nu)-((1/2)*TLE.ecc*nu^3))+((3/4)*(1-theta^2)*...
    ((2*nu^2)-(TLE.ecc*nu)-(TLE.ecc*nu^3))*cos(2*TLE.argPer)))));
C5=(2*QOMS*zeta^4*a0pp*beta0^2*((1-nu^2)^(-7/2)))*(1+((11/4)*nu*(nu+TLE.ecc))+...
    (TLE.ecc*nu^3));

D2=4*a0pp*zeta*C1^2;
D3=(4/3)*a0pp*zeta^2*((17*a0pp)+S)*C1^3;
D4=((2/3)*a0pp*zeta^3)*((221*a0pp)+(31*S))*C1^4;

%% Main loop
tvec=ts:dt:te;
rf=zeros(3,length(tvec));
vf=zeros(3,length(tvec));
toe=TLE.epoch;
for c=1:length(tvec)
    t=tvec(c);
    tse=(t-toe)*86400;  %convert decimal day to seconds
    %% Secular drag and gravitation effects (f(t))
    MDF=TLE.meanAnomaly+((1+((3*k2*(3*theta^2-1))/(2*a0pp^2*beta0^3))+((3*k2^2*...
        (13-(78*theta^2)+(137*theta^4)))/(16*a0pp^4*beta0^7)))*(n0pp*tse));
    wDF=TLE.argPer+((((-3*k2*(1-(5*theta^2)))/(2*a0pp^2*beta0^4))+((3*k2^2*...
        (7-(114*theta^2)+(395*theta^4)))/(16*a0pp^4*beta0^8))+((5*k4*...
        (3-(36*theta^2)+(49*theta^4)))/(4*a0pp^4*beta0^8)))*(n0pp*tse));
    RAANDF=TLE.RAAN+((((-3*k2*theta)/(a0pp^2*beta0^4))+((3*k2^2*(4*theta-19*theta^3))/...
        (2*a0pp^4*beta0^8))+((5*k4*theta*(3-7*theta^2))/(2*a0pp^4*beta0^8)))*(n0pp*tse));
    delw=TLE.BSTAR*C3*cos(TLE.argPer)*tse;
    delM=(-2/3)*QOMS*TLE.BSTAR*zeta^4*(RE/(TLE.ecc*nu))*(((1+(nu*cos(MDF)))^3)...
        -((1+(nu*cos(TLE.meanAnomaly)))^3));
    Mp=MDF+delw+delM;
    w=wDF-delw-delM;
    Omega=RAANDF-((21/2)*((n0pp*k2*theta)/(a0pp^2*beta0^2))*C1*tse^2);
    e=TLE.ecc-(TLE.BSTAR*C4*tse)-(TLE.BSTAR*C5*(sin(Mp)-sin(TLE.meanAnomaly)));
    %should add check for perigee height at epoch
    a=a0pp*(1-(C1*tse)-(D2*tse^2)-(D3*tse^3)-(D4*tse^4))^2;
    IL=Mp+w+Omega+(n0pp*((3/2)*C1*tse^2+((D2+2*C1^2)*tse^3)+...
        ((1/4)*(3*D3+12*C1*D2+10*C1^3)*tse^4)+...
        ((1/5)*(3*D4+12*C1*D3+6*D2^2+30*C1^2*D2+15*C1^4)*tse^5)));
    beta=sqrt(1-e^2);
    n=ke/(a^(3/2));

    %% Long-period periodics
    axN=e*cos(w);
    ILL=((A30*sin(TLE.incl))/(8*k2*a*beta^2))*(e*cos(w))*((3+5*theta)/(1+theta));
    ayNL=(A30*sin(TLE.incl))/(4*k2*a*beta^2);
    ILt=IL+ILL;
    ayN=(e*sin(w))+ayNL;
    U=ILt-Omega;

    %% Solve Kepler's equation
    tol=1e-6;
    EW1=U;
    err=1;
    while(err>tol)
        delEW=(U-(ayN*cos(EW1))+(axN*sin(EW1))-EW1)/((-ayN*sin(EW1))-(axN*cos(EW1))+1);
        EW2=EW1+delEW;
        err=abs(EW2-EW1);
        EW1=EW2;
    end
    EW=EW2;

    %% Short-period periodics
    ecE=(axN*cos(EW))+(ayN*sin(EW));
    esE=(axN*sin(EW))-(ayN*cos(EW));
    el=sqrt(axN^2+ayN^2);
    pl=a*(1-el^2);
    r=a*(1-ecE);
    rdot=ke*(sqrt(a)/r)*esE;
    rfdot=ke*(sqrt(pl)/r);
    cu=(a/r)*(cos(EW)-axN+((ayN*esE)/(1+sqrt(1-el^2))));
    su=(a/r)*(sin(EW)-ayN-((axN*esE)/(1+sqrt(1-el^2))));
    u=atan2(su,cu);
    delr=(k2/(2*pl))*(1-theta^2)*cos(2*u);
    delu=(-k2/(4*pl^2))*(7*theta^2-1)*sin(2*u);
    delOm=((3*k2*theta)/(2*pl^2))*sin(2*u);
    deli=((3*k2*theta)/(2*pl^2))*sin(TLE.incl)*cos(2*u);
    delrdot=((-k2*n)/pl)*(1-theta^2)*sin(2*u);
    delrfdot=((k2*n)/pl)*(((1-theta^2)*cos(2*u))-((3/2)*(1-3*theta^2)));

    %% Compute osculating values
    rk=(r*(1-((3/2)*k2*(sqrt(1-el^2)/(pl^2))*(3*theta^2-1))))+delr;
    uk=u+delu;
    Omk=Omega+delOm;
    ik=TLE.incl+deli;
    rdotk=rdot+delrdot;
    rfdotk=rfdot+delrfdot;

    %% Compute unit vectors
    M=[-sin(Omk)*cos(ik);cos(Omk)*cos(ik);sin(ik)];
    N=[cos(Omk);sin(Omk);0];
    U=M*sin(uk)+N*cos(uk);
    V=M*cos(uk)-N*sin(uk);

    %% Compute position and velocity
    rf(:,c)=rk*U;
    vf(:,c)=rdotk*U+rfdotk*V;

end

end


%--------------------------------
% Sets the reference trajectory for the current position and mode
%--------------------------------
% Franklin Hinckley
% 30 January 2016
%--------------------------------
% Inputs:
%   SV (struct): satellite properties
%   SIM (struct): simulation properties
% Outputs:
%   SV (struct): satellite properties
%   SIM (struct): simulation properties
%   mode (string): defines current mode
%--------------------------------

function [MRPRef, omegaRef, omegaDotRef] = ...
    modeRefSet(mode, sunVec, COMvec, R)

%% Body-frame definitions
% Solar panel axis
SPaxis = [-sqrt(2)/2; -sqrt(2)/2; 0];

% Transmit antenna
TXaxis = [1; 0; 0];

%% Reference MRP
% Check mode and set appropriate targets
switch mode
    case 0 % Detumble
        BT = eye(3);
        NT = eye(3);
    case 1 % Bus priority (sun-pointing)
        % Set targets
        priTarg = sunVec;
        secTarg = [0; 1; 0];

        % Compute inertial-frame DCM for desired axes
        NT = targetFrame(priTarg, secTarg);
        
        % Compute body-frame DCM for desired axes
        BT = targetFrame(SPaxis, TXaxis);
        
        % Form reference DCM [BN]
        %Cref = BT*NT'; 
        
    case 2 % S-Band COM
        % Set targets
        priTarg = COMvec - R;
        secTarg = sunVec;

        % Compute inertial-frame DCM for desired axes
        NT = targetFrame(priTarg, secTarg);
        
        % Compute body-frame DCM for desired axes
        BT = targetFrame(TXaxis, SPaxis);
        
        % Form reference DCM
        %Cref = BT*NT';   
        
    case {3,4} % Science and pre-science
        % Set targets
        priTarg = -R;
        secTarg = sunVec;

        % Compute inertial-frame DCM for desired axes
        NT = targetFrame(priTarg, secTarg);
        
        % Form reference DCM
        %C1 = rotMatrix(2,72.72); % Pitch up so beam points nadir
        C2 = rotMatrix(3,45); % Roll so boom is aligned limb
        %Cref = C2*(C1*NT');   
        Cref = C2*NT';
        
        BT = eye(3);
        
end   

% Compute DCM for reference trajectory relative to inertial
Cref = BT*NT';

% Convert to MRP
MRPRef = dcm2mrp(Cref);
    
%% Reference angular velocity and angular acceleration
omegaRef = zeros(3,1);
omegaDotRef = zeros(3,1);

end 
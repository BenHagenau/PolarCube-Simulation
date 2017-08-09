%--------------------------------
% Sets the parameters for the current mode
%--------------------------------
% Franklin Hinckley
% 12 January 2016
%--------------------------------
% Inputs:
%   SV (struct): satellite properties
%   SIM (struct): simulation properties
% Outputs:
%   SV (struct): satellite properties
%   SIM (struct): simulation properties
%   mode (string): defines current mode
%--------------------------------

function [SV, SIM] = modeParamSet(SV, SIM, mode)

% Switch block for current mode
switch mode
    case 0 % Detumble
        % Set gains
        SV.K = 0;
        SV.P = diag([0.003; 0.003; 0.003]);
        SV.Ki = diag([0.001; 0.001; 0.001]);
        SV.Km = 0.0003;
        %SV.Km = 0;
        % Set flags
        SIM.RWdyn  = 1;
        SIM.mirrorDyn = 1;
        %SIM.RWmode  = 1;
        SIM.integralControl = 1;
        % Set mirror properties
        SV.muDot = @(t)0;
        SV.muDoubDot = @(t)0;
        % Enable control smoothing
        SV.contFilt = 0;
        % Set control deadband
        %SV.deadband = 0.01/1000;
        % SV.deadband = 1E-12;
        
    case 1 % Bus priority (sun-pointing)
        % Set gains
        SV.K = 0.005;
        SV.P = diag([0.0075; 0.0075; 0.0075]);
        SV.Ki = diag([0.0025; 0.0025; 0.0025]);
        SV.Km = 0.0003;
        %SV.Km = 0;
        % Set flags
        SIM.RWdyn  = 1;
        SIM.mirrorDyn = 1;
        %SIM.RWmode  = 1;
        SIM.integralControl = 1;
        % Set mirror properties
        SV.muDot = @(t)0;
        SV.muDoubDot = @(t)0;
        % Enable control smoothing
        SV.contFilt = 0;
        % Set control deadband
        % SV.deadband = 1E-12;

    case 2 % S-Band COM
        % Set gains
        SV.K = 0.005;
        SV.P = diag([0.0075; 0.0075; 0.0075]);
        SV.Ki = diag([0; 0; 0]);
        SV.Km = 0.0003;  
        %SV.Km = 0;
        % Set flags
        SIM.RWdyn  = 1;
        SIM.mirrorDyn = 1;
        %SIM.RWmode  = 1;
        SIM.integralControl = 1;
        % Set mirror properties
        SV.muDot = @(t)0;
        SV.muDoubDot = @(t)0;
        % Enable control smoothing
        SV.contFilt = 0;
        % Set control deadband
        % SV.deadband = 1E-12;

    case 3 % Science
        % Set gains
        SV.K = 0.0075;
        SV.P = diag([0.01; 0.01; 0.0125]);
        SV.Ki = diag([0.0025;0.0025;0.0025]);
        %SV.Ki = diag([0; 0; 0]);
        SV.Km = 0; 
        % Set flags
        SIM.RWdyn  = 1;
        SIM.mirrorDyn = 1;
        %SIM.RWmode  = 1;
        SIM.integralControl = 1;
        % Set mirror properties (including spin-up)
        SV.muDot = @(t)(0.13*2*pi*t)*(t <= 10) + (1.3*2*pi)*(t > 10);
        SV.muDoubDot = @(t)0.13*2*pi*(t <= 10);
%         SV.muDot = @(t)(0.13*2*pi*(t-0))  *((0 < t)   && (t <= 10))  + (1.3*2*pi)*((10 < t)  && (t <= 60))  - (0.13*2*pi*(t-50)) *((50 < t)  && (t <= 60))  + ...
%                        (0.13*2*pi*(t-90)) *((90 < t)  && (t <= 100)) + (1.3*2*pi)*((100 < t) && (t <= 150)) - (0.13*2*pi*(t-140))*((140 < t) && (t <= 150)) + ...
%                        (0.13*2*pi*(t-180))*((180 < t) && (t <= 190)) + (1.3*2*pi)*((190 < t) && (t <= 240)) - (0.13*2*pi*(t-230))*((230 < t) && (t <= 240)) + ...
%                        (0.13*2*pi*(t-270))*((270 < t) && (t <= 280)) + (1.3*2*pi)*((280 < t) && (t <= 330)) - (0.13*2*pi*(t-320))*((320 < t) && (t <= 330));
%         SV.muDoubDot = @(t)(0.13*2*pi)*((0 < t)   && (t <= 10))  - (0.13*2*pi)*((50 < t)  && (t <= 60))  + ...
%                            (0.13*2*pi)*((90 < t)  && (t <= 100)) - (0.13*2*pi)*((140 < t) && (t <= 150)) + ...
%                            (0.13*2*pi)*((180 < t) && (t <= 190)) - (0.13*2*pi)*((230 < t) && (t <= 240)) + ...
%                            (0.13*2*pi)*((270 < t) && (t <= 280)) - (0.13*2*pi)*((320 < t) && (t <= 330));
% %         SV.muDot = @(t)(0.025*2*pi*t)*(t <= 10) + (0.25*2*pi)*(t > 10);
%         SV.muDoubDot = @(t)(0.025*2*pi)*(t <= 10);
        % Enable control smoothing
        SV.contFilt = 1;
        % Set control deadband
        % SV.deadband = 1E-12; 
        
    case 4 % Pre-science
        % Set gains
        SV.K = 0.005;
        SV.P = diag([0.0075; 0.0075; 0.0075]);
        SV.Ki = diag([0; 0; 0]);
        SV.Km = 0; 
        % Set flags
        SIM.RWdyn  = 1;
        SIM.mirrorDyn = 1;
        %SIM.RWmode  = 1;
        SIM.integralControl = 1;
        % Set mirror properties (including spin-up)
        SV.muDot = @(t)0;
        SV.muDoubDot = @(t)0;
        % Enable control smoothing
        SV.contFilt = 0;
        % Set control deadband
        % SV.deadband = 1E-12; 
        
    case 5 % Mode transitions
        % Set gains
        SV.K = 0.005;
        SV.P = diag([0.0075; 0.0075; 0.0075]);
        SV.Ki = diag([0; 0; 0]);
        SV.Km = 0;  
        % Set flags
        SIM.RWdyn  = 1;
        SIM.mirrorDyn = 1;
        %SIM.RWmode  = 1;
        SIM.integralControl = 1;
        % Set mirror properties 
        SV.muDot = @(t)0;
        SV.muDoubDot = @(t)0;
        % Enable control smoothing
        SV.contFilt = 0;
        % Set control deadband
        % SV.deadband = 1E-12; 
end   

end 
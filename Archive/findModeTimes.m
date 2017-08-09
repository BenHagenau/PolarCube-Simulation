%------------------------------
% Runs the orbit propagator to identify times for COM and SCI modes
%------------------------------
% Franklin Hinckley
% 15 January 2016
%------------------------------
% Inputs:
%   
% Outputs:
%   
%------------------------------

%% Get simulation settings
SIM = simSettings();

%% Run propagator
SIM = orbitPropagator(SIM);

%% Iterate through data to find COM and SCI modes
[~,len] = size(SIM.prop.R);
mode = zeros(len,0);
tvec = 0 : 1 : SIM.simTime;

for ii = 1:len
    % Determine current Julian Date
    JD = SIM.JD0 + (ii/86400);
        
    % Evaluate potential pointing targets
    % Determine the sun vector
    sunVec = sun_vector(JD);

    % Determine COM and science targets
    COMvec = get_target_vector(SIM.latCOM, SIM.latCOM, JD);
    SCIvec = get_target_vector(SIM.latSCI, SIM.latSCI, JD);
    
    %% Determine current mode
    % Check ranges to targets
    COMrange = norm(SIM.prop.R(:,ii) - COMvec);
    SCIrange = norm(SIM.prop.R(:,ii) - SCIvec);
    
    % Check detumble time
    if (ii > SIM.detTime)
        % Check target ranges
        if COMrange < SIM.COMrngLim
            mode(ii) = 2;
        elseif SCIrange < SIM.SCIrngLim
            mode(ii) = 3;
        else
            mode(ii) = 1;
        end
    else
        mode(ii) = 0;
    end
end

%% Plots
figure
plot(tvec,mode(1:length(tvec)))
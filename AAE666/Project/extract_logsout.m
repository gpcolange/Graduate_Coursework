function  extract_logsout(out)
%%% Inputs: logsout - SimulationOutput object containing logged signals
%%% Outputs: Writes data in SimulationOutput object to base workspace

if ~isa(out,'Simulink.SimulationOutput')
    error('Function input must be a SimulationOutput object')
end

% Loop through each logged signal
for j = 1:numElements(out.logsout)
    % Assign data to variable 
    data_temp = out.logsout{j}.Values.Data;
    name      = out.logsout{j}.Name;
    assignin('base', name, data_temp)
end

% Write time vector to base workspace
assignin('base','tout', out.tout)
end
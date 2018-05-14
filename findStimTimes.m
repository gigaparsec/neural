% Author: Jack Tchimino

%% Find Stimulation Timings
% parse the stimulation waveform, see where the stimulation takes place,
% return a vector with time instants
% for now it is simply an application of the findpeaks function, but it
% might have things added to it, depending on the data provided
function refPoints = findStimTimes(data)
    
    [~, refPoints] = findpeaks(data);
    
    % this is not right. the stimulation is denoted with pulses, the spikes
    % are only for the oxidation of the electrodes, not to be taken into
    % consideration in the analysis
    
end
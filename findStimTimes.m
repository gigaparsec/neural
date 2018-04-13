% Author: Jack Tchimino

%% Find Stimulation Timings
% parse the stimulation waveform, see where the stimulation takes place,
% return a vector with time instants
% for now it is simply an application of the findpeaks function, but it
% might have things added to it, depending on the data provided
function refPoints = findStimTimes(data)
    
    [~, refPoints] = findpeaks(data);
    
end
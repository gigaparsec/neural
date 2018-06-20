% Author: Jack Tchimino

%% create signal struct 
% with user defined parameters. A shorter signal has a different frequency
% vector. The time vector can be easily computed outside the function.

function signalStruct = createSignalStruct(signal,fs,tV)
    
    N = length(signal);
    
    time = tV(end);

    fV = 0:1/time:fs;
    if isequal(length(fV),N)==0
        fV = 0:1/time:fs-1/time;
    end
    
    signalStruct = struct('signal',double(signal),'SamplingFrequency',fs,'TimeVector',tV,'FrequencyVector',fV);
    
end
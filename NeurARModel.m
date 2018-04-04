% Author: Jack Tchimino

%% Neural Signal Modeling with AR Processes

function ARmod = NeurARModel(data,order)

    signal = data.signal;
    tV     = data.TimeVector;
    fV     = data.FrequencyVector;
    fs     = data.SamplingFrequency;
    
    
    temp = iddata(double(signal),[],1/fs);
    ARmod = ar(temp,order);

    
    
end
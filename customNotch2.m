function notched = customNotch2(comb_b,comb_a,data)

fV = data.FrequencyVector;
tV = data.TimeVector;
fs = data.SamplingFrequency;
FileDir = data.FilePath;

numOfNotches = length(comb_b);

sig = data.signal;

for i=1:numOfNotches
    sig = filtfilt(comb_b{i},comb_a{i},sig);
end

notched = struct('signal',sig,'SamplingFrequency',fs,'TimeVector',tV,'FrequencyVector',fV,'FilePath',FileDir);
% Author: Jack Tchimino
%%
function filtered = filterHF(data,fc)


signal = data.signal;
timeVector = data.TimeVector;
frequencyVector = data.FrequencyVector;
fs = data.SamplingFrequency;
FileDir = data.FilePath;

WHP = fc/(0.5*fs);

[b,a] = butter(5,WHP,'low');
% [b,a] = cheby2(2,0.1,WHP,'high');
% [b,a] = ellip(2,1,20,WHP,'high');

% figure
% freqz(b,a,frequencyVector,fs)

filtered = filtfilt(b,a,signal);
filtered = struct('signal',filtered,'SamplingFrequency',fs,'TimeVector',timeVector,'FrequencyVector',frequencyVector,'FilePath',FileDir);

% Author: Jack Tchimino
%%
function detrended = removeDC(data)

fHP = 0.1;

signal = data.signal;
timeVector = data.TimeVector;
frequencyVector = data.FrequencyVector;
fs = data.SamplingFrequency;

WHP = fHP/(0.5*fs);

[b,a] = butter(2,WHP,'high');
% [b,a] = cheby2(2,0.1,WHP,'high');
% [b,a] = ellip(2,1,20,WHP,'high');

% figure
% freqz(b,a,frequencyVector,fs)

filtered = filter(b,a,signal);
detrended = struct('signal',filtered,'SamplingFrequency',fs,'TimeVector',timeVector,'FrequencyVector',frequencyVector);

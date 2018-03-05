%% Author: Jack Tchimino

% This function resamples the signal if the coeff variable is other than
% one. Generates and returns the resampled signal, new sampling frequency
% and time and frequency vectors.

function output = prepareSignal(signal,fs,coeff)

if coeff<1
    disp("Error: Oversampling not possible (prepareSignal)");
    return;
end

% coeff is the resampling coefficient
fs = ceil(fs/coeff);

signal = signal(1:coeff:end);

N = length(signal);     % number of samples

time = N/fs;

timeVector = 0:time/N:time;
if isequal(length(timeVector),N)==0
    timeVector = 0:time/N:time-time/N;
end

frequencyVector = 0:1/time:fs;
if isequal(length(frequencyVector),N)==0
    frequencyVector = 0:1/time:fs-1/time;
end

output = struct('signal',double(signal),'SamplingFrequency',fs,'TimeVector',timeVector,'FrequencyVector',frequencyVector);
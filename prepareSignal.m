%% Author: Jack Tchimino

% This function resamples the signal if the coeff variable is other than
% one. Generates and returns the resampled signal, new sampling frequency
% and time and frequency vectors.

% INPUTS:
% signal: 1xN time series data array
% fs    : data sampling frequency
% coeff : resampling coefficient
% filepath : directory of source file, can be omitted

function output = prepareSignal(signal,fs,coeff,varargin)

if coeff<1
    error("Error: Oversampling not possible");
%     return;
end

% coeff is the resampling coefficient
fs = ceil(fs/coeff);

signal = signal(1:coeff:end);

N = length(signal);     % number of samples

time = N/fs; % total recording duration

timeVector = 0:time/N:time;
if isequal(length(timeVector),N)==0
    timeVector = 0:time/N:time-time/N;
end

frequencyVector = 0:1/time:fs;
if isequal(length(frequencyVector),N)==0
    frequencyVector = 0:1/time:fs-1/time;
end

if isempty(varargin)
    FilePath = 'N/A';
else
    FilePath = varargin{1};
end

output = struct('signal',double(signal),'SamplingFrequency',fs,'TimeVector',timeVector,'FrequencyVector',frequencyVector,'FilePath',FilePath);
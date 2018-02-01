%% Author: Jack Tchimino

%Run everything

%%load signal
load('C:\Users\Jack Tchimino\Dropbox\BME\internship\data\signal\14053-05 - 5.mat')
signalHJ = RAW.ewave;
fs = RAW.eHz;

%% get shorter signal

seconds=10;
% numOfSamples = 300*fs;
% tenMinSignal = signal(1:numOfSamples);
% tV = [300/numOfSamples:300/numOfSamples:300];
% fV = [1/300:1/300:fs];
numOfSamples = ceil(seconds*fs);
shortSignal = signalHJ(1:numOfSamples+1);

prepped = prepareSignal(shortSignal,fs,1);

plot(prepped.TimeVector,prepped.signal)

spikes = detectSpikes(prepped.signal,prepped.SamplingFrequency,25,25,0.25,'down');

% create averaging function
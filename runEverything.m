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

%% Spike Detection

spikes = detectSpikes(prepped.signal,prepped.SamplingFrequency,25,25,0.25,'down');

%% Averaging

averaged = averaging(spikes,'Y');

%% Spectrogram/fourier of the whole waveform/of each spike
% not sure if spectrogram of each spike is useful:
%   too few samples
%   can get better results with simple fft
% spectrogram better at monitoring transient changes in longer signal(?)
spectrogram(prepped.signal,6000,2000,'power')
view(-45,65)
% use detrend in spectrogram to remove DC offset and spike at 0Hz
spectrogram(detrend(prepped.signal),600,400,600,fs,'yaxis')

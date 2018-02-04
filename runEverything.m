%% Author: Jack Tchimino

%Run everything

%% load cerebellar signal
clear all
clc
prompt = 'input 1 for cerebellar, 2 for epilepsy: \n';
selection = input(prompt);

if selection==1
% load cerebellar signal
    load('C:\Users\Jack Tchimino\Dropbox\BME\internship\data\signal\14053-05 - 5.mat')
    signal = RAW.ewave;
    fs = RAW.eHz;
    clear RAW S SET TRIG
elseif selection==2
% load epilepsy signal
    load('C:\Users\Jack Tchimino\Dropbox\BME\Thesis\data\seizure data\2011_08_24_TH_THAL_3ch_0008.mat')
    signal = block.segments{1, 1}.analogsignals{1, 1}.signal;
    fs = block.segments{1, 1}.analogsignals{1, 1}.sampling_rate;
    seizureTimes = block.segments{1, 1}.events{1, 1}.times;
else
    disp("Error: Try again");
    return;
end
clear prompt selection
clc
%% get shorter signal

seconds=20;
% numOfSamples = 300*fs;
% tenMinSignal = signal(1:numOfSamples);
% tV = [300/numOfSamples:300/numOfSamples:300];
% fV = [1/300:1/300:fs];
numOfSamples = ceil(seconds*fs);
shortSignal = signal(1:numOfSamples+1);

prepped = prepareSignal(shortSignal,fs,1);

%plot(prepped.TimeVector,prepped.signal)

%% Spike Detection

spikes = detectSpikes(prepped.signal,prepped.SamplingFrequency,25,25,0.25,'down');

%% Averaging

averaged = averaging(spikes,'N');

%% Spectrogram/fourier of the whole waveform/of each spike
% not sure if spectrogram of each spike is useful:
%   too few samples
%   can get better results with simple fft
% spectrogram better at monitoring transient changes in longer signal(?)

% use detrend in spectrogram to remove DC offset and spike at 0Hz
spectrogram(detrend(prepped.signal),ceil(0.5*fs),ceil(0.25*fs),ceil(0.5*fs),fs,'yaxis')
view(-45,65)

%% adaptive segmentation
% must define: -reference window
%              -test window
%              -dissimilarity measure

segments = adaptiveSegmentation(prepped.signal,prepped.SamplingFrequency,0.5,100,5,'Y','SEM');






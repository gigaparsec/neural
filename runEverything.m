%% Author: Jack Tchimino

%Run everything

%% load cerebellar signal
clear all
clc
prompt = 'input 1 for cerebellar, 2 for epilepsy, 3 for 16/32 channel ephys: \n';
selection = input(prompt);

if selection==1
% load cerebellar signal
    load('C:\Users\Jack Tchimino\Dropbox\BME\internship\data\signal\14053-05 - 5.mat')
    raw = RAW.ewave;
    fs = RAW.eHz;
    clear RAW S SET TRIG
elseif selection==2
% load epilepsy signal
    load('C:\Users\Jack Tchimino\Dropbox\BME\Thesis\data\seizure data\2011_08_24_TH_THAL_3ch_0008.mat')
    raw = block.segments{1, 1}.analogsignals{1, 1}.signal;
    fs = block.segments{1, 1}.analogsignals{1, 1}.sampling_rate;
    seizureTimes = block.segments{1, 1}.events{1, 1}.times;
elseif selection==3
    [raw,timestamps,info] = load_open_ephys_data_faster("C:\Users\Jack Tchimino\Documents\17-MI22288-01f_VTA_2018-02-08_16-20-21\100_CH2_2.continuous");
    fs = info.header.sampleRate;
    clear timestamps info
else
    disp("Error: Try again");
    return;
end
clear prompt selection
clc
%% get shorter signal

seconds=10;
% numOfSamples = 300*fs;
% tenMinSignal = signal(1:numOfSamples);
% tV = [300/numOfSamples:300/numOfSamples:300];
% fV = [1/300:1/300:fs];
numOfSamples = ceil(seconds*fs);
shortSignal = raw(1:numOfSamples+1);

signal = prepareSignal(shortSignal,fs,1);

plotted = signal.signal;
plot(signal.TimeVector,plotted)
hold on
plot(signal.TimeVector(1:5:end),plotted(1:5:end))

%% Remove DC offset

signal = removeDC(signal);

%% Power spectrum plot

makePowerSpectrum(signal);

%% Notch filters



%% Spike Detection

spikes = detectSpikes(signal.signal,signal.SamplingFrequency,25,25,0.25,'down');

%% Averaging

averaged = averaging(spikes,'N');

figure
subplot(2,2,[1, 2])
plot(signal.TimeVector,signal.signal);
title('Cerebellar Signal: spike train')
xlabel('Time [sec]');ylabel('Voltage [mV]');
subplot(2,2,3)
plot(spikes{1}.TimeVector,spikes{1}.signal);
hold all
for i=2:length(spikes)
    plot(spikes{i}.TimeVector,spikes{i}.signal);
end
title('Aligned Spikes')
xlabel('Time [sec]');ylabel('Voltage [mV]');
xlim([0,averaged.TimeVector(end)])
subplot(2,2,4)
plot(averaged.TimeVector,averaged.signal,'LineWidth',2);
title('Averaged Clear Spike')
xlabel('Time [sec]');ylabel('Voltage [mV]');
xlim([0,averaged.TimeVector(end)])
suptitle('Averaging Demonstration')

%% Spectrogram/fourier of the whole waveform/of each spike
% not sure if spectrogram of each spike is useful:
%   too few samples
%   can get better results with simple fft
% spectrogram better at monitoring transient changes in longer signal(?)

% use detrend in spectrogram to remove DC offset and spike at 0Hz
spectrogram(detrend(signal.signal),ceil(0.5*fs),ceil(0.25*fs),ceil(0.5*fs),fs,'yaxis')
view(-45,65)

%% adaptive segmentation
% must define: -reference window
%              -test window
%              -dissimilarity measure

segments = adaptiveSegmentation(signal.signal,signal.SamplingFrequency,0.5,100,5,'Y','SEM');






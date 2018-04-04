%% Author: Jack Tchimino

%Run everything

%% load cerebellar signal
% clear all
clc
prompt = 'input 1 for cerebellar, 2 for epilepsy, 3 for 16/32 channel ephys, 4 for ensemble 32-channel: \n';
selection = input(prompt);

if selection==1
% load cerebellar signal
    load('C:\Users\Jack Tchimino\Dropbox\BME\internship\data\signal\14053-05 - 5.mat')
    raw = RAW.ewave;
    fs = RAW.eHz;
    clear RAW S SET TRIG
elseif selection==2
% load epilepsy signal
    load('C:\Users\Jack Tchimino\Dropbox\BME\Thesis\data\seizure data\2011_08_24_TH_THAL_3ch_0011.mat') %11-> 7 seizures
    raw = block.segments{1, 1}.analogsignals{1, 1}.signal;
    fs = block.segments{1, 1}.analogsignals{1, 1}.sampling_rate;
    seizureTimes = block.segments{1, 1}.events{1, 1}.times;
elseif selection==3
    [a,b] = uigetfile({'*.*',  'All Files (*.*)'},'Pick recording file','C:\Users\Jack Tchimino\Documents\Recordings');
    [raw,timestamps,info] = load_open_ephys_data_faster([b,a]);
%     [raw,timestamps,info] = load_open_ephys_data_faster("C:\Users\Jack Tchimino\Documents\17-MI22288-01f_VTA_2018-02-08_16-20-21\100_CH2_2.continuous");
%     [raw,timestamps,info] = load_open_ephys_data_faster("");
    fs = info.header.sampleRate;
    clear timestamps info c
    
elseif selection==4
    prompt2 = 'how many seconds? \n';
    RecTime = input(prompt2);
    Channels = importEnsemble(RecTime);
    clear RecTime prompt2 
else
    disp("Error: Try again");
    return;
end
clear prompt selection
% clc
%%
% get shorter signal

seconds=100;
% numOfSamples = 300*fs;
% tenMinSignal = signal(1:numOfSamples);
% tV = [300/numOfSamples:300/numOfSamples:300];
% fV = [1/300:1/300:fs];
numOfSamples = ceil(seconds*fs);
shortSignal = raw(1:numOfSamples+1);

OriginalSignal = prepareSignal(shortSignal,fs,1);
signal1 = prepareSignal(shortSignal,fs,1);
signal2 = prepareSignal(shortSignal,fs,4);
signal3 = prepareSignal(shortSignal,fs,10);
signal4 = prepareSignal(shortSignal,fs,30);

figure(2)
plotted = OriginalSignal.signal;
plot(OriginalSignal.TimeVector,plotted)
ax = gca;
ax.FontSize=14;
xlabel('Time [sec]');ylabel('Voltage [mV]')
title('\fontsize{20}Prefrontal Cortex')

% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% drawnow;
% set(get(handle(gcf),'JavaFrame'),'Maximized',1);
% y1 = get(gca,'ylim');
% x1 = get(gca,'xlim');
% hold all
% plot([50 50],y1);
% plot(x1,[2 2])
% 
% plot(signal.TimeVector(1:5:end),plotted(1:5:end))

%% Remove DC offset DO THIS ALWAYS
signal = OriginalSignal;
signal = removeDC(signal);

%% lowpass

fsignal = filterHF(signal,1000);
figure(1)
plot(signal.TimeVector,OriginalSignal.signal)
hold on
plot(fsignal.TimeVector,fsignal.signal,'LineWidth',2)

%% Power spectrum plot

makePowerSpectrum(signal);

%% Notch filters
clc
% fn = [100 200 300 400 500 600 700];% 1500 2000 2300 2700 2900 3500 4000 4500 6000 6500 7500 7600 7700 7800 7900 8000 8100 8200 8400 8500 11000 11500 13000 13400 13500 13600 13700 13800 13900 14000 14100 14200 14300 14400 14999];
fn=50;
clean = notch(signal,fn,50,'N');
%% plot notch powerspectrums
figure
subplot(2,1,1)
h1=makePowerSpectrum(signal);
xlim([10^-2 fs/2])
ax = gca;
for i=1:length(fn)
h2=line([fn(i) fn(i)],[ax.YLim(1) ax.YLim(2)]);
h2.Color='r';
h2.LineWidth=2;
end
ax.FontSize = 14;
for i=1:length(fn) uistack(h1); end
title('Original')
subplot(2,1,2)
h1 = makePowerSpectrum(clean);
xlim([10^-2 fs/2])
ax = gca;
for i=1:length(fn)
h2=line([fn(i) fn(i)],[ax.YLim(1) ax.YLim(2)]);
h2.Color='r';
h2.LineWidth=2;
end
ax.FontSize = 14;
for i=1:length(fn) uistack(h1); end
title('Notched')
suptitle('\fontsize{20}Power Spectrums of Original and Notched Signal')
%% plot notched
% close all

figure
subplot(2,1,1)
plot(OriginalSignal.TimeVector,OriginalSignal.signal)
xlabel('Time [sec]');ylabel('Voltage [mV]');
ax=gca;
ax.FontSize = 14;
subplot(2,1,2)
plot(clean.TimeVector,clean.signal)
xlabel('Time [sec]');ylabel('Voltage [mV]');
ax=gca;
ax.FontSize = 14;
suptitle('\fontsize{20}Original and notched signal')

figure
plot(OriginalSignal.TimeVector,OriginalSignal.signal)
hold all
plot(signal.TimeVector,fsignal.signal)
plot(clean.TimeVector,clean.signal)
xlabel('Time [sec]');ylabel('Voltage [mV]');
ax=gca;
ax.FontSize = 14;
legend('Original','LPF','LPF & Notched');
title('\fontsize{20}Original, filtered and notched signal')
%% Spike Detection

spikes = detectSpikes(signal.signal,signal.SamplingFrequency,25,25,0.25,'down');

%% Averaging

averaged = averaging(spikes,'N');

figure
subplot(2,2,[1, 2])
plot(signal.TimeVector,signal.signal);
title('\fontsize{16}Cerebellar Signal: spike train')
xlabel('\fontsize{14}Time [sec]');ylabel('\fontsize{14}Voltage [mV]');
subplot(2,2,3)
plot(spikes{1}.TimeVector,spikes{1}.signal);
hold all
for i=2:length(spikes)
    plot(spikes{i}.TimeVector,spikes{i}.signal);
end
title('\fontsize{16}Aligned Spikes')
xlabel('\fontsize{14}Time [sec]');ylabel('\fontsize{14}Voltage [mV]');
xlim([0,averaged.TimeVector(end)])
subplot(2,2,4)
plot(averaged.TimeVector,averaged.signal,'LineWidth',2);
title('\fontsize{16}Averaged Clear Spike')
xlabel('\fontsize{14}Time [sec]');ylabel('\fontsize{14}Voltage [mV]');
xlim([0,averaged.TimeVector(end)])
suptitle('\fontsize{20}Averaging Demonstration')

%% Spectrogram/fourier of the whole waveform/of each spike
% not sure if spectrogram of each spike is useful:
%   too few samples
%   can get better results with simple fft
% spectrogram better at monitoring transient changes in longer signal(?)

% use detrend in spectrogram to remove DC offset and spike at 0Hz
% spectrogram(detrend(signal.signal),ceil(0.5*fs),ceil(0.25*fs),ceil(0.5*fs),fs,'yaxis')
close all
% 
temp = seizureTimes/60; % seizure times in minutes

figure
spectrogram(detrend(clean.signal),ceil(1*fs),ceil(0.9*fs),ceil(1*fs),fs,'yaxis')
view(-45,65)
ax=gca;
ax.FontSize = 14;
for i=1:length(temp)
h=line([temp(i), temp(i)],[0 0],[-100,0]);
h.Color='r';
h.LineWidth=2;
end
title('\fontsize{20}Spectrogram of epileptic neural signal notched at 50Hz')
% figure
% spectrogram(detrend(OriginalSignal.signal),ceil(0.01*fs),ceil(0.005*fs),ceil(1*fs),fs,'yaxis')
% view(-45,65)
% ax=gca;
% ax.FontSize = 14;
% for i=1:length(temp)
% h=line([temp(i), temp(i)],[0 0],[-100,0]);
% h.Color='r';
% h.LineWidth=2;
% end
% suptitle('\fontsize{20}Spectrogram of epileptic neural signal')

%% adaptive segmentation
% must define: -reference window
%              -test window
%              -dissimilarity measure
clc
segments = adaptiveSegmentation(clean,0.5,100,4,'Y','SEM');

%% Segment and spectrogram

fixedSegments = fixedSegmentation(clean,30);
numOfSegs = length(fixedSegments);

WINDOW = ceil(0.01*fs);
NOVERLAP = ceil(0.009);
NFFT = ceil(0.01*fs);

for i=1:numOfSegs
set(0,'DefaultFigureWindowStyle','docked')
figure(i)
spectrogram(detrend(fixedSegments{i}.signal),WINDOW,NOVERLAP,NFFT,fs,'yaxis');
% xticks(fixedSegments{i}.TimeVector);
view(-45,65)
ax=gca;
ax.FontSize = 14;
title(['Segment ',num2str(i)])
end






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
%     [raw,timestamps,info] = load_open_ephys_data_faster([b,a]);
    [raw,~,info] = load_open_ephys_data_faster([b,a]);
%     [raw,timestamps,info] = load_open_ephys_data_faster("C:\Users\Jack Tchimino\Documents\17-MI22288-01f_VTA_2018-02-08_16-20-21\100_CH2_2.continuous");
%     [raw,timestamps,info] = load_open_ephys_data_faster("");
    fs = info.header.sampleRate;
    clear timestamps info c
    
elseif selection==4
    prompt2 = 'how many seconds? \n';
    RecTime = input(prompt2);
    [Channels, Aux, ADC] = importEnsemble(RecTime);
    clear RecTime prompt2 
else
    disp("Error: Try again");
    return;
end
clear prompt %selection
% clc

% sig30k = prepareSignal(raw,fs,1);
% sig6k = prepareSignal(raw,fs,5);
% sig3k = prepareSignal(raw,fs,10);
original = prepareSignal(raw,fs,1,[b,a]);
makePowerSpectrum(original,'x');
startTime = 6.458e06;
endTime = 1.155e07;
endTime = 6.458e06+160*fs;
original2 = prepareSignal(original.signal(startTime:endTime),fs,1,original.FilePath);

figure
plot(original.TimeVector,original.signal)
figure
plot(sig1k.TimeVector,sig1k.signal)

makePowerSpectrum(original,'x')
makePowerSpectrum(sig1k,'x')
original161036_2 = prepareSignal(raw,fs,1,[b,a]);
sTime = 1.445e7;
eTime = 1.654e7;
original3 = prepareSignal(w.*original.signal(sTime:sTime+50*fs),fs,1,[b,a]);
%%
% get shorter signal

seconds=10;
% numOfSamples = 300*fs;
% tenMinSignal = signal(1:numOfSamples);
% tV = [300/numOfSamples:300/numOfSamples:300];
% fV = [1/300:1/300:fs];
numOfSamples = ceil(seconds*fs);
shortSignal = raw(1:numOfSamples+1);

OriginalSignal = prepareSignal(shortSignal,fs,1);
% signal1 = prepareSignal(shortSignal,fs,1);
% signal2 = prepareSignal(shortSignal,fs,4);
% signal3 = prepareSignal(shortSignal,fs,10);
% signal4 = prepareSignal(shortSignal,fs,30);

% figure(2)
% plotted = OriginalSignal.signal;
% plot(OriginalSignal.TimeVector,plotted)
% ax = gca;
% ax.FontSize=14;
% xlabel('Time [sec]');ylabel('Voltage [mV]')
% title('\fontsize{20}Prefrontal Cortex')

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

%%
sig30kW = removeDC(sig30k,0.5);
sig6kW = removeDC(sig6k,0.5);
sig3kW = removeDC(sig3k,0.5);
sig1kW = removeDC(sig1k,0.5);

original = removeDC(original,0.5);
%% Remove DC offset DO THIS ALWAYS

if selection ~= 4
signal = OriginalSignal;
signal = removeDC(signal,0.5);
end

if selection ==4
    for i=1:length(Channels)
         Channels{i} = removeDC(Channels{i},0.5);
    end
end


%% lowpass
fc_LowPass = 200;
if selection ==4
    for i=1:length(Channels)
         Channels{i} = filterHF(Channels{i},fc_LowPass);
    end
    disp('end')
else
fsignal = filterHF(signal,fc_LowPass);
% figure(1)
% plot(signal.TimeVector,OriginalSignal.signal)
% hold on
% plot(fsignal.TimeVector,fsignal.signal,'LineWidth',2)
end

fsig30k = filterHF(sig30kW,fc_LowPass);
fsig6k = filterHF(sig6kW,fc_LowPass);
fsig3k = filterHF(sig3kW,fc_LowPass);
fsig1k = filterHF(sig1kW,fc_LowPass);

foriginal = filterHF(original,fc_LowPass);

makePowerSpectrum(foriginal,'x');
foriginal_sub = prepareSignal(foriginal.signal,foriginal.SamplingFrequency,30);
makePowerSpectrum(foriginal_sub,'x')

%% Power spectrum plot
% signal = Channels{1};
figure
makePowerSpectrum(sig30kW);
title([b,a,' raw signal FFT, highPassed at 0.5Hz'])
figure
makePowerSpectrum(sig1kW)
title([b,a,'   subsampled at 1kHz FFT, highPassed at 0.5Hz'])

figure
makePowerSpectrum(fsig30k);
title([b,a,' raw signal FFT, highPassed at 0.5Hz, lowpassed at 200Hz'])
figure
makePowerSpectrum(fsig1k)
title([b,a,' subsampled at 1kHz FFT, highPassed at 0.5Hz, lowpassed at 200Hz'])

%% Notch filters
clc
% fn = [100 200 300 400 500 600 700];% 1500 2000 2300 2700 2900 3500 4000 4500 6000 6500 7500 7600 7700 7800 7900 8000 8100 8200 8400 8500 11000 11500 13000 13400 13500 13600 13700 13800 13900 14000 14100 14200 14300 14400 14999];
fn=[50 100 150 200 250 300 350 400];
nOriginal = notch(foriginal,fn,50,'N');
original2 = prepareSignal(nOriginal.signal(startTime:endTime),fs,1,[b,a]);%,original.FilePath);

spectrogram(detrend(original2.signal),ceil(1*fs),ceil(0.9*fs),ceil(1*fs),fs,'yaxis')

spectrogram(detrend(original2.signal(1:30:end)),ceil(1*fs/30),ceil(0.9*fs/30),ceil(1*fs/30),fs/30,'yaxis')

original2_sub = prepareSignal(original2.signal,fs,30,[b,a]);
subfs = original2_sub.SamplingFrequency;
spectrogram((original2_sub.signal),ceil(1*subfs),ceil(0.9*subfs),ceil(1*subfs),subfs,'yaxis')

n30k = notch(fsig30k,fn,50,'N');
n6k = notch(fsig6k,fn,50,'N');
n3k = notch(fsig3k,fn,50,'N');
n1k = notch(fsig1k,fn,50,'N');

plot(sig30k.TimeVector,sig30k.signal)
hold all
plot(n30k.TimeVector,n30k.signal,'LineWidth',2)
% plot(n6k.TimeVector,n6k.signal)
% plot(n3k.TimeVector,n3k.signal)
plot(n1k.TimeVector,n1k.signal,'LineWidth',2)
% legend('Original','filtered (0.5-200) notched @ 100,200,300','subsampled @6k','subsampled @3k','subsampled @1k')
legend('Original','filtered (0.5-200) notched @ 100,200,300','subsampled @1k')


Hn30k = removeDC(n30k,4);
Hn1k  = removeDC(n1k,4);
plot(sig30k.TimeVector,sig30k.signal)
hold all
plot(n30k.TimeVector,Hn30k.signal,'LineWidth',2)
% plot(n6k.TimeVector,n6k.signal)
% plot(n3k.TimeVector,n3k.signal)
plot(n1k.TimeVector,Hn1k.signal,'LineWidth',2)
% legend('Original','filtered (0.5-200) notched @ 100,200,300','subsampled @6k','subsampled @3k','subsampled @1k')
legend('Original','filtered (0.5-200) notched @ 100,200,300','subsampled @1k')
title('highpassed at 4Hz')

between0and4 = filterBP(sig1k,0.1,4);
between2and4 = filterBP(sig1k,2,4);
between0and10 = filterBP(sig1k,0.1,10);
between4and10 = filterBP(sig1k,4,10);





%%
temp = removeDC(fsignal,4);
plot(OriginalSignal.TimeVector,OriginalSignal.signal)
hold on
plot(fsignal.TimeVector,temp.signal)

makePowerSpectrum(OriginalSignal)
figure
makePowerSpectrum(removeDC(fsignal,4))

%% fixedSegmentation

segmented = fixedSegmentation(signal,1);


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
h1 = makePowerSpectrum(signal);
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
plot(signal.TimeVector,signal.signal)
xlabel('Time [sec]');ylabel('Voltage [mV]');
ax=gca;
ax.FontSize = 14;
suptitle('\fontsize{20}Original and notched signal')

figure
plot(OriginalSignal.TimeVector,OriginalSignal.signal)
hold all
plot(signal.TimeVector,fsignal.signal)
plot(signal.TimeVector,signal.signal)
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
spectrogram(detrend(signal.signal),ceil(1*fs),ceil(0.9*fs),ceil(1*fs),fs,'yaxis')
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
segments = adaptiveSegmentation(signal,0.5,100,4,'Y','SEM');

%% Segment and spectrogram

fixedSegments = fixedSegmentation(signal,30);
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






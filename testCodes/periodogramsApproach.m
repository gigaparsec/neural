% N = length(relNums);
% 
% % for k=1:N
% %     sum=0;
% %     for n=1:N-k
% %         sum = sum+(relNums(n+k));
% %     end
% %     rx(k) = sum/N;
% % end
% 
% sigF = flip(sig);
% 
% rx = conv(sig,sigF)/N;
% 
% S = (((abs(fft(rx(1:end))))))/N;
% figure
% plot(S(1:ceil(length(S)/2)))
% title('custom')
% 
% pxx = periodogram(sig);
% figure
% plot(pxx)
% title('matlab')

%%
epochsRaw = cell(1,35);

for i=1:length(stims)
    epochsRaw{i} = raw(stims(i)-5*fs:stims(i)+25*fs);
end

epochs = cell(1,35);
epochsSub = cell(1,35);


fn = 50:50:500;
for i=1:35
    temp = prepareSignal(epochsRaw{i},fs,1,[b,a]);
    temp = customNotch(temp,fn,100,'N');
    tempN = customNotch(temp,[2.433,2.333],10,'N');
    tempFN = filterHF(tempN,300);
    epochsSub{i}=prepareSignal(tempFN.signal,tempFN.SamplingFrequency,30,tempFN.FilePath);
    epochs{i}=prepareSignal(tempFN.signal,tempFN.SamplingFrequency,1,tempFN.FilePath);
end

sum = zeros(1,length(epochs{1}.signal));

for i=1:35
    sig = epochs{i}.signal;
    sum = sum + sig;
end

sum = (epochs{1}.signal+epochs{2}.signal+epochs{3}.signal+epochs{4}.signal+epochs{5}.signal...
    +epochs{5}.signal+epochs{6}.signal+epochs{7}.signal+epochs{8}.signal+epochs{9}.signal...
    +epochs{10}.signal+epochs{11}.signal+epochs{12}.signal+epochs{13}.signal+epochs{14}.signal...
    +epochs{15}.signal+epochs{16}.signal+epochs{17}.signal+epochs{18}.signal+epochs{19}.signal...
    +epochs{20}.signal+epochs{21}.signal+epochs{22}.signal+epochs{23}.signal+epochs{23}.signal...
    +epochs{25}.signal+epochs{26}.signal+epochs{27}.signal+epochs{28}.signal+epochs{29}.signal...
    +epochs{30}.signal+epochs{31}.signal+epochs{32}.signal+epochs{33}.signal+epochs{34}.signal...
    +epochs{35}.signal)/35;

plot(epochs{1}.TimeVector,sum)


plot(epochs{1}.TimeVector,epochs{1}.signal)

makePowerSpectrum(epochs{1},'x')

test = customNotch(epochs{1},[2.433,2.333],10,'N');

makePowerSpectrum(test,'x')
plot(test.TimeVector,test.signal)

%%
sig = epochs{1}.signal;
fsam = epochs{1}.SamplingFrequency;
[pxx2,freqV] = periodogram(sig,hamming(length(sig)),ceil(length(sig)/2),fsam);
[pxx2p,freqVp] = periodogram(sig,hamming(length(sig)),ceil(length(sig)/2),fsam,'power');

plot(freqV,10*log10(pxx2))
hold on
plot(freqVp,pxx2p)


%% compare periodograms subsampled

fs = epochsSub{1}.SamplingFrequency;
fp = epochsSub{1}.FilePath;

stimDur = 5*fs:6*fs;
noStimDur = 10*fs:11*fs;

stimEpoch = prepareSignal(epochsSub{1}.signal(stimDur),fs,1,fp);

noStimEpoch = prepareSignal(epochsSub{1}.signal(noStimDur),fs,1,fp);

sigS = stimEpoch.signal;
sigNS = noStimEpoch.signal;
winS = hamming(length(sigS));
lenS = ceil(length(sigS)/2);
winNS = hamming(length(sigNS));
lenNS = ceil(length(sigNS)/2);


[pxxStim,freqStim] = periodogram(sigS,winS,lenS,fs);
[pxxNStim,freqNStim] = periodogram(sigNS,winNS,lenNS,fs);


[pxxStimW,freqStimW] = pwelch(sigS,winS,lenS,length(winS),fs);
[pxxNStimW,freqNStimW] = pwelch(sigNS,winNS,lenNS,length(winS),fs);

figure
plot(freqStim,smooth(10*log10(pxxStim)));
hold all
plot(freqNStim,smooth(10*log10(pxxNStim)));
plot(freqStimW,smooth(10*log10(pxxStimW)));
plot(freqNStimW,smooth(10*log10(pxxNStimW)));
legend('periodogram stim','periodogram not stim','welch periodogram stim','welch periodogram not stim')
title('subsampled')
%% compare periodograms not subsampled

fs = epochs{1}.SamplingFrequency;
fp = epochs{1}.FilePath;

stimDur = 5*fs:6*fs;
noStimDur = 10*fs:11*fs;

stimEpoch = prepareSignal(epochs{1}.signal(stimDur),fs,1,fp);

noStimEpoch = prepareSignal(epochs{1}.signal(noStimDur),fs,1,fp);

sigS = stimEpoch.signal;
sigNS = noStimEpoch.signal;
winS = hamming(length(sigS));
lenS = ceil(length(sigS)/2);
winNS = hamming(length(sigNS));
lenNS = ceil(length(sigNS)/2);


[pxxStim,freqStim] = periodogram(sigS,winS,lenS,fs);
[pxxNStim,freqNStim] = periodogram(sigNS,winNS,lenNS,fs);


[pxxStimW,freqStimW] = pwelch(sigS,winS,lenS,length(winS),fs);
[pxxNStimW,freqNStimW] = pwelch(sigNS,winNS,lenNS,length(winS),fs);
figure
plot(freqStim,smooth(10*log10(pxxStim)));
hold all
plot(freqNStim,smooth(10*log10(pxxNStim)));
plot(freqStimW,smooth(10*log10(pxxStimW)));
plot(freqNStimW,smooth(10*log10(pxxNStimW)));
legend('periodogram stim','periodogram not stim','welch periodogram stim','welch periodogram not stim')
title('original')

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
% PERIODOGRAM COMPARISON
%% 
clc
clear all
load stims.mat;
importSignals
epochsRaw = cell(1,35);

for i=1:length(stims)
    epochsRaw{i} = raw(stims(i)-5*fs:stims(i)+25*fs);
end

epochs = cell(1,35);
epochsSub = cell(1,35);


fn = 50:50:500;
for i=1:35
    temp = prepareSignal(epochsRaw{i},30000,1,[b,a]);
    temp2 = customNotch(temp,fn,250,'N');
%     tempN = customNotch(temp,[2.433,2.333],10,'N');
    tempFN = filterHF(temp2,300);
    epochsSub{i}=prepareSignal(tempFN.signal,tempFN.SamplingFrequency,2,tempFN.FilePath);
    epochs{i}=prepareSignal(tempFN.signal,tempFN.SamplingFrequency,2,tempFN.FilePath);
end

fs = epochsSub{1}.SamplingFrequency;
fp = epochsSub{1}.FilePath;

stimDur = 4*fs:8*fs;
noStimDur = 10*fs:14*fs;

PXXStimcell = cell(1,35);
PXXNoStimcell = cell(1,35);

sumPXX = 0;
sumNPXX = 0;

for i=1:35
    stimEpoch = prepareSignal(epochsSub{i}.signal(stimDur),fs,1,fp);
    noStimEpoch = prepareSignal(epochsSub{i}.signal(noStimDur),fs,1,fp);
    sigS = stimEpoch.signal;
    sigNS = noStimEpoch.signal;
    winS = hamming(length(sigS));
    lenS = ceil(length(sigS)/2);
    winNS = hamming(length(sigNS));
    lenNS = ceil(length(sigNS)/2);
    
    [pxxStim,freqStim] = periodogram(sigS,winS,lenS,fs);
    [pxxNStim,freqNStim] = periodogram(sigNS,winNS,lenNS,fs);
    
    sumPXX = sumPXX + pxxStim;
    sumNPXX = sumNPXX + pxxNStim;
    
    PXXStimcell{i} = pxxStim;
    PXXNoStimcell{i} = pxxNStim;
end

clear sigS sigNS winS winNS lenS lenNS

avgPXX = 10*log10(sumPXX)/35;
avgPXXN = 10*log10(sumNPXX)/35;

plot(avgPXXN)
hold on
plot(avgPXX)

% figure
% plot(freqStim,smooth(10*log10(PXXStimcell{1})))
% hold on
% plot(freqStim,smooth(10*log10(PXXNoStimcell{1})))
% legend('stim','no stim');
% title(a)

%% ====================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
% FFT COMPARISON
clc
clear all
load stims.mat;
importSignals
epochsRaw = cell(1,35);

for i=1:length(stims)
    epochsRaw{i} = raw(stims(i)-5*fs:stims(i)+25*fs);
end

epochs = cell(1,35);
epochsSub = cell(1,35);


fn = 50:50:10000;
for i=1:35
    temp = prepareSignal(epochsRaw{i},30000,1,[b,a]);
    temp2 = customNotch(temp,fn,250,'N');
%     tempN = customNotch(temp,[2.433,2.333],10,'N');
    tempFN = filterHF(temp2,10000);
    epochsSub{i}=prepareSignal(tempFN.signal,tempFN.SamplingFrequency,2,tempFN.FilePath);
    epochs{i}=prepareSignal(tempFN.signal,tempFN.SamplingFrequency,2,tempFN.FilePath);
end

fs = epochsSub{1}.SamplingFrequency;
fp = epochsSub{1}.FilePath;

stimDur = 4*fs:8*fs;
noStimDur = 10*fs:14*fs;

PXXStimcell = cell(1,35);
PXXNoStimcell = cell(1,35);

for i=1:35
stimEpoch = prepareSignal(epochsSub{i}.signal(stimDur),fs,1,fp);
noStimEpoch = prepareSignal(epochsSub{i}.signal(noStimDur),fs,1,fp);
h = figure(1);
makePowerSpectrum(stimEpoch,'x');
hold on
makePowerSpectrum(noStimEpoch,'x');
legend('stim','nostim')
pause
close(h)
end

%%
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
% SEGMENTATION

clc
clear all
load stims.mat;
importSignals
epochsRaw = cell(1,35);

for i=1:length(stims)
    epochsRaw{i} = raw(stims(i)-5*fs:stims(i)+25*fs);
end

epochs = cell(1,35);
epochsSub = cell(1,35);


fn = 50:50:500;
for i=1:35
    temp = prepareSignal(epochsRaw{i},30000,1,[b,a]);
    temp2 = customNotch(temp,fn,250,'N');
%     tempN = customNotch(temp,[2.433,2.333],10,'N');
    tempFN = filterHF(temp2,300);
%     epochsSub{i}=prepareSignal(tempFN.signal,tempFN.SamplingFrequency,30,tempFN.FilePath);
    epochs6khz{i}=prepareSignal(tempFN.signal,tempFN.SamplingFrequency,5,tempFN.FilePath);
end



disp('Running...........................')
for k=1:5
tic
[segmentBounds,SEM] = adaptiveSegm(epochsSub{k},1000,1000,'periodogram','fixed');
toc

SB{k} = segmentBounds;
SEMc{k} = SEM;

figure(1)
plot(epochsSub{k}.TimeVector,SEM)
hold all


end
disp('end')


%%
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
% PCA/SVD
clc
clear all
load stims.mat;
importSignals
epochsRaw = cell(1,35);

for i=1:length(stims)
    epochsRaw{i} = raw(stims(i)-5*fs:stims(i)+25*fs);
end

epochs = cell(1,35);
epochsSub = cell(1,35);


fn = 50:50:500;
for i=1:35
    temp = prepareSignal(epochsRaw{i},30000,1,[b,a]);
    temp2 = customNotch(temp,fn,250,'N');
%     tempN = customNotch(temp,[2.433,2.333],10,'N');
    tempFN = filterHF(temp2,300);
    epochsSub{i}=prepareSignal(tempFN.signal,tempFN.SamplingFrequency,30,tempFN.FilePath);
%     epochs6khz{i}=prepareSignal(tempFN.signal,tempFN.SamplingFrequency,5,tempFN.FilePath);
end


response = zeros(length(epochsSub{i}.signal),35);

for i=1:35
    response(:,i) = epochsSub{i}.signal';
end

[U,S,V] = svdecon(response);

recon = U*S*V';

% get one component

for i=1:35
    Stemp = zeros(size(S));

    Stemp(i,i) = S(i,i);

    partial_reconstr{i} =  U*Stemp*V';
   
end

figure
for i=1:18
subplot(6,3,i)
plot(partial_reconstr{3}(:,i))
title(['eigenvalue',num2str(i)]);
end

%%
load('stims.mat');

FStim = stims(1);

firstEpochs = cell(1,16);
stimEpochs = cell(1,16);
baseEpochs = cell(1,16);

fn = 50:50:1000;
L = 1000*3+1;
matrStim = zeros(L,16);
matrBase = zeros(L,16);

for i=1:16
    directory = ['C:\Users\Jack Tchimino\Documents\Recordings\2018-05-09_19-29-57\100_CH',num2str(i),'.continuous'];
    [raw,~,info] = load_open_ephys_data_faster(directory);
    fs = info.header.sampleRate;
    
    stimEpoch = raw(FStim-fs:FStim+2*fs);
    temp = prepareSignal(stimEpoch,fs,1,directory);
    tempF = filterHF(temp,300);
    tempFN = customNotch(tempF,fn,250,'N');
    stimEpochs{i} = prepareSignal(tempFN.signal,tempFN.SamplingFrequency,30,directory);
    matrStim(:,i) = stimEpochs{i}.signal';
    
    baseEpoch = raw(FStim+fs*20:FStim+23*fs);
    temp = prepareSignal(baseEpoch,fs,1,directory);
    tempF = filterHF(temp,300);
    tempFN = customNotch(tempF,fn,250,'N');    
    baseEpochs{i} = prepareSignal(tempFN.signal,tempFN.SamplingFrequency,30,directory);
    matrBase(:,i) = baseEpochs{i}.signal';
end

% [coeff,score,latent] = pca(matrStim);
% [coeff,score,latent] = pca(matrBase);

muS = mean(matrStim);
muB = mean(matrBase);
%%
% close(h)
[coeffS,scoreS,latentS] = pca(matrStim,'numComponents',10,'Algorithm','svd');
[coeffB,scoreB,latentB] = pca(matrBase,'numComponents',10,'Algorithm','svd');
%%
eiSt = 3;
eiEnd = 3;


xS = scoreS(:,eiSt:eiEnd)*coeffS(:,eiSt:eiEnd)'+muS;
xB = scoreB(:,eiSt:eiEnd)*coeffB(:,eiSt:eiEnd)'+muB;
h=figure(1);
plot(baseEpochs{1}.TimeVector,xB(:,1))
hold on
% plot(stimEpochs{1}.TimeVector,firstEpochs{1}.signal)
plot(stimEpochs{1}.TimeVector,xS(:,1))
% 
% test = prepareSignal(x(:,1),1000,1,'test');
% figure
% makePowerSpectrum(test,'x')
% title([num2str(eiSt),'-',num2str(eiEnd)])

% semilogx(firstEpochs{1}.FrequencyVector,abs(fft(x(:,1))))


%%

fs = firstEpochs{1}.SamplingFrequency;

figure(2)
spectrogram(detrend(firstEpochs{12}.signal),ceil(0.01*fs),ceil(0.001*fs),ceil(0.1*fs),fs,'yaxis')

%% PCA change detection

load('stims.mat');

fn = 50:50:1000;

range = (stims(1)-30000*5):(stims(10)+30000*5-1);

signals = zeros(length(range)/30,16);

disp('Running..............')
for i=1:16
    disp(num2str(i))
    directory = ['C:\Users\Jack Tchimino\Documents\Recordings\2018-05-09_19-29-57\100_CH',num2str(i),'.continuous'];
    [raw,~,info] = load_open_ephys_data_faster(directory);
    fs = info.header.sampleRate; 
%     sig = prepareSignal(raw(range),fs,1,directory);
%     
%     temp = prepareSignal(raw(range),fs,1,directory);
%     tempF = filterHF(temp,300);
%     tempFN = customNotch(tempF,fn,250,'N');
%     
%     sig = prepareSignal(tempFN.signal,fs,30,'N/A');
    
    sig = preprocessing300(raw(range),fs,fn,directory);

    signals(:,i) = sig.signal;
    
end
disp('end')

%%
sigSF = sig.SamplingFrequency;
sigTV = sig.TimeVector;
sigFV = sig.FrequencyVector;
sigL  = length(sig.signal);

pcaWindowSize = sigSF/10;

start=pcaWindowSize/2+1;

iRange = start:(sigL-start);

lambdasMatrix = zeros(16,length(iRange));
k=0;
disp('Running PCA loop.....')
for i=start:(sigL-start)
    pcaWindowRange = (i-0.5*pcaWindowSize):(i+0.5*pcaWindowSize);
    pcaWindowSig = signals(pcaWindowRange,:);
    [~,~,latent] = pca(pcaWindowSig,'Algorithm','svd');
%     k=k+1;
    lambdasMatrix(:,i) = latent;
end
disp('End of PCA loop')

plot(sigTV(1:length(lambdasMatrix)),log10(lambdasMatrix))
ax=gca;
% ax.FontSize = 14;
for i=5:30:30*9+5
h=line([(i),(i)],ax.YLim);
h.Color='r';
h.LineWidth=2;
uistack(h,'bottom');
end

%%

% clear all

fifties = 50:100:6550;
fn1 = [fifties+3.13 , fifties-3.13];
fn2 = 50:50:7000;
fn3 = 30:30:6900;

fn = [fn1, fn2, fn3, 615, 1230];

epochs = cell(1,16);

% signals = zeros(30001,16);

for i=1:16
    directory = ['C:\Users\Jack Tchimino\Documents\Recordings\2018-05-09_16-56-09\100_CH',num2str(i),'.continuous'];
    [raw,~,info] = load_open_ephys_data_faster(directory);
    fs = info.header.sampleRate;
    raw = raw(85.1*30000:86.6*30000);
    sig = prepareSignal(raw,fs,1,directory);
    epochs{i} = preprocessSignal(sig,10,7000,fn,2);
    signals(:,i)=epochs{i}.signal;
end

figure
plot(epochs{1}.TimeVector,epochs{1}.signal)
ylim([-400 400])
% grid on
xlabel('Time [sec]');
ylabel('Voltage [mV]');
title('Spontaneous Firing of Dopaminergic Neurons');

% hold on
% for k=2:16
% plot(epochs{k}.signal)
% end
% 
% [U,S,V] = svdecon(signals);
% 
% recon = U*S*V';
% 
% plot(recon(:,1))
% hold on
% plot(signals(:,1))
% 
% for i=1:16
%     Stemp = zeros(size(S));
%     
%     Stemp(i,i) = S(i,i);
%     
%     partialRecon{i} = U*Stemp*V';
% end
% 
% figure
% for i=1:16
%     subplot(4,4,i)
%     plot(partialRecon{i}(:,1));
%     ylim([-100 100])
%     title(['Channel 1, component ',num2str(i)])
% end

%%

sig = epochs{1};

sigSF = sig.SamplingFrequency;
sigTV = sig.TimeVector;
sigFV = sig.FrequencyVector;
sigL  = length(sig.signal);


pcaWindowSize = sigSF/30;

start=(pcaWindowSize/2+1);

iRange = start:(sigL-start);

lambdasMatrix = zeros(16,length(iRange));
k=0;
disp('Running PCA loop.....')
for i=start:(sigL-start)
    pcaWindowRange = (i-0.5*pcaWindowSize):(i+0.5*pcaWindowSize);
    pcaWindowSig = signals(pcaWindowRange,:);
    [~,~,latent] = pca(pcaWindowSig,'Algorithm','svd');
%     k=k+1;
    lambdasMatrix(:,i) = latent;
end
disp('End of PCA loop')
% plot(sigTV(1:length(lambdasMatrix)),log10(lambdasMatrix))

plot(sigTV(1:length(lambdasMatrix)),log10(lambdasMatrix(1,:)))
title('Magnitude of First Eigenvalue')
xlabel('Time [sec]');
ylabel('Log-Magnitude');










































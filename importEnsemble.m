% Author: Jack Tchimino

% IMPORT ENSEMBLE

function [Channels, Aux, ADC]= importEnsemble(seconds)
folder = uigetdir('C:\Users\Jack Tchimino\Documents\Recordings');

Channels = cell(32,1);
tic
for i=1:1:32
    recPath = [folder,'\100_CH',num2str(i),'.continuous'];
    [raw,~,info] = load_open_ephys_data_faster(recPath);
    fs = info.header.sampleRate;
    clear info
    numOfSamples = ceil(seconds*fs);
    shortSignal = raw(1:numOfSamples);
    Channels{i} = prepareSignal(shortSignal,fs,1);
    clear raw
end


Aux = cell(3,1);
for i=1:3
    recPath = [folder,'\100_AUX',num2str(i),'.continuous']; 
    [raw,~,info] = load_open_ephys_data_faster(recPath);
    fs = info.header.sampleRate;
    clear info
    numOfSamples = ceil(seconds*fs);
    shortSignal = raw(1:numOfSamples);
    Aux{i} = prepareSignal(shortSignal,fs,1);
    clear raw
end

ADC = cell(8,1);
for i=1:1:8
    recPath = [folder,'\100_ADC',num2str(i),'.continuous']; 
    [raw,~,info] = load_open_ephys_data_faster(recPath);
    fs = info.header.sampleRate;
    clear info
    numOfSamples = ceil(seconds*fs);
    shortSignal = raw(1:numOfSamples);
    ADC{i} = prepareSignal(shortSignal,fs,1);
    clear raw
end
toc
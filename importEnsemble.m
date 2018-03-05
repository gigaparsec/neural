% Author: Jack Tchimino

% IMPORT ENSEMBLE

function Channels = importEnsemble(seconds)
folder = uigetdir('C:\Users\Jack Tchimino\Documents\Recordings');

Channels = cell(32,1);
tic
for i=1:1:32
    recPath = [folder,'\100_CH',num2str(i),'.continuous'];
    [raw,timestamps,info] = load_open_ephys_data_faster(recPath);
    fs = info.header.sampleRate;
    clear timestamps info
    numOfSamples = ceil(seconds*fs);
    shortSignal = raw(1:numOfSamples);
    Channels{i} = prepareSignal(shortSignal,fs,1);
    clear raw
end
toc
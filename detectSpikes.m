%% Author: Jack Tchimino

%% Peak Detection and segmentation
% This code receives a raw signal as input, detects peaks based on a specific threshold
% aligns the peaks based on predetermined limits and returns an array
% containing all the spikes. upOrDown specifies whether the detected
% spikes will be upward or downward, can take values 'up' or 'down'.

function spikes = detectSpikes(signal,fs,samplesBefore,samplesAfter,threshold,upOrDown)
flag=0;
if strcmp(upOrDown,'down')==1
    signal = -signal;
    flag=1;
end
len = length(signal);
[peak,loc] = findpeaks(signal,'MinPeakHeight',threshold); % loc contains the indices of the most prominent peaks
%maybe put minpeakheight as an input variable

% simple averaging (Ragnayyan, Biomedical signal Analysis p.143-4)
% the peak is regarded as a reference point
N = samplesBefore+samplesAfter+1;
numOfSpikes = length(loc);
spikesMat = zeros(numOfSpikes,N); % each peak has been observed to be ~50 samples in length

for i=1:1:numOfSpikes % first and last spikes can be omitted, in case the window exceeds matrix dimensions
    if loc(i)<samplesBefore
        continue
    end
    if loc(i)>len-samplesAfter
        break
    end
    spikesMat(i,:) = signal(loc(i)-samplesBefore:loc(i)+samplesAfter);
end

spikesMat(~any(spikesMat,2),:) = [];  % delete all nonzero rows

if flag==1
    spikesMat=-spikesMat;
end

spikes = cell(size(spikesMat,1),1);

for i=1:size(spikesMat,1)
    spikes{i} = prepareSignal(spikesMat(i,:),fs,1);
end

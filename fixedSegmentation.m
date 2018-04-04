% Author: Jack Tchimino

% FIXED SEGMENTATION: this script will segment the waveform in segments of
% fixed length, starting from the beginning. The script produces only as
% many full durations fit in the signal. the goal is the better
% visualization of spectrograms, because a spectrogram of the entire length
% of the signal would be unintelligible

function fixedSegments = fixedSegmentation(data,duration)

signal = data.signal;
fs = data.SamplingFrequency;
fV = data.FrequencyVector;
tV = data.TimeVector;

rem = mod(tV(end),duration);

numOfSegs = tV(end-rem)/duration;

durationSamples = duration*fs;

fixedSegments = cell(numOfSegs,1);


segFV = 0:1/duration:fs;
if isequal(length(segFV),durationSamples)==0
    segFV = 0:1/duration:fs-1/duration;
end

for i=0:numOfSegs-1
    beginning = 1+i*durationSamples;
    endofSeg  = beginning+durationSamples-1;
    segment = signal(beginning:endofSeg);
    segTV   = tV(beginning:endofSeg);
    fixedSegments{i+1} = struct('signal',double(segment),'SamplingFrequency',fs,'TimeVector',segTV,'FrequencyVector',segFV);
end
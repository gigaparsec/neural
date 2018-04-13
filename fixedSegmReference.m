% Author: Jack Tchimino

%% Fixed Segmentation around predetermined reference points
% refPoints is a 1D array, given by the findStimulations script

function segments = fixedSegmReference(data, samplesBefore, samplesAfter, refPoints)

    signal = data.signal;
    fs = data.SamplingFrequency;
    tV = data.TimeVector;

    % check if the last segment exceeds the length of the remaining signal,
    % if so, decrease the number of segments by 1
    if refPoints(end)>length(signal)-samplesAfter
        numOfSegs = length(refPoints)-1;
    else
        numOfSegs = length(refPoints);
    end

    %check if the first segment will begin earlier than the first sample of
    %the signal. if so, start the iteration from the second segment
    if refPoints(1)<=samplesBefore
        iterationStart=2;
    else
        iterationStart=1;
    end
    
    % preallocate for the segments
    segments = cell(numOfSegs,1);
    
    for i=iterationStart:1:numOfSegs
        TDsegment = signal(refPoints(i)-samplesBefore:refPoints(i)+samplesAfter);
        segments{i} = createSignalStruct(TDsegment,fs,tV(refPoints(i)-samplesBefore:refPoints(i)+samplesAfter));
    end
    
    disp(['Successful segmentation. Number of segments: ',num2str(length(segments))]);
    
end
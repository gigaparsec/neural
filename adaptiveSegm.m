% Author: Jack Tchimino
%% Adaptive Segmentation algorithms as they are implemented in Sornmo's book
% data: signal struct
% threshold: SEM threshold value for segmentation
% N: window length in samples
% method: 'periodogram' or 'whitening'
function [segmentBoundaries,SEM] = adaptiveSegm(data,threshold,N,method)

if strcmp(method,'periodogram')==1
    signal = data.signal;
    fs = data.SamplingFrequency;
    tV = data.TimeVector;
    
    % SORNMO FORMULA 3.194 PERIODOGRAM APPROACH
    startPoint = 1;
    movingPoint = 2;
    
    numOfSamples = length(signal);
    SEM = zeros(1,numOfSamples);
    segmentNum = 1;
    segmentBoundaries = zeros(1,numOfSamples);
    i=1;
    while movingPoint+N<numOfSamples
        fixedWindow = signal(startPoint:startPoint+N);
        movingWindow = signal(movingPoint:movingPoint+N);
        ACFfixed = xcorr(fixedWindow);% MIGHT HAVE TO USE AUTOCORR - MUST FIND OUT THE DIFFERENCE
        ACFmoving = xcorr(movingWindow);
        num = sum((ACFmoving-ACFfixed).^2);
        den = ACFfixed(ceil(end/2))*ACFmoving(ceil(end/2));
        Delta1 = num/den;
        SEM(i)=Delta1;
        i=i+1;
        if Delta1<threshold
            movingPoint = movingPoint+1;
        elseif Delta1>=threshold
            segmentBoundaries(segmentNum)=movingPoint;
            segmentNum = segmentNum + 1;
            startPoint = movingPoint;
        end
    end
    segmentBoundaries = nonzeros(segmentBoundaries);
    
elseif strcmp(method,'whitening')==1
    
    
    
    % whitening approach is sensitive to increase in power. must first
    % check the spectrograms to make sure it's worth looking into. the
    % power might decrease during DA release
    
    
    
    
else
    disp('Error, try again. "periodogram" or "whitening" in method');
    segmentBoundaries=[];
    SEM=[];
    return;
end
end
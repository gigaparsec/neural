% Author: Jack Tchimino
%% Adaptive Segmentation algorithms as they are implemented in Sornmo's book
% data: signal struct
% threshold: SEM threshold value for segmentation
% N: window length in samples
% method: 'periodogram' or 'whitening'
function [segmentBoundaries,SEM] = adaptiveSegm(data,threshold,N,method,refWinMode)
    signal = data.signal;
    fs = data.SamplingFrequency;
    tV = data.TimeVector;
if strcmp(method,'periodogram')==1

    if strcmp(refWinMode,'fixed')==1
    % SORNMO FORMULA 3.194 PERIODOGRAM APPROACH, fixed reference window
    startPoint = 1;
    movingPoint = 2;
    
    numOfSamples = length(signal);
    SEM = zeros(1,numOfSamples);
    segmentNum = 1;
    segmentBoundaries = zeros(1,numOfSamples);
    i=1;
    while movingPoint+N<numOfSamples
        referenceWindow = signal(startPoint:startPoint+N);
        movingWindow = signal(movingPoint:movingPoint+N);
        ACFfixed = xcorr(referenceWindow);% MIGHT HAVE TO USE AUTOCORR - MUST FIND OUT THE DIFFERENCE
        ACFmoving = xcorr(movingWindow);
        num = sum((ACFmoving-ACFfixed).^2);
        den = ACFfixed(ceil(end/2))*ACFmoving(ceil(end/2));
        Delta1 = num/den;
%         SEM(i)=Delta1;
%         i=i+1;
%         if Delta1<threshold
%             movingPoint = movingPoint+1;
%         elseif Delta1>=threshold
%             segmentBoundaries(segmentNum)=movingPoint;
%             segmentNum = segmentNum + 1;
%             startPoint = movingPoint;
%         end

        if Delta1<threshold
            movingPoint = movingPoint+1;
            SEM(i)=Delta1;
            i=i+1;
        elseif Delta1>=threshold
            segmentBoundaries(segmentNum)=movingPoint;
            segmentNum = segmentNum + 1;
            startPoint = movingPoint;
            i=i+N;
            SEM(i)=Delta1;
            i=i+1;
        end
    end
    segmentBoundaries = nonzeros(segmentBoundaries);
    
    elseif strcmp(refWinMode,'growing')
        numOfSamples = length(signal);
        SEM = zeros(1,numOfSamples);
        segmentNum=1;
        segmentBoundaries = zeros(1,numOfSamples);
        i=1;
        startPoint = 1;
        movingPoint = N+1;
        referenceAdder = 0;
        while movingPoint+N<numOfSamples
            referenceWindow = signal(startPoint:startPoint+N+referenceAdder);
            movingWindow = signal(movingPoint:movingPoint+N);
            ACFfixed = xcorr(referenceWindow,N);% MIGHT HAVE TO USE AUTOCORR - MUST FIND OUT THE DIFFERENCE
            ACFmoving = xcorr(movingWindow,N); % decleare the lags in this cse, because the two vectors will be of different lengths in due time
            num = sum((ACFmoving-ACFfixed).^2);
            den = ACFfixed(ceil(end/2))*ACFmoving(ceil(end/2));
            Delta1 = num/den;
            SEM(i)=Delta1;
            i=i+1;
            if Delta1<threshold
                movingPoint = movingPoint+1;
                referenceAdder = referenceAdder + 1;
            elseif Delta1>threshold
                segmentBoundaries(segmentNum)=movingPoint; 
                segmentNum = segmentNum + 1;
                startPoint = movingPoint; % -N???
                referenceAdder = 0;
                movingPoint = movingPoint + N; % and don't change?
            end
        end  
        segmentBoundaries = nonzeros(segmentBoundaries);
    else
        error('Error. Input "growing" or "fixed" reference window length')
%         segmentBoundaries=[];
%         SEM=[];
%         return;
    end
elseif strcmp(method,'whitening')==1
    

    
    
    % whitening approach is only sensitive to increase in power. must first
    % check the spectrograms to make sure it's worth looking into. the
    % power might decrease during DA release
    segmentBoundaries=[];
    SEM=[];
    return;    
    
else
    error('Error, try again. "periodogram" or "whitening" in method');
%     segmentBoundaries=[];
%     SEM=[];
%     return;
end
end
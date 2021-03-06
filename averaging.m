%% Author: Jack Tchimino

% Basic Averaging Function. Takes a "spikes" variable and averages over
% time. plotYN can be omitted, but if it is Y, the function will also plot
% all the spikes in one figure and their average in another.

function averaged = averaging(spikes,plotYN)

N = size(spikes,1); %number of spikes

sampleNum = size(spikes{1}.signal,2);

summed = zeros(1,sampleNum);

for i=1:size(spikes,1)
    summed = summed+spikes{i}.signal;
end

averageSpike = summed/N;

if nargin==2 && strcmp(plotYN,'Y')==1
    figure
    plot(spikes{i}.TimeVector,spikes{i}.signal)
    hold all
    for i=2:size(spikes,1)
    plot(spikes{i}.TimeVector,spikes{i}.signal)
    end
    xlabel('Time [sec]');ylabel('Amplitude [V]');
    xlim([0,spikes{i}.TimeVector(end)]);
    title('Spikes Detected in the Signal');
    
    figure
    plot(spikes{1}.TimeVector,averageSpike)
    xlabel('Time [sec]');ylabel('Amplitude [V]');
    xlim([0,spikes{i}.TimeVector(end)]);
    title('Average of All Spikes'); 
end

averaged = struct('signal',averageSpike,'SamplingFrequency',spikes{1}.SamplingFrequency,'TimeVector',spikes{1}.TimeVector,'FrequencyVector',spikes{1}.FrequencyVector);
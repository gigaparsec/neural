% Author: Jack Tchimino

% IMPORT STIMULATION TIMINGS

function stimulationTimes = importStimulation(seconds)
folder = uigetdir('C:\Users\Jack Tchimino\Documents\Recordings');

[a,b] = uigetfile({'*.*',  'All Files (*.*)'},'Pick stimulation file','C:\Users\Jack Tchimino\Documents\Recordings');
[raw,timestamps,info] = load_open_ephys_data_faster([b,a]);

fs = info.header.sampleRate;
clear timestamps info
numOfSamples = ceil(seconds*fs);
shortSignal = raw(1:numOfSamples);

stimulationTimes = prepareSignal(shortSignal,fs,1);
% MUST FIND A WAY TO DETECT THE STIMULATION ONSET 
% find, isolate the 1st element, count 1-2 seconds, change elements to
% zero, find again...
% 0000011000001100000...
%find = 6 7 13 14
% first = 6
% new = 000000000000110000.. reiterate (?)
clear raw
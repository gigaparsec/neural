clc
prompt = 'input 1 for cerebellar, 2 for epilepsy, 3 for 16/32 channel ephys, 4 for ensemble 32-channel: \n';
selection = input(prompt);

if selection==1
% load cerebellar signal
    load('C:\Users\Jack Tchimino\Dropbox\BME\internship\data\signal\14053-05 - 5.mat')
    raw = RAW.ewave;
    fs = RAW.eHz;
    clear RAW S SET TRIG
elseif selection==2
% load epilepsy signal
    load('C:\Users\Jack Tchimino\Dropbox\BME\Thesis\data\seizure data\2011_08_24_TH_THAL_3ch_0011.mat') %11-> 7 seizures
    raw = block.segments{1, 1}.analogsignals{1, 1}.signal;
    fs = block.segments{1, 1}.analogsignals{1, 1}.sampling_rate;
    seizureTimes = block.segments{1, 1}.events{1, 1}.times;
elseif selection==3
    [a,b] = uigetfile({'*.*',  'All Files (*.*)'},'Pick recording file','C:\Users\Jack Tchimino\Documents\Recordings');
%     [raw,timestamps,info] = load_open_ephys_data_faster([b,a]);
    [raw,timestamps,info] = load_open_ephys_data_faster([b,a]);
%     [raw,timestamps,info] = load_open_ephys_data_faster("C:\Users\Jack Tchimino\Documents\17-MI22288-01f_VTA_2018-02-08_16-20-21\100_CH2_2.continuous");
%     [raw,timestamps,info] = load_open_ephys_data_faster("");
    fs = info.header.sampleRate;
    clear  info c
    
elseif selection==4
    prompt2 = 'how many seconds? \n';
    RecTime = input(prompt2);
    [Channels, Aux, ADC] = importEnsemble(RecTime);
    clear RecTime prompt2 
else
    disp("Error: Try again");
    return;
end
clear prompt %selection
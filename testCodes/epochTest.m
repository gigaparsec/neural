%%
importSignals

% figure
% subplot(2,1,1)
% plot(timestamps,sign)
% subplot(2,1,2)
% plot(timestamps,raw)

% figure
subplot(2,1,1)
% plot(original.TimeVector(6.534e06:1.169e07),sign(6.534e06:1.169e07))
plot(original.TimeVector,sign)
title('CH1')
subplot(2,1,2)
% plot(original.TimeVector(6.534e06:1.169e07),raw(6.534e06:1.169e07))
plot(original.TimeVector,raw)
% ylim([min(raw(6.534e06:1.169e07))-0.1,max(raw(6.534e06:1.169e07))+0.1])
title(a)

%%

stim = [stim1.DataIndex,stim2.DataIndex,stim3.DataIndex,stim4.DataIndex];

%%

% responses = cell(32*length(stim),1);
responsesFiltAt5Hz = cell(32,length(stim));
% k=1;
for i=1:16 % only 16 workable signals
    Directory = ['C:\Users\Jack Tchimino\Documents\Recordings\2018-04-19_14-55-45\100_CH',num2str(i),'.continuous'];
    [raw,~,info] = load_open_ephys_data_faster(Directory);
    fs = info.header.sampleRate;
    clear info;
    original = prepareSignal(raw,fs,1,Directory); % create original
%     original = removeDC(original,0.5);  %remove DC
    fc_high=10;
    originalFiltAt5Hz = removeDC(original,fc_high);  %remove DC
    fc_low = 300;
%     filteredOriginal = filterHF(original,fc_low); %Lowpass filter
    filteredOriginal5Hz = filterHF(originalFiltAt5Hz,fc_low);
    fn = [100 200 300 333 400];
%     notchedOriginal = notch(filteredOriginal,fn,50,'N');
    notchedOriginal5Hz = notch(filteredOriginal5Hz,fn,50,'N');
    tempFn = [46.9 53.1];
    fn2 = [tempFn,50+tempFn,100+tempFn,150+tempFn,200+tempFn,250+tempFn,300+tempFn,350+tempFn,400+tempFn,450+tempFn];
    notchedOriginal5Hz = notch(notchedOriginal5Hz,fn2,100,'N'); % notch at the harmonics around the multiples of 50
    for j=1:length(stim)
%         tempSignal = notchedOriginal.signal(stim(j)-original.SamplingFrequency*5:stim(j)+original.SamplingFrequency*15);
%         responses{k} = prepareSignal(tempSignal,fs,30,[original.FilePath]);%subsample the epoch and save to responses array
        
        tempSignal = notchedOriginal5Hz.signal(stim(j)-original.SamplingFrequency*5:stim(j)+original.SamplingFrequency*15);
        responsesFiltAt5Hz{i,j} = prepareSignal(tempSignal,fs,30,[original.FilePath]);%subsample the epoch and save to responses array
        
%         k=k+1;
    responsesFiltAt5Hz{i,j}.Description = ['HPF @ ',num2str(fc_high),'Hz, LPF @ ',num2str(fc_low),'Hz, Notched at [',num2str(fn2),' ',num2str(fn),'Hz'];
    end
    
    disp(['File ',num2str(i),' parsed'])
end
% responsesClean = cell(size(responses));
% for i=1:length(responses)
%     responsesClean{i} = notch(responses{i},333,50,'N');
% end
%%
figure
makePowerSpectrum(responsesFiltAt5Hz{1,4},'x');

% figure
% plot(responses{1}.TimeVector,responses{1}.signal)
% hold all
% for i=2:length(responses)
%     plot(responses{i}.TimeVector,responses{i}.signal)
% end

figure
hold all
for i=1:4:length(responses)
    plot(responses{i}.TimeVector,responses{i}.signal)
end
title('1st stimulation responses')

figure
hold all
for i=2:4:length(responses)
    plot(responses{i}.TimeVector,responses{i}.signal)
end
title('2nd stimulation responses')

figure
hold all
for i=3:4:length(responses)
    plot(responses{i}.TimeVector,responses{i}.signal)
end
title('3rd stimulation responses')

figure
hold all
for i=1:32
    plot(responses_2018_04_19_14_55_45{i,4}.TimeVector,responses_2018_04_19_14_55_45{i,4}.signal)
end
title('4th stimulation responses')

figure
makePowerSpectrum(original,'x')
% spectrogram(responses{1}.signal,)

figure
plot(responses{4}.TimeVector,responses{4}.signal)
hold on
plot(responsesFiltAt5Hz{4}.TimeVector,responsesFiltAt5Hz{4}.signal)

%%
figure
subplot(2,1,1)
plot(responsesFiltAt5Hz{4}.TimeVector,responsesFiltAt5Hz{4}.signal)
subplot(2,1,2)
spectrogram(responsesFiltAt5Hz{4}.signal,10,9,200,1000,'yaxis')

%% file 2018_05_09_19-29-57 ADC7 -> stimulations

% first_stim = cursor_info.DataIndex;
first_stim = 1160501;
stims = zeros(1,35);

stims(1) = first_stim;

for i=2:length(stims)
    stims(i) = 30*fs+stims(i-1);
end

responses_2018_05_09_19_29_57 = cell(35,length(stims));
%%
disp('Running...')
for i=1:17
%     i
    % import recording files
    if i<17
    Directory = ['C:\Users\Jack Tchimino\Documents\Recordings\2018-05-09_19-29-57\100_CH',num2str(i),'.continuous'];
    [raw,~,info] = load_open_ephys_data_faster(Directory);
    fs = info.header.sampleRate;
    clear info;
    original = prepareSignal(raw,fs,1,Directory); % create original
    fc_high=10;
    originalHPF = removeDC(original,fc_high);  %remove DC
    fc_low = 300;
    BPfiltered = filterHF(originalHPF,fc_low); %Lowpass filter
    % notch
    fn = [30 50 60 90 100 120 ...
          150 180 200 210 240 ...
          250 270 300 330 350 ...
          360 390 400 420 450 ...
          480 500 510 540 550 ...
          580 600 630 650 660 690 700];
      notched = notch(BPfiltered,fn,100,'N');
      
      for j = 1:length(stims)
%           j
          response = notched.signal(stims(j)-fs*5:stims(j)+fs*20);
          responses_2018_05_09_19_29_57{i,j} = prepareSignal(response,fs,15,original.FilePath);
          responses_2018_05_09_19_29_57{i,j}.Description = ['HPF @ ',num2str(fc_high),...
                                       'Hz, LPF @ ',num2str(fc_low),'Hz, Notched at ['...
                                       ,num2str(fn),'Hz'];
      end
      disp(['File ',num2str(i),' parsed'])
    else
        Directory = 'C:\Users\Jack Tchimino\Documents\Recordings\2018-05-09_19-29-57\100_ADC7.continuous';
        [raw,~,info] = load_open_ephys_data_faster(Directory);
        fs = info.header.sampleRate;
        clear info;
%         stimulation = prepareSignal(raw,fs,1,Directory); % create original
        for j = 1:length(stims)
            stimulation_pulse = raw(stims(j)-fs*5:stims(j)+fs*20);
            responses_2018_05_09_19_29_57{i,j} = stimulation_pulse(1:15:end);%prepareSignal(stimulation_pulse,fs,15,Directory);
        end
    end
end


%%

figure
subplot(3,3,1)
plot(responses_2018_05_09_19_29_57{1,1}.TimeVector,responses_2018_05_09_19_29_57{1,1}.signal)
subplot(3,3,2)
plot(responses_2018_05_09_19_29_57{1,2}.TimeVector,responses_2018_05_09_19_29_57{1,2}.signal)
subplot(3,3,3)
plot(responses_2018_05_09_19_29_57{1,3}.TimeVector,responses_2018_05_09_19_29_57{1,3}.signal)

subplot(3,3,4)
plot(responses_2018_05_09_19_29_57{2,1}.TimeVector,responses_2018_05_09_19_29_57{2,1}.signal)
subplot(3,3,5)
plot(responses_2018_05_09_19_29_57{2,2}.TimeVector,responses_2018_05_09_19_29_57{2,2}.signal)
subplot(3,3,6)
plot(responses_2018_05_09_19_29_57{2,3}.TimeVector,responses_2018_05_09_19_29_57{2,3}.signal)

subplot(3,3,7)
plot(responses_2018_05_09_19_29_57{3,1}.TimeVector,responses_2018_05_09_19_29_57{17,1})
subplot(3,3,8)
plot(responses_2018_05_09_19_29_57{3,2}.TimeVector,responses_2018_05_09_19_29_57{17,2})
subplot(3,3,9)
plot(responses_2018_05_09_19_29_57{3,3}.TimeVector,responses_2018_05_09_19_29_57{17,3})


%%

figure

subplot(3,3,1)
spectrogram(responses_2018_05_09_19_29_57{1,1}.signal,100,5,100,1000,'yaxis')
subplot(3,3,2)
plot(responses_2018_05_09_19_29_57{1,2}.TimeVector,responses_2018_05_09_19_29_57{1,2}.signal)
subplot(3,3,3)
plot(responses_2018_05_09_19_29_57{1,3}.TimeVector,responses_2018_05_09_19_29_57{1,3}.signal)

subplot(3,3,4)
plot(responses_2018_05_09_19_29_57{2,1}.TimeVector,responses_2018_05_09_19_29_57{2,1}.signal)
subplot(3,3,5)
plot(responses_2018_05_09_19_29_57{2,2}.TimeVector,responses_2018_05_09_19_29_57{2,2}.signal)
subplot(3,3,6)
plot(responses_2018_05_09_19_29_57{2,3}.TimeVector,responses_2018_05_09_19_29_57{2,3}.signal)


subplot(3,3,7)
plot(responses_2018_05_09_19_29_57{3,1}.TimeVector,responses_2018_05_09_19_29_57{17,1})
subplot(3,3,8)
plot(responses_2018_05_09_19_29_57{3,2}.TimeVector,responses_2018_05_09_19_29_57{17,2})
subplot(3,3,9)
plot(responses_2018_05_09_19_29_57{3,3}.TimeVector,responses_2018_05_09_19_29_57{17,3})












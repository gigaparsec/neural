function preprocessedSignal = preprocessSignal(data,HighPassCutoff,LowPassCutoff,fNotch,subsamplingFactor)
% set(0,'DefaultFigureVisible','off');
signal = data.signal;
% tV = data.TimeVector;
fV = data.FrequencyVector;
fs = data.SamplingFrequency;
Directory = data.FilePath;
Q = 70;
if HighPassCutoff>=LowPassCutoff
    error('High Pass filter cutoff cannot be higher than Low Pass Cutoff');
end
if subsamplingFactor<1 || subsamplingFactor>30000
    error('Invalid subsampling factor');
end

WLP = LowPassCutoff/(0.5*fs);
[b,a] = butter(5,WLP,'low');

% h=figure(1000);
% subplot(2,2,1)
% [h,f]=freqz(b,a,fV,fs);
% plot(f,20*log10(abs(h)),'LineWidth',2)
% title('LowPass Butterworth')
% xlabel('Freqency [Hz]');ylabel('Magnitude [dB]');
% xlim([1,LowPassCutoff*2])

preprocessed = filtfilt(b,a,signal);

if HighPassCutoff>0
    WHP = HighPassCutoff/(0.5*fs);
    [b,a] = butter(4,WHP,'high');
    preprocessed = filtfilt(b,a,preprocessed);
%     h=figure(1000);
%     subplot(2,2,2)
%     [h,f]=freqz(b,a,fV,fs);
%     plot(f,20*log10(abs(h)),'LineWidth',2)
%     title('HighPass Butterworth')
%     xlabel('Freqency [Hz]');ylabel('Magnitude [dB]');
%     xlim([0.1,HighPassCutoff*2])
end



if ~isempty(fNotch)
    nyquist = fs/2;
    w0 = fNotch/nyquist;
    bw = w0/Q;
    for i=1:length(fNotch)
        [b,a] = iirnotch(w0(i),bw(i));%,20);
        preprocessed = filtfilt(b,a,preprocessed);
%         h=figure(1000);
%         subplot(2,2,[3,4])
%         [h,f]=freqz(b,a,fV,fs);
%         plot(f,20*log10(abs(h)),'LineWidth',2,'color',[0    0.4470    0.7410])
%         hold on
%         title('Comb Filter')
%         xlabel('Freqency [Hz]');ylabel('Magnitude [dB]');
    end
%     xlim([1,fNotch(end)*2])
end

% h=figure(1000);
% suptitle('Preprocessing Filters')


preprocessedSignal = prepareSignal(preprocessed,fs,subsamplingFactor,Directory);
% close(h)
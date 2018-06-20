function cleanSignal = customNotch(data,freq,Q,plot)
% close all

fV = data.FrequencyVector;
tV = data.TimeVector;
fs = data.SamplingFrequency;
y = data.signal;
FileDir = data.FilePath;

nyquist = fs/2;
w0 = freq/nyquist;
% Q=5;
bw = w0/Q;
% y=signal;
% if strcmp(plot,'Y')==1
% for i=1:length(freq)
% [b,a] = iirnotch(w0(i),bw(i));%,20);
% figure(1)
% freqz(b,a,fV,fs)
% hold all
% y = filtfilt(b,a,y);
% end
% else
    for i=1:length(freq)
        tic
    [b,a] = iirnotch(w0(i),bw(i));%,20);
    toc
%     figure(1)
%     freqz(b,a,fV,fs)
%     hold all
tic
    y = filtfilt(b,a,y);
    toc
    end
% end    
    
    
cleanSignal = struct('signal',y,'SamplingFrequency',fs,'TimeVector',tV,'FrequencyVector',fV,'FilePath',FileDir);

% figure
% plot(t,signal);hold on
% plot(t,signal);
% % % % % 
% % % % % signal = data.signal;
% % % % % fs = data.SamplingFrequency;
% % % % % timeVector = data.TimeVector;
% % % % % frequencyVector = data.FrequencyVector;
% % % % % 
% % % % % fpar = 50;
% % % % % % figure
% % % % % 
% % % % % % for i=1:10
% % % % % %     f = i*fpar;
% % % % % % Q=2;
% % % % % Wo = fpar/(fs/2);
% % % % % % BW = Wo/Q;
% % % % % % [b,a]=iirnotch(Wo,BW);
% % % % % % notched = filter(b,a,signal);
% % % % % % freqz(b,a,frequencyVector,fs);
% % % % % % hold all
% % % % % % end
% % % % % 
% % % % % for i=1:1
% % % % %     fpar = i*50;
% % % % %     Wo = fpar/(fs/2);
% % % % %     [b,a] = ellip(6,5,40,Wo);
% % % % %     notched = filter(b,a,signal);
% % % % %     figure
% % % % %     plot(timeVector(1:fs),notched(1:fs))
% % % % %     title(['i=',num2str(i)])
% % % % % end
% % % % % cleanSignal = struct('signal',notched,'SamplingFrequency',fs,'TimeVector',timeVector,'FrequencyVector',frequencyVector);
% % % % % 
% % % % % % figure
% % % % % % makePowerSpectrum(data);hold on
% % % % % % makePowerSpectrum(cleanSignal)
% % % % % 
% % % % % % figure
% % % % % % plot(timeVector(1:fs),signal(1:fs))
% % % % % % hold on
% % % % % plot(timeVector(1:fs),notched(1:fs))
function cleanSignal = notch(data)
close all

signal = data.signal;
fs = data.SamplingFrequency;
timeVector = data.TimeVector;
frequencyVector = data.FrequencyVector;

fpar = 150;
figure

% for i=1:10
%     f = i*fpar;
% Q=2;
Wo = fpar/(fs/2);
% BW = Wo/Q;
% [b,a]=iirnotch(Wo,BW);
% notched = filter(b,a,signal);
% freqz(b,a,frequencyVector,fs);
% hold all
% end

for i=1:10
    fpar = i*50;
    Wo = fpar/(fs/2);
    [b,a] = ellip(6,5,40,Wo);
    notched = filter(b,a,signal);
    figure
    plot(timeVector(1:fs),notched(1:fs))
    title(['i=',num2str(i)])
end
cleanSignal = struct('signal',notched,'SamplingFrequency',fs,'TimeVector',timeVector,'FrequencyVector',frequencyVector);

% figure
% makePowerSpectrum(data);hold on
% makePowerSpectrum(cleanSignal)

% figure
% plot(timeVector(1:fs),signal(1:fs))
% hold on
plot(timeVector(1:fs),notched(1:fs))
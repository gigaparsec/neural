function filtered = filterBP(data,f_highpass,f_lowpass)

signal = data.signal;
timeVector = data.TimeVector;
frequencyVector = data.FrequencyVector;
fs = data.SamplingFrequency;

% highpass
WHP = f_highpass/(0.5*fs);

[b,a] = butter(2,WHP,'high');

first = filter(b,a,signal);

%lowpass
WLP = f_lowpass/(0.5*fs);

[b,a] = butter(5,WLP,'low');

second = filtfilt(b,a,first);

filtered = struct('signal',second,'SamplingFrequency',fs,'TimeVector',timeVector,'FrequencyVector',frequencyVector);

figure
plot(timeVector,second)
title(['BandPassed (',num2str(f_highpass),',',num2str(f_lowpass),')'])
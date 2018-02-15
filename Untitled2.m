%% 
close all
fpar = 100;
% 
% 
fs = signal.SamplingFrequency;
timeVector = signal.TimeVector;
frequencyVector = signal.FrequencyVector;
sig = signal.signal;

%% create sine waves
fs = 1000;                    % Sampling frequency (samples per second)
dt = 1/fs;                   % seconds per sample
StopTime = 1;             % seconds
t = (0:dt:StopTime-dt)';     % seconds
F1 = 5;                      % Sine wave frequency (hertz)
F2 = 100;
data = sin(2*pi*F1*t)+sin(2*pi*F2*t);

plot(t,data)
%% notch for sines

% nyquist = fs/2;
% w0 = 100/nyquist;
% Q=20;
% bw = w0/Q;
% [b,a] = iirnotch(w0,bw,20);
% y = filtfilt(b,a,data);
% plot(t,data);hold on
% plot(t,y);

%% notch for data
close all

nyquist = fs/2;
filtered=sig;
for i=1:200
    fc = i*50;
    w0 = fc/nyquist;
    Q=5;
    bw = w0/Q;
    [b,a] = iirnotch(w0,bw,Q);
    filtered = filtfilt(b,a,filtered);
end

plot(timeVector,sig);
hold on
plot(timeVector,filtered)

figure
freqz(b,a,frequencyVector,fs)

figure
semilogx(frequencyVector,mag2db(abs(fft(sig))))
hold on
semilogx(frequencyVector,mag2db(abs(fft(filtered))))
xlim([0,fs/2])
%% notch for data attempt2
close all
fn = fs/2;
cf = linspace(100,1500,150);
fil = cell(length(cf),1);
for i=1:length(cf)-1
    [b,a] = butter(7,[cf(i)-1 cf(i)+1]/(0.5*fs),'stop');
    fil{i}=[b,a];
end
figure(1)
freqz(fil{1},frequencyVector,fs);
hold on
for i=2:length(cf)
    freqz(fil{i},frequencyVector,fs)
end
hold off

% filter
l = length(fil{1})/2;
for i=1:length(fil)-1
    sig = filtfilt(fil{i}(1:l),fil{i}(l:end),sig);
end

plot(sig)

%% comb
close all
q=100;
bw = 100/(fs/2);
[bC,aC] = iircomb(fs/100,bw,'notch');
figure(1)
freqz(bC,aC,frequencyVector,fs);
filtered = filtfilt(bC,aC,sig);
figure(2)
plot(timeVector,filtered)
 % NOT WORKING
%%
% 
% i=1;
% T = 20*pi;
% time = 0:0.001:T;
% fV = [1/T:1/T:length(time)/time(end)];
% 
% testS = sin(time) + sin(10*time);
% figure(i)
% plot(time,testS)
% hold all
% plot(time,sin(time))
% plot(time,sin(10*time))
% i=i+1;
% Wnotch = fpar/(0.5*(length(time)/time(end)));
% Q = 50;
% BW = Wnotch/Q;
% [b,a]=iirnotch(Wnotch,BW,Q);
% 
% figure(i)
% freqz(b,a,fV,length(time)/time(end))
% i=i+1;
% figure(i)
% semilogx(fV,mag2db(abs(fft(testS))))
% i=i+1;
% filted = filter(b,a,testS);
% figure(i)
% semilogx(fV,mag2db(abs(fft(filted))))
% i=i+1;
% figure(i)
% plot(testS)
% hold on
% plot(filted)
% i=i+1
% % 
% % % filted = filter(b,a,testS);
% % cleanSignal = struct('signal',filted,'SamplingFrequency',fs,'TimeVector',timeVector,'FrequencyVector',frequencyVector);
% % 
% % figure
% % plot(shortSignal)
% % hold on
% % plot(filted)
% % legend('orig','detrended')
% % 
% % figure
% % makePowerSpectrum(prepped);
% % hold on
% % makePowerSpectrum(cleanSignal);
% 
%% lowpass
close all

%lowpass at 300, per Kuhn et al. 2005
[bl,al] = butter(5,200/fs);
% [bH,aH] = butter(2,1/fs,'high');
[bH,aH] = ellip(5,1,100,10/fs,'high');
filtered = filtfilt(bl,al,sig);
filtered = filtfilt(bH,aH,filtered);
figure(1)
freqz(bl,al,frequencyVector,fs);
must also notch at 100 and 200

[bN,aN] = butter(3,[99 101]/(0.5*fs),'stop');



filtered = filtfilt(bN,aN,filtered);
% filtered = filtfilt(bN,aN,sig);
figure(4)
freqz(bN,aN,frequencyVector,fs)

figure(2)
% plot(timeVector,sig,'b','linewidth',2)
% hold on
plot(timeVector,filtered,'r','linewidth',2)
legend('Original','Filtered')
xlabel('Time [sec]');ylabel('Potential [\mu V]')
ylim([-2000,2000])
figure(3)
semilogx(frequencyVector,mag2db(abs((fft(detrend(filtered))))))
xlim([0 15000])
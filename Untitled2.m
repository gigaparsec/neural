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
%% notch for sines WORKS. MAKE FREQ VECTOR IN SINES TO GET A GOOD FREQZ, NOTCH MUST BE VERY NARROW

nyquist = fs/2;
w0 = 100/nyquist;
Q=5;
bw = w0/Q;
[b,a] = iirnotch(w0,bw,20);
freqz(b,a)
y = filtfilt(b,a,data);
figure
plot(t,data);hold on
plot(t,y);

%% notch for data
close all


fs1=signal1.SamplingFrequency;
fV1 = signal1.FrequencyVector;
nyquist1 = fs1/2;
% filtered=sig;

fc = 2000;
w01 = fc/nyquist1;
Q=20;
bw1 = w01/Q;
[b1,a1] = iirnotch(w01,bw1,Q);
figure(1)
freqz(b1,a1,fV1,fs1)

%%%%%

fs2 = signal2.SamplingFrequency;
fV2 = signal2.FrequencyVector;
nyquist2 = fs2/2;
% filtered=sig;
% 
% fc = 30;
w02 = fc/nyquist2;
% Q=5;
bw2 = w02/Q;
[b2,a2] = iirnotch(w02,bw2,Q);
figure(2)
freqz(b2,a2,fV2,fs2)
%%
% 
% for i=1:200
%     fc = i*50;
%     w0 = fc/nyquist;
%     Q=5;
%     bw = w0/Q;
%     [b,a] = iirnotch(w0,bw,Q);
%     filtered = filtfilt(b,a,filtered);
% end

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

%% code to plot vertical and horizontal lines in the plot
plot(signal.TimeVector,plotted)
y1 = get(gca,'ylim');
x1 = get(gca,'xlim');
hold all
plot([50 50],y1);
plot(x1,[2 2])

%% wavelet plots
close all
clc
figure
subplot(3,2,1)
[psi,xval] = wavefun('mexh',10);
plot(xval,psi,'LineWidth',2); title('Mexican Hat');
ax = gca;
ax.FontSize = 14;

subplot(3,2,2)
[phi,psi,xval] = wavefun('haar',10);
plot(xval,psi,'LineWidth',2); title('Daubechies-1 or Haar');
hold on
plot(xval,phi,'LineWidth',2);
ax = gca;
ax.FontSize = 14;
legend('\psi','\phi')

subplot(3,2,3)
[phi,psi,xval] = wavefun('db5',10);
plot(xval,psi,'LineWidth',2); title('Daubechies-5');
hold on
plot(xval,phi,'LineWidth',2);
ax = gca;
ax.FontSize = 14;
legend('\psi','\phi')

subplot(3,2,4)
[phi,psi,xval] = wavefun('db40',10);
plot(xval,psi,'LineWidth',2); title('Daubechies-40');
hold on
plot(xval,phi,'LineWidth',2);
ax = gca;
ax.FontSize = 14;
legend('\psi','\phi')

subplot(3,2,5)
[psi,xval] = wavefun('morl',10);
plot(xval,psi,'LineWidth',2); title('Morlet');
ax = gca;
ax.FontSize = 14;

subplot(3,2,6)
[psi,xval] = wavefun('gaus1',10);
plot(xval,psi,'LineWidth',2); title('Gaussian');
ax = gca;
ax.FontSize = 14;

suptitle('\fontsize{20} Mother Wavelet Examples')

%% PCA

pcaMAT = zeros(300000,32);

for i=1:32
    pcaMAT(:,i) = Channels{i}.signal;
end

[coeff,score,latent] = pca(pcaMAT);

[eigenvectors, scores] = pca(pcaMAT);

% reconstructedChannels = scores*eigenvectors' + mean of signal

%% segmentation


clc
data = OriginalSignal.signal;
tV = OriginalSignal.TimeVector;
fV = OriginalSignal.FrequencyVector;
fS = OriginalSignal.SamplingFrequency;
windowLength = ceil(fS/2);
refWindow = data(1:windowLength);

windowDuration = windowLength/fS;

windowTimeVector = 0:windowDuration/windowLength:windowDuration;

frequencyVector = 0:1/windowDuration:fs;

refPer = periodogram(refWindow);%,windowLength,frequencyVector,fS);
i=1;
for n = 2:length(data)-windowLength
    slidWindow = data(n:n+windowLength);
    slidPer = periodogram(slidWindow);%,windowLength,frequencyVector,fS);
    num = sum((refPer-slidPer).^2)/(2*pi); % spectral calculation of SEM
    den = sum((refPer.*slidPer).^2)/(4*pi^2);
    SEM(n) = num/den;
    if SEM(n)>400
        refWindow = slidWindow;
        refPer = slidPer;
        ind(i)=n;
        i=i+1;
    end
    
end
%%
data = OriginalSignal.signal;
tV = OriginalSignal.TimeVector;
fV = OriginalSignal.FrequencyVector;
fS = OriginalSignal.SamplingFrequency;
windowLength = ceil(fS/2);
refWindow = data(1:windowLength);

windowDuration = windowLength/fS;

windowTimeVector = 0:windowDuration/windowLength:windowDuration;

frequencyVector = 0:1/windowDuration:fs;

refPer = periodogram(refWindow);%,windowLength,frequencyVector,fS);
i=1;
n=2;
clear ind
%     figure(1)
%     plot(data)
%     ax=gca;
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     h3=line([1 1],ax.YLim);
%     h4=line([windowLength windowLength],ax.YLim);
%     h3.Color='r';
%     h4.Color='r';
while n<length(data)-windowLength
    slidWindow = data(n:n+windowLength);
%     figure(1)
%     delete(h1)
%     delete(h2)
%     h1=line([n n],ax.YLim);
%     h2=line([n+windowLength n+windowLength],ax.YLim);
    slidPer = periodogram(slidWindow);%,windowLength,frequencyVector,fS);
    num = sum((refPer-slidPer).^2)/(2*pi); % spectral calculation of SEM
    den = sum((refPer.*slidPer).^2)/(4*pi^2);
    SEM(n) = num/den;
    if SEM(n)>200
        refWindow = data(n+windowLength:n+2*windowLength);
        refPer = periodogram(refWindow);
        ind(i) = n+windowLength;
        i=i+1;
        n=n+windowLength+1;
        disp('segment');
%         figure(1)
%         h3=line([n+windowLength n+windowLength],ax.YLim);
%         h4=line([n+2*windowLength n+2*windowLength],ax.YLim);
%         h3.Color='r';
%         h4.Color='r';
    else
        n=n+1;
    end
    
end


plot(tV(76:end),SEM)
%% SEM segmentation, time domain approach

data = OriginalSignal.signal;
tV = OriginalSignal.TimeVector;
fV = OriginalSignal.FrequencyVector;
fS = OriginalSignal.SamplingFrequency;
windowLength = ceil(fS/2);
refWindow = data(1:windowLength);

windowDuration = windowLength/fS;

windowTimeVector = 0:windowDuration/windowLength:windowDuration;

windowFrequencyVector = 0:1/windowDuration:fs;

refPer = periodogram(refWindow);%,windowLength,frequencyVector,fS);
refWinACF = autocorr(refWindow);
n=2;
while n<length(data)-windowLength
    slidWindow = data(n:n+windowLength);
    slidWinACF = autocorr(slidWindow);
    num = sum((slidWinACF-refWinACF)^2);
%     den = 
end
%%
% test
times=0;
n=0;
while n<10
    n=n+1;
    times=times+1;
    if n==3
        n=7;
    end
end
% implement SEM
%%
% acf

N = length(slidWindow);
r_s = zeros(1,N);
temp=0;
for k = 1:N
    for n=1:N-k
    temp=temp+slidWindow(n+k)*slidWindow(n);
    end
    r_s(k)=temp/N;
end
% implement ACF, not sure this is right
r_r = zeros(1,N);
temp=0;
for k = 1:N
    for n=1:N-k
    temp=temp+refWindow(n+k)*refWindow(n);
    end
    r_r(k)=temp/N;
end

plot(r_r)
hold on
plot(r_s,'r')

[r_s2,lagsrs] = xcorr(slidWindow,slidWindow);
[r_r2,lagsrr] = xcorr(refWindow,refWindow);

figure
plot(r_r2)
hold on
plot(r_s2,'r')

%%
% data = OriginalSignal.signal;
% tV = OriginalSignal.TimeVector;
% fV = OriginalSignal.FrequencyVector;
% fS = OriginalSignal.SamplingFrequency;

data = clean.signal;
tV = clean.TimeVector;
fV = clean.FrequencyVector;
fS = clean.SamplingFrequency;
windowLength = ceil(fS/2);
refWindow = data(1:windowLength+1);

windowDuration = windowLength/fS;

windowTimeVector = 0:windowDuration/windowLength:windowDuration;

frequencyVector = 0:1/windowDuration:fs;

refACF = xcorr(refWindow,refWindow);
i=1;
n=2;
clear ind
k0 = ceil(length(refACF)/2);
while n<length(data)-windowLength
    slidWindow = data(n:n+windowLength);
    slidACF = xcorr(slidWindow,slidWindow);
    num = sum((refACF-slidACF).^2);
    den = refACF(k0)*slidACF(k0);
    SEM(n) = num/den;
    if SEM(n)>200
        refWindow = data(n+windowLength:n+2*windowLength);
        refACF = xcorr(refWindow,refWindow);
        ind(i) = n+windowLength;
        i=i+1;
        n=n+windowLength+1;
        disp('segment');
    else
        n=n+1;
    end
    
end

plot(tV(1:end-151),SEM)
ax=gca;
for i=1:length(ind)
    h=line([tV(ind(i)) tV(ind(i))],ax.YLim);
end
for j=1:length(seizureTimes)
    h=line([seizureTimes(j) seizureTimes(j)] , ax.YLim);
    h.Color='r';
end

% 
% length(refACF)
% length(slidACF)
% length(refWindow)
% length(slidWindow)

%% 
    bound = ceil(length(clean.signal)/4);
    maxLag = fs;
    sigACF = xcorr(clean.signal(1:bound),maxLag);
    sigACF = sigACF./max(sigACF);
    plot(-maxLag:1:maxLag,sigACF)
    

    train = iddata(double(clean.signal),[],1/fs);
%     temp = double(clean.signal);
    ARmod = ar(train,50,'Ts',1/fs);
    
    valid = iddata(double(clean.signal(bound:end)),[],1/fs);
    
    [YH, FIT, X0] = compare(valid,ARmod);
    compare(valid,ARmod);
    
%%

signal1 = removeDC(signal1);
fsignal1 = filterHF(signal1,1000);
fn = [100 200 300 400 500 600 700];
clean = notch(fsignal1,fn,50,'N');

fsignal2 = prepareSignal(clean.signal,clean.SamplingFrequency,2);%fs=15k
fsignal3 = prepareSignal(clean.signal,clean.SamplingFrequency,10);%fs=3k
fsignal4 = prepareSignal(clean.signal,clean.SamplingFrequency,30);%fs=1k

%signal,fs,coeff
figure(1)
plot(OriginalSignal.TimeVector,OriginalSignal.signal)
hold all
plot(fsignal1.TimeVector,fsignal1.signal,'LineWidth',2)
plot(fsignal2.TimeVector,fsignal2.signal,'LineWidth',2)
plot(fsignal3.TimeVector,fsignal3.signal,'LineWidth',2)
plot(fsignal4.TimeVector,fsignal4.signal,'LineWidth',2)
%%
figure(2)
subplot(5,1,1)
plot(OriginalSignal.TimeVector(1:30000),OriginalSignal.signal(1:30000))
ylabel({'Original Signal';'f_s=30kHz'})
set(get(gca,'YLabel'),'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14)
subplot(5,1,2)
plot(fsignal1.TimeVector(1:30000),fsignal1.signal(1:30000),'LineWidth',2)
ylabel({'Filtered Signal';'f_s=30kHz'})
set(get(gca,'YLabel'),'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14)
subplot(5,1,3)
plot(fsignal2.TimeVector(1:15000),fsignal2.signal(1:15000),'LineWidth',2)
ylabel('f_s = 15kHz')
set(get(gca,'YLabel'),'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14)
subplot(5,1,4)
plot(fsignal3.TimeVector(1:3000),fsignal3.signal(1:3000),'LineWidth',2)
ylabel('f_s = 3kHz')
set(get(gca,'YLabel'),'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14)
subplot(5,1,5)
plot(fsignal4.TimeVector(1:1000),fsignal4.signal(1:1000),'LineWidth',2)
ylabel('f_s=1kHz')
set(get(gca,'YLabel'),'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14)
xlabel('Time [sec]','FontSize',14)
suptitle('\fontsize{20} Downsampling of the Neural Recording')

%%

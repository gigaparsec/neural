%% MATCHED FILTERS. APPEND TO ARMODELS.M
FFTsNoWindow = {};
FFTsWindow = {};
Lsig = length(preSig.signal);
% directory = 'C:\Users\Jack Tchimino\Documents\testFigures_MF\noWindow';
directory = 'C:\Users\Jack Tchimino\Documents\testFigures_MF\13_Sept\12_05';
% ON = 20; % 16 is the best
progress = 0;
f = waitbar(progress,'Processing..');
pause(1e-05)
% filter function, flip template, normalize
% goodOscs = [11 12 20 33 34 35 47 50 52 57 60 64 74 75 77 94 95];
goodOscs = oscillations;
for i=1:length(goodOscs)

% ON=goodOscs(i);
ON=i;
close all
newDir = [directory,'\',num2str(ON)];
mkdir(newDir);
cd(newDir);
wind = window(@hann,length(oscillations{ON}.signal));
oscSig = oscillations{ON}.signal.*wind;

%normalize oscillation signal
oscSig = oscSig/norm(oscSig);

hs=1/fss;
filterTemplate = K*flip(oscSig); % per theory: it must be flipped
diffFilterTemplate = diff(filterTemplate)/hs;
diffFilterTemplate = diffFilterTemplate/norm(diffFilterTemplate);
figure;plot(oscSig);hold on;plot(filterTemplate);title(['Oscillation',num2str(ON)]);

r = filtfilt(filterTemplate,1,preSig.signal);
dr = filtfilt(diffFilterTemplate,1,preSig.signal);
[drUpper,~] = envelope(dr,10,'peak');
dr1=dr;
figure
h(1)=subplot(4,1,1);
plot(preSig.TimeVector,preSig.signal)
hold on
plot(preSig.TimeVector,plotted)
title('Original Waveform')
h(2)=subplot(4,1,2);
plot(preSig.TimeVector,r)
% plot(r)
% ylim([-6e04 , 6e04])
title('Filtered Waveform')
h(3)=subplot(4,1,3);
plot(preSig.TimeVector,dr)
% plot(dr);
hold on
plot(preSig.TimeVector(1:end),drUpper)
ylim([-1000 1000])
% plot(drUpper)
% ylim([-2.5e07,2.5e07]);
title('derivative of filtered')
subplot(4,1,4);
% plot(preSig.TimeVector(3:end),ddr)
plot(oscillations{ON}.TimeVector,filterTemplate)
title(['Oscillation',num2str(ON),' length=',num2str(length(oscillations{ON}.signal))])
linkaxes(h,'x')
suptitle('Matched filter output. Normalized, flipped template. filter function')
savefig(['Outputs_',num2str(ON),'.fig']);

rawStr = prepareSignal(oscillations{ON}.signal,fss,1,'n/a');
hannStr = prepareSignal(oscSig,fss,1,'n/a');


% % % % % % % % % % % % % % % % 
%        PAD THE SIGNALS      %
% % % % % % % % % % % % % % % % 

len1 = length(rawStr.signal);
len2 = length(hannStr.signal);

rem1 = 500-len1;
rem2 = 500-len2;

if mod(rem1,2)==0
    newRawStrSig = [zeros(rem1/2,1);rawStr.signal;zeros(rem1/2,1)];
else
    newRawStrSig = [zeros(ceil(rem1/2),1);rawStr.signal;zeros(floor(rem1/2),1)];
end

if mod(rem2,2)==0
    newHannStrSig = [zeros(rem2/2,1);rawStr.signal;zeros(rem2/2,1)];
else
    newHannStrSig = [zeros(ceil(rem2/2),1);rawStr.signal;zeros(floor(rem2/2),1)];
end

paddedRawStr = prepareSignal(newRawStrSig,fss,1,'n/a');
paddedHannStr = prepareSignal(newHannStrSig,fss,1,'n/a');

diffPRS = prepareSignal(diff(paddedRawStr.signal),fss,1,'n/a');
diffPHS = prepareSignal(diff(paddedHannStr.signal),fss,1,'n/a');

% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % 

figure
hF(1)=subplot(2,1,1);
semilogx(preSig.FrequencyVector,abs(fft(r)),'DisplayName','filtered')
xlim([1 500])
title('FFT of filtered')
hF(2)=subplot(2,1,2);
semilogx(preSig.FrequencyVector(1:end),abs(fft(dr)),'DisplayName','d-filtered')
xlim([1 500])
title(['FFT of d-filtered, Oscillation ',num2str(ON)])
linkaxes(h,'x')
savefig(['FFTs',num2str(ON),'.fig']);


figure
subplot(2,2,1)
h=plot(oscillations{ON}.TimeVector,oscillations{ON}.signal,'LineWidth',2);
title('Raw oscillation')
xlabel('Time [sec]');ylabel('Voltage [mV]');
axis tight
subplot(2,2,2)
plot(oscillations{ON}.TimeVector,oscSig,'LineWidth',2)
hold on
plot(oscillations{ON}.TimeVector,max(oscSig)*wind)
title('The oscillation after the application of a Hann window')
xlabel('Time [sec]');ylabel('Voltage [mV]');axis tight
subplot(2,2,3)
h=makePowerSpectrum(paddedRawStr,'-');
h.LineWidth=2;
hold on
h2 = makePowerSpectrum(diffPRS,'-');
h2.LineWidth = 2;h2.Color='r';
title('Spectrum of raw oscillation')
xlabel('Frequency [Hz]');ylabel('Magnitude [dB]');xlim([1 500]);ylim([0 80])
legend('Original','Diff')
subplot(2,2,4)
h=makePowerSpectrum(paddedHannStr,'-');
h.LineWidth=2;
hold on
h2 = makePowerSpectrum(diffPHS,'-');
h2.LineWidth = 2;h2.Color='r';
title('Spectrum of oscillation multiplied by Hann window')
xlabel('Frequency [Hz]');ylabel('Magnitude [dB]');xlim([1 500]);ylim([0 80])
suptitle(['The oscillation template, ', num2str(ON)])
legend('Original','Diff')
savefig(['Templates_',num2str(ON),'.fig']);


temp = rawStr.signal;
L1 = length(temp);
L2 = 500-length(temp);
if mod(L2,2)==0
    newS = [zeros(L2/2,1);temp;zeros(L2/2,1)];
else
    newS = [zeros(ceil(L2/2),1);temp;zeros(floor(L2/2),1)];
end
% FFTsNoWindow{ON} = (fft(newS));
FFTsNoWindow{i} = (fft(newS));

temp = hannStr.signal;
L1 = length(temp);
L2 = 500-length(temp);
if mod(L2,2)==0
    newS = [zeros(L2/2,1);temp;zeros(L2/2,1)];
else
    newS = [zeros(ceil(L2/2),1);temp;zeros(floor(L2/2),1)];
end
% FFTsWindow{ON} = (fft(newS));
FFTsWindow{i} = (fft(newS));







figure
stem(flip(hannStr.signal))
xlabel('Samples');ylabel('Absolute Magnitude');
title(['Matched Filter Impulse Response, Oscillation ', num2str(ON)])
savefig(['FIRtaps_',num2str(ON),'.fig']);


progress = progress+1/length(goodOscs);
waitbar(progress,f,'Processing..')
pause(1e-05)
end

close(f)
close all
cd('C:\Users\Jack Tchimino\Documents\GitHub\neural')
% 23200 - 
disp('finished')


%% THRESHOLD AR13

directory = 'C:\Users\Jack Tchimino\Documents\testFigures_MF\13_Sept\18_10';
% ON = 20; % 16 is the best
progress = 0;
f = waitbar(progress,'Processing..');
pause(1e-05)

goodOscs = [4 6 7 8 9 10 12 15 17 18 21 24 30 32 33 42 43 54 55 56];

threshold = 450;
minDuration = 0.02;
minDurationInSamples = 0.05*preSig.SamplingFrequency;

for i=1:length(goodOscs)
% for i=1:1
    ON = goodOscs(i);
    newDir = [directory,'\',num2str(ON)];
    mkdir(newDir);
    cd(newDir);

    
%     ON = 12;
    wind = window(@hann,length(oscillations{ON}.signal));
    oscSig = oscillations{ON}.signal.*wind;
    
    %normalize oscillation signal
    oscSig = oscSig/norm(oscSig);

    hs=1/fss;
    filterTemplate = K*flip(oscSig); % per theory: it must be flipped
    diffFilterTemplate = diff(filterTemplate)/hs;
    diffFilterTemplate = diffFilterTemplate/norm(diffFilterTemplate);
    
    r = filtfilt(filterTemplate,1,preSig.signal);
    dr = filtfilt(diffFilterTemplate,1,preSig.signal);
    [drUpper,~] = envelope(dr,10,'peak');
    dr1=dr;
    
    MFranges = find(drUpper>=threshold);
    MFranges2 = zeros(size(drUpper));
    MFranges2(MFranges) = 1;
    
    MFranges2_2 = MFranges2;
    counter = 0;
    oscCounterNew = 1;
    oscillationsNew = {};
    
    for j=2:length(MFranges2)
        if MFranges2(j)==0 && MFranges2(j-1)==0
            continue
        elseif MFranges2(j)==1 && MFranges2(j-1)==0
            counter = counter+counter+1;
        elseif MFranges2(j)==1 && MFranges2(j-1)==1
            counter = counter+1;
        elseif MFranges2(j)==0 && MFranges2(j-1)==1
            if counter<minDurationInSamples
                MFranges2_2(j-counter-1:j) = 0;
            else
                oscillationsNew{oscCounterNew}.signal = preSig.signal(j-counter-1:j);
                oscillationsNew{oscCounterNew}.TimeVector = preSig.TimeVector(j-counter-1:j);
                oscCounterNew = oscCounterNew + 1;
            end
            counter = 0;
        end
    end
    
    plottedMF = zeros(1,length(preSig.TimeVector));
    plottedMF(1:end) = NaN;
    plottedMF(find(MFranges2_2)) = preSig.signal(find(MFranges2_2));
    
    plotThreshold = threshold*ones(size(preSig.TimeVector));
    
    figure
    h(1)=subplot(4,1,1);
    plot(preSig.TimeVector,preSig.signal)
    hold on
    plot(preSig.TimeVector,plotted)
    plot(preSig.TimeVector,plottedMF,'g')
    title('Original Waveform')
    h(2)=subplot(4,1,2);
    plot(preSig.TimeVector,r)
    % plot(r)
    % ylim([-6e04 , 6e04])
    title('Filtered Waveform')
    h(3)=subplot(4,1,3);
    plot(preSig.TimeVector,dr)
    % plot(dr);
    hold on
    plot(preSig.TimeVector(1:end),drUpper)
    plot(preSig.TimeVector(1:end),plotThreshold,'g')
    ylim([-1000 1000])
    % plot(drUpper)
    % ylim([-2.5e07,2.5e07]);
    title('derivative of filtered')
    subplot(4,1,4);
    % plot(preSig.TimeVector(3:end),ddr)
    plot(oscillations{ON}.TimeVector,filterTemplate)
    title(['Oscillation',num2str(ON),' length=',num2str(length(oscillations{ON}.signal))])
    linkaxes(h,'x')
    suptitle(['Matched filter output. Threshold=',num2str(threshold),'. MinDuration=',num2str(minDuration)])
    savefig(['Outputs_',num2str(ON),'.fig']);
    
    progress = progress+1/length(goodOscs);
    waitbar(progress,f,'Processing..')
    pause(1e-05)
    
    close all
end


close(f)
close all
cd('C:\Users\Jack Tchimino\Documents\GitHub\neural')
disp('finished')


















































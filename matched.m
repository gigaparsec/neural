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
clear oscillationsNew oscillationsSampleRange
directory = 'C:\Users\Jack Tchimino\Documents\testFigures_MF\15_Sept\18_23';
% ON = 20; % 16 is the best
% progress = 0;
% f = waitbar(progress,'Processing..');
% pause(1e-05)

% goodOscs = [4 6 7 8 9 10 12 15 17 18 21 24 30 32 33 42 43 54 55 56];
goodOscs = [12];

threshold = 440;
minDuration = 0.035;
minDurationInSamples = minDuration*preSig.SamplingFrequency;

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
    
    oscillationsSampleRange = cell(1,length(oscillationsNew));
    for j=1:length(oscillationsNew)
        oscillationsSampleRange{j} = find(tV==oscillationsNew{j}.TimeVector(1)):find(tV==oscillationsNew{j}.TimeVector(end));
    end
    
    plottedMF = zeros(1,length(preSig.TimeVector));
    plottedMF(1:end) = NaN;
    plottedMF(find(MFranges2_2)) = preSig.signal(find(MFranges2_2));
    
    plotThreshold = threshold*ones(size(preSig.TimeVector));

    figure
    h(1)=subplot(4,1,1);
    plot(preSig.TimeVector,preSig.signal)
    hold on
    plot(preSig.TimeVector,manualMarker,'LineWidth',2,'Color','k')
    plot(preSig.TimeVector,plotted)
    plot(preSig.TimeVector,plottedMF,'g')
    title('Original Waveform')
    set(gca, 'FontName', 'Myriad Pro')
    h(2)=subplot(4,1,2);
    plot(preSig.TimeVector,r)
    % plot(r)
    % ylim([-6e04 , 6e04])
    title('Filtered Waveform')
    set(gca, 'FontName', 'Myriad Pro')
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
    set(gca, 'FontName', 'Myriad Pro')
    subplot(4,1,4);
    % plot(preSig.TimeVector(3:end),ddr)
    plot(oscillations{ON}.TimeVector,filterTemplate)
    title(['Oscillation',num2str(ON),' length=',num2str(length(oscillations{ON}.signal))])
    linkaxes(h,'x')
    st = suptitle(['Matched filter output. Threshold=',num2str(threshold),'. MinDuration=',num2str(minDuration)]);
    st.FontName='Myriad Pro';
%     savefig(['Outputs_',num2str(ON),'Th=',num2str(threshold),'MD=',num2str(minDuration),'.fig']);
    
%     progress = progress+1/length(goodOscs);
%     waitbar(progress,f,'Processing..')
%     pause(1e-05)
    
    


%     close all
end
%%

% close(f)
% close all
cd('C:\Users\Jack Tchimino\Documents\GitHub\neural')
disp('filtering and plot creation finished')

disp('Starting performance metric calculation')

% manual -> oscRanges
% algorithm -> oscillationsSampleRange

flagsDetectedManuals = zeros(1,length(oscRanges));
flagsDetectedAlgorithm = zeros(1,length(oscillationsSampleRange));

TP = 0; % true positives
FP = 0; % false positives (erroneous detections)
FN = 0; % false negative (oscillations that were not detected)

for i=1:length(oscillationsSampleRange)
    
    %take a detected oscillation
    detectedStart = oscillationsSampleRange{i}(1);
    detectedEnd = oscillationsSampleRange{i}(end);
    
    % look for it in the manually detected events
    for j=1:length(oscRanges)
        %take a manually detected event
        manualStart = oscRanges{j}(1);
        manualEnd = oscRanges{j}(end);
        
        if manualStart<=detectedStart && manualEnd>=detectedEnd %if the manual overlaps the detected
            TP=TP+1;
            flagsDetectedManuals(j)=1;
            flagsDetectedAlgorithm(i)=1;
            break;
        elseif manualStart>=detectedStart && manualEnd<=detectedEnd% if the detected overlaps the manual
            TP=TP+1;
            flagsDetectedManuals(j)=1;
            flagsDetectedAlgorithm(i)=1;
            break;
        elseif manualStart<detectedStart && manualEnd>detectedStart && manualEnd<detectedEnd% if the detected overlaps on the right
            TP=TP+1;
            flagsDetectedManuals(j)=1;
            flagsDetectedAlgorithm(i)=1;
            break;
        elseif manualStart>detectedStart && manualStart<detectedEnd && manualEnd>detectedEnd% if the detected overlaps on the left
            TP=TP+1;
            flagsDetectedManuals(j)=1;
            flagsDetectedAlgorithm(i)=1;
            break;
            % if no overlap is detected, continue
        end
        
    end
    
    
end

% now, flagsDetectedManuals has 1 wherever a manually designates
% oscillation was correctly identified by the algorithm, 0 wherever an
% event was not detected (False Negative). flagsDetectedAlgorithm has 1 
% wherever the detection was correct, 0 where the detection was erroneous 
% (false positive)

FN = length(find(flagsDetectedManuals==0));
FP = length(find(flagsDetectedAlgorithm==0));

% RUN IT ON SUNDAY, COUNT THE RESULTS. COMPARE TO EXCEL FILE, DEBUG

Precision = TP/(TP+FP);
Recall    = TP/(TP+FN);

disp('performance metrics calculated')

table(TP,FP,FN, Precision, Recall , threshold,minDuration)

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% PARAMETRIC ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% directory = 'C:\Users\Jack Tchimino\Documents\testFigures_MF\15_Sept\18_23';
% ON = 20; % 16 is the best
% progress = 0;
% f = waitbar(progress,'Processing..');
% pause(1e-05)
tV = BaselineRecording.TimeVector;
goodOscs = [4 6 7 8 9 10 12 15 17 18 21 24 30 32 33 42 43 54 55 56];
% goodOscs = [4];
% goodOscs = [12];
thresholdValues = 360:5:450;
minDurationValues = 0.01:0.005:0.06;
% goodOscs = [8];

% threshold = 460;
% minDuration = 0.04;
% minDurationInSamples = minDuration*preSig.SamplingFrequency;
K=1;
totalNumOfIters = length(goodOscs)*length(thresholdValues)*length(minDurationValues);

MetricsLog = cell(1,totalNumOfIters);
metricsCounter=1;
Fmeasure = zeros(length(thresholdValues),length(minDurationValues));
maxFmeasure = zeros(1,length(goodOscs));
tempMax = 0;
maxFmeasures = zeros(1,length(goodOscs));
FmeasureArray = cell(1,length(goodOscs));
for i=1:length(goodOscs)
% for i=1:1
%     i
    ON = goodOscs(i);
    
%     newDir = [directory,'\',num2str(ON)];
%     mkdir(newDir);
%     cd(newDir);

    
%     ON = 12;
    wind = window(@hann,length(oscillations{ON}.signal));
    oscSig = oscillations{ON}.signal.*wind;
    
    %normalize oscillation signal
    oscSig = oscSig/norm(oscSig);

    hs=1/fss;
    filterTemplate = K*flip(oscSig); % per theory: it must be flipped
    diffFilterTemplate = diff(filterTemplate)/hs;
    diffFilterTemplate = diffFilterTemplate/norm(diffFilterTemplate);
    
    r = filtfilt(filterTemplate,1,BaselineRecording.signal);
    dr = filtfilt(diffFilterTemplate,1,BaselineRecording.signal);
    [drUpper,~] = envelope(dr,10,'peak');
    dr1=dr;
   
    for j=1:length(thresholdValues)
        threshold = thresholdValues(j);
%         j
        for k=1:length(minDurationValues)
        minDuration = minDurationValues(k);
        minDurationInSamples = minDuration*BaselineRecording.SamplingFrequency;
%         k
    MFranges = find(drUpper>=threshold);
    MFranges2 = zeros(size(drUpper));
    MFranges2(MFranges) = 1;
    
    MFranges2_2 = MFranges2;
    counter = 0;
    oscCounterNew = 1;
    oscillationsNew = {};
    
    for jjj=2:length(MFranges2)
        if MFranges2(jjj)==0 && MFranges2(jjj-1)==0
            continue
        elseif MFranges2(jjj)==1 && MFranges2(jjj-1)==0
            counter = counter+counter+1;
        elseif MFranges2(jjj)==1 && MFranges2(jjj-1)==1
            counter = counter+1;
        elseif MFranges2(jjj)==0 && MFranges2(jjj-1)==1
            if counter<minDurationInSamples
                MFranges2_2(jjj-counter-1:jjj) = 0;
            else
                oscillationsNew{oscCounterNew}.signal = BaselineRecording.signal(jjj-counter-1:jjj);
                oscillationsNew{oscCounterNew}.TimeVector = BaselineRecording.TimeVector(jjj-counter-1:jjj);
                oscCounterNew = oscCounterNew + 1;
            end
            counter = 0;
        end
    end
    
    oscillationsSampleRange = cell(1,length(oscillationsNew));
    for jjj=1:length(oscillationsNew)
        oscillationsSampleRange{jjj} = find(tV==oscillationsNew{jjj}.TimeVector(1)):find(tV==oscillationsNew{jjj}.TimeVector(end));
    end
    
    plottedMF = zeros(1,length(BaselineRecording.TimeVector));
    plottedMF(1:end) = NaN;
    plottedMF(find(MFranges2_2)) = BaselineRecording.signal(find(MFranges2_2));
    
    plotThreshold = threshold*ones(size(BaselineRecording.TimeVector));

    flagsDetectedManuals = zeros(1,length(oscRanges));
    flagsDetectedAlgorithm = zeros(1,length(oscillationsSampleRange));

    TP = 0; % true positives
    FP = 0; % false positives (erroneous detections)
    FN = 0; % false negative (oscillations that were not detected)

    for ii=1:length(oscillationsSampleRange)

        %take a detected oscillation
        detectedStart = oscillationsSampleRange{ii}(1);
        detectedEnd = oscillationsSampleRange{ii}(end);

        % look for it in the manually detected events
        for jj=1:length(oscRanges)
            %take a manually detected event
            manualStart = oscRanges{jj}(1);
            manualEnd = oscRanges{jj}(end);

            if manualStart<=detectedStart && manualEnd>=detectedEnd %if the manual overlaps the detected
                TP=TP+1;
                flagsDetectedManuals(jj)=1;
                flagsDetectedAlgorithm(ii)=1;
                break;
            elseif manualStart>=detectedStart && manualEnd<=detectedEnd% if the detected overlaps the manual
                TP=TP+1;
                flagsDetectedManuals(jj)=1;
                flagsDetectedAlgorithm(ii)=1;
                break;
            elseif manualStart<detectedStart && manualEnd>detectedStart && manualEnd<detectedEnd% if the detected overlaps on the right
                TP=TP+1;
                flagsDetectedManuals(jj)=1;
                flagsDetectedAlgorithm(ii)=1;
                break;
            elseif manualStart>detectedStart && manualStart<detectedEnd && manualEnd>detectedEnd% if the detected overlaps on the left
                TP=TP+1;
                flagsDetectedManuals(jj)=1;
                flagsDetectedAlgorithm(ii)=1;
                break;
                % if no overlap is detected, continue
            end

        end


    end

    FN = length(find(flagsDetectedManuals==0));
    FP = length(find(flagsDetectedAlgorithm==0));

    % RUN IT ON SUNDAY, COUNT THE RESULTS. COMPARE TO EXCEL FILE, DEBUG

    Precision = TP/(TP+FP);
    Recall    = TP/(TP+FN);
    Fmeasure(j,k) = 2*(Precision*Recall)/(Precision+Recall);
    MetricsLog{metricsCounter}.Oscillation = ON;
    MetricsLog{metricsCounter}.TruePositives = TP;
    MetricsLog{metricsCounter}.FalsePositives = FP;
    MetricsLog{metricsCounter}.FalseNegatives = FN;
    MetricsLog{metricsCounter}.Precision = Precision;
    MetricsLog{metricsCounter}.Recall = Recall;
    MetricsLog{metricsCounter}.Threshold = threshold;
    MetricsLog{metricsCounter}.MinimumDuration = minDuration;
    MetricsLog{metricsCounter}.Fmeasure = Fmeasure(j,k);
    metricsCounter = metricsCounter+1;
        end
    end
    FmeasureArray{i}=Fmeasure;
    maxFmeasures(i) = max(max(Fmeasure));
end
disp('finished')

figure
stem(goodOscs,maxFmeasures)
clear BestFscores ThVal MDVal
for O=1:length(goodOscs)

[th,md]=find(FmeasureArray{O}==max(max(FmeasureArray{O})));
    
BestFscores(O)=max(max(FmeasureArray{O}));
ThVal(O) = thresholdValues(th);
MDVal(O) = minDurationValues(md);

% figure
% surf(minDurationValues,thresholdValues,FmeasureArray{O})
% colormap(winter)
% xlabel('Minimum Duration [sec]');ylabel('Amplitude');zlabel('F-measure');
% title(['F-measure values for a template number',num2str(goodOscs(O)),' max(Fmeas)=',num2str(max(max(FmeasureArray{O})))]);
% set(gca, 'FontName', 'Myriad Pro')
% pause
% close all
end
final =table(goodOscs',BestFscores',ThVal',MDVal');
final.Properties.VariableNames = {'Oscillation_Number','F_measure','Threshold','Min_Duration'}
disp('finished')


%%
ON=32;
wind = window(@hann,length(oscillations{ON}.signal));
    oscSig = oscillations{ON}.signal.*wind;
    
    %normalize oscillation signal
    oscSig = oscSig/norm(oscSig);

    hs=1/fss;
    filterTemplate = K*flip(oscSig); % per theory: it must be flipped
    diffFilterTemplate = diff(filterTemplate)/hs;
    diffFilterTemplate = diffFilterTemplate/norm(diffFilterTemplate);


figure
plot(oscillations{12}.TimeVector,filterTemplate,'LineWidth',2)
hold on
plot(oscillations{12}.TimeVector(2:end),diffFilterTemplate,'LineWidth',2)
xlabel('Time [sec]');ylabel('Normalized Amplitude');
legend('Original','Differentiated')
title('The final selection for the matched filter template');
set(gca, 'FontName', 'Myriad Pro')


% temp = cellfun(@(x) x.Precision,MetricsLog,'UniformOutput',false);
% precisions = cell2mat(temp);
% 
% temp = cellfun(@(x) x.Recall,MetricsLog,'UniformOutput',false);
% recalls = cell2mat(temp);
% 
% clc
% goodPrec = find(precisions>0.85);
% goodRec  = find(recalls>0.7);
% % table(length(goodPrec),length(goodRec))
% 
% goodRecInd = [];
% counter=1;
% for i = 1:length(goodRec)
%     
%     ind = find(goodPrec==goodRec(i));
%     if ~isempty(ind)
%         goodRecInd(counter)=i;
%         counter=counter+1;
%     end
%     
% end
% % goodRecInd
% 
% MetricsLog{goodRec(goodRecInd)}





% close(f)
% % close all
% cd('C:\Users\Jack Tchimino\Documents\GitHub\neural')
% disp('filtering and plot creation finished')
% 
% disp('Starting performance metric calculation')
% 
% % manual -> oscRanges
% % algorithm -> oscillationsSampleRange
% 
% 
% % now, flagsDetectedManuals has 1 wherever a manually designates
% % oscillation was correctly identified by the algorithm, 0 wherever an
% % event was not detected (False Negative). flagsDetectedAlgorithm has 1 
% % wherever the detection was correct, 0 where the detection was erroneous 
% % (false positive)
% 
% 
% disp('performance metrics calculated')
% 
% 
% % table(TP,FP,FN, Precision, Recall , threshold,minDuration)
% 
% 
% 








































%%baseline
%% manual characterization of oscillations

TimeInstantsJack = BaselineOscManual2;

% TimeInstantsFarnaz = FarnazMarkingBaseline;
            % maybe 53.8
            
SampleNum = TimeInstantsJack;
SampleNumStart = SampleNum-50;
SampleNumEnd = SampleNum+50;

% SampleNumFarnaz = TimeInstantsFarnaz;
% SampleNumStartFarnaz = SampleNumFarnaz;
% SampleNumEndFarnaz = SampleNumFarnaz+100;

numOfOscs = length(TimeInstantsJack);
% numOfOscsFarnaz = length(TimeInstantsFarnaz);

oscRanges = cell(1,numOfOscs);
% oscRangesF= cell(1,numOfOscsFarnaz);

for i=1:numOfOscs
    oscRanges{i} = SampleNumStart(i):SampleNumEnd(i);
end

% for i=1:numOfOscsFarnaz
%     oscRangesFarnaz{i} = SampleNumStartFarnaz(i):SampleNumEndFarnaz(i);
% end

manualMarker = zeros(1,length(BaselineRecording.TimeVector));

% manualMarkerF = zeros(1,length(BaselineRecording.TimeVector));

manualMarker(1:end) = NaN;
% manualMarkerF(1:end) = NaN;

for i=1:numOfOscs
    manualMarker(oscRanges{i})=BaselineRecording.signal(oscRanges{i});
end

% for i=1:numOfOscsFarnaz
%     manualMarkerF(oscRangesFarnaz{i})=BaselineRecording.signal(oscRangesFarnaz{i});
% end

figure
% h(1)=subplot(2,1,1);
plot(BaselineRecording.signal)
hold on
plot(manualMarker,'LineWidth',2,'Color','k');ylim([-150 150])
ax = gca;
ax.XAxis.Exponent = 0;
grid on
title('Jack')
% h(2)=subplot(2,1,2);
% plot(BaselineRecording.signal)
% hold on
% plot(manualMarkerF,'LineWidth',2,'Color','k');ylim([-150 150])
% ax = gca;
% ax.XAxis.Exponent = 0;
% grid on
% title('Farnaz')
% linkaxes(h,'x')

%% stimulation

%% manual characterization of oscillations

% TimeInstantsJack = StimOscManual;
TimeInstantsJack = StimNew2;

% TimeInstantsFarnaz = FarnazMarkingStim;
            % maybe 53.8
            
SampleNum = TimeInstantsJack;
SampleNumStart = SampleNum-50;
SampleNumEnd = SampleNum+50;

% SampleNumFarnaz = TimeInstantsFarnaz;
% SampleNumStartFarnaz = SampleNumFarnaz;
% SampleNumEndFarnaz = SampleNumFarnaz+100;

numOfOscs = length(TimeInstantsJack);
% numOfOscsFarnaz = length(TimeInstantsFarnaz);

oscRangesStim = cell(1,numOfOscs);
% oscRangesF= cell(1,numOfOscsFarnaz);

for i=1:numOfOscs
    oscRangesStim{i} = SampleNumStart(i):SampleNumEnd(i);
end

for i=1:numOfOscsFarnaz
    oscRangesF{i} = SampleNumStartFarnaz(i):SampleNumEndFarnaz(i);
end

manualMarkerStimNew = zeros(1,length(StimulationRecording.TimeVector));

manualMarkerF = zeros(1,length(StimulationRecording.TimeVector));

manualMarkerStimNew(1:end) = NaN;
manualMarkerF(1:end) = NaN;

for i=1:numOfOscs
    manualMarkerStimNew(oscRangesStim{i})=StimulationFilteredAt150.signal(oscRangesStim{i});
end

for i=1:numOfOscsFarnaz
    manualMarkerF(oscRangesF{i})=StimulationRecording.signal(oscRangesF{i});
end

figure
h(1)=subplot(2,1,1);
plot(StimulationRecording.signal)
hold on
plot(manualMarkerStimNew,'LineWidth',2,'Color','k');ylim([-150 150])
ax = gca;
ax.XAxis.Exponent = 0;
grid on
title('Jack')
h(2)=subplot(2,1,2);
plot(StimulationRecording.signal)
hold on
plot(manualMarkerF,'LineWidth',2,'Color','k');ylim([-150 150])
ax = gca;
ax.XAxis.Exponent = 0;
grid on
title('Farnaz')
linkaxes(h,'x')


























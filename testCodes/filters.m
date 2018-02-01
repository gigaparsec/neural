% Author: Jack Tchimino

%% peak detection and averaging

inverted = tenMinSignal*(-1);

[peak,loc] = findpeaks(inverted,'MinPeakProminence',5e-4); % loc contains the indices of the most prominent peaks
figure
plot(inverted);findpeaks(inverted,'MinPeakProminence',5e-4);

% simple averaging (Ragnayyan, Biomedical signal Analysis p.143-4)
% the peak is regarded as a reference point
summed = zeros(51); % each peak has been observed to be ~50 samples in length
figure
for i=2:1:length(loc-1) %dont grab the first and second, in case they're too clos to the edges
    summed = summed+tenMinSignal(loc(i)-25:loc(i)+25);
    plot(tenMinSignal(loc(i)-25:loc(i)+25)) % plot all the peaks in one figure
    hold all
end

averaged = summed/(length(loc)-2);
figure
plot(averaged)

summed = zeros(1,51);

for i=1:size(spikes,1)
    summed = summed+spikes(i,:);
    plot(spikes(i,:))
    hold all
end

figure
plot(summed/1153)

%%
% moving average filter
a = 1;
filterLength = 5;
b = [ones(filterLength,1)]/filterLength;
% the default filter function implements a moving average filter
clean = filter(b,a,tenMinSignal);

plot(tenMinSignal);
hold on
plot(clean,'r');

%% Hann filter
b = [1 2 1]/4;
a=1;

clearHann = filter(b,a,tenMinSignal);

plot(tenMinSignal);
hold on
plot(clearHann,'r');


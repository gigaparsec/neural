%% averaging display

clear all
close all
x = [0:0.001:2*pi];

y = sin(x);

numOfTests = 1000;

for i=1:numOfTests
noise(i,:) = awgn(y,5);
end
% 
% figure(1)
% for i=1:1:numOfTests
% plot(noise(i,:));
% hold all
% end
% hold all
% plot(y)

avged = zeros(size(noise(1,:)));

for i=1:numOfTests
    avged = avged+noise(i,:);
end
avged = avged/numOfTests;
figure(2)
plot(avged)
hold on
plot(y)
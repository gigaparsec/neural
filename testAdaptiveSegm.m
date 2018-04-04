% x = linspace(0,8);
% y = sin(x);
% figure
% h = plot(x,y);
% pause
% % h.XDataSource = 'x';
% h.YDataSource = 'y';
% 
% y = sin(x.^3);
% refreshdata
%% test adaptive seg
% [SB,SEM30] = adaptiveSegm(OriginalSignal,50000,30);
figure(2);plot(OriginalSignal.TimeVector,SEM30);
% [SB,SEM60] = adaptiveSegm(OriginalSignal,50000,60);
figure(3);plot(OriginalSignal.TimeVector,SEM60)
% [SB,SEM150] = adaptiveSegm(OriginalSignal,50000,150);
figure(4);plot(OriginalSignal.TimeVector,SEM150)
% [SB,SEM300] = adaptiveSegm(OriginalSignal,50000,300);
figure(5);plot(OriginalSignal.TimeVector,SEM300)
% [SB,SEM600] = adaptiveSegm(OriginalSignal,50000,600);
figure(6);plot(OriginalSignal.TimeVector,SEM600)
% [SB,SEM900] = adaptiveSegm(OriginalSignal,50000,900);
figure(7);plot(OriginalSignal.TimeVector,SEM900)
% [SB,SEM1200] = adaptiveSegm(OriginalSignal,50000,1200);
figure(8);plot(OriginalSignal.TimeVector,SEM1200)

legend('30','60','150','300','600','900','1200')

%%
close(hfig)
[SB,SEM180] = adaptiveSegm(OriginalSignal,180,150);
tV = OriginalSignal.TimeVector;
hfig = figure(10);plot(tV,SEM180);hold all
ax = gca;
for i=1:length(SB)
    figure(10)
    h1=line([tV(SB(i)) tV(SB(i))],[ax.YLim(1) ax.YLim(2)]);
    if mod(i,2)==0
        h1.Color='g';
    else
        h1.Color='r';
    end
    h1.LineWidth=2;
end
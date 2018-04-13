sig = clean.signal;
SL = length(sig);

level = 6;

rem = mod(2^level,SL);

fl = floor(SL/(2^level));

[SWA, SWD] = swt(sig(1:64*466),5,'coif5');

[cA,cD] = dwt(sig,'db44');
% %%
figure(1)
subplot(7,1,1)
plot(sig)
subplot(7,1,2)
plot(SWA(1,:))
subplot(7,1,3)
plot(SWA(2,:))
subplot(7,1,4)
plot(SWA(3,:))
subplot(7,1,5)
plot(SWA(4,:))
subplot(7,1,6)
plot(SWA(5,:))
subplot(7,1,7)
% plot(SWA(6,:))


figure(2)
subplot(7,1,1)
plot(sig)
subplot(7,1,2)
plot(SWD(1,:))
subplot(7,1,3)
plot(SWD(2,:))
subplot(7,1,4)
plot(SWD(3,:))
subplot(7,1,5)
plot(SWD(4,:))
subplot(7,1,6)
plot(SWD(5,:))
subplot(7,1,7)
% plot(SWD(6,:))


% subplot(8,1,1)
% plot(sig)
% 
% subplot(8,1,2)
% plot(swc(1,:))
% ylabel('D1')
% set(gca,'ytick',[])
% 
% subplot(8,1,3)
% plot(swc(2,:))
% ylabel('D3')
% set(gca,'ytick',[])
% 
% subplot(8,1,4)
% plot(swc(3,:))
% ylabel('D3')
% set(gca,'ytick',[])
% 
% subplot(8,1,5)
% plot(swc(4,:))
% ylabel('A4')
% set(gca,'ytick',[])
% 
% subplot(8,1,6)
% plot(swc(5,:))
% ylabel('D3')
% set(gca,'ytick',[])
% 
% subplot(8,1,7)
% plot(swc(6,:))
% ylabel('D3')
% set(gca,'ytick',[])
% 
% subplot(8,1,8)
% plot(swc(7,:))
% ylabel('D3')
% set(gca,'ytick',[])

%%

[c,l] = wavedec(sig,3,'db44');
approx = appcoef(c,l,'db44');
[cd1,cd2,cd3] = detcoef(c,l,[1 2 3]);

subplot(4,1,1)
plot(approx)
title('Approximation Coefficients')
subplot(4,1,2)
plot(cd3)
title('Level 3 Detail Coefficients')
subplot(4,1,3)
plot(cd2)
title('Level 2 Detail Coefficients')
subplot(4,1,4)
plot(cd1)
title('Level 1 Detail Coefficients')

%%

[c,l] = wavedec(sig,6,'db44');
[c,l] = dwt(sig,'db44');
approx = appcoef(c,l,'db44');
[cd1,cd2,cd3,cd4,cd5,cd6] = detcoef(c,l,[1 2 3 4 5 6]);

subplot(7,1,1)
plot(approx)
title('Approximation Coefficients')
subplot(7,1,2)
plot(cd6)
title('Level 6 Detail Coefficients')
subplot(7,1,3)
plot(cd5)
title('Level 5 Detail Coefficients')
subplot(7,1,4)
plot(cd4)
title('Level 4 Detail Coefficients')
subplot(7,1,5)
plot(cd3)
title('Level 3 Detail Coefficients')
subplot(7,1,6)
plot(cd2)
title('Level 2 Detail Coefficients')
subplot(7,1,7)
plot(cd1)
title('Level 1 Detail Coefficients')


%% decomposition

% create filters for daubechies 44
[LO_D,HI_D,LO_R,HI_R] = wfilters('db44');

N = fix(log2(length(sig)));

% deconstruct 
[c,l] = wavedec(sig,6,LO_D,HI_D);

approx = appcoef(c,l,'db44');
[cd1,cd2,cd3,cd4,cd5,cd6] = detcoef(c,l,[1:1:6]);

% plot

subplot(7,1,1);plot(approx);title('Approximation Coefficients/Level 6 Lowpass coefficients');
subplot(7,1,2);plot(cd6);title('Level 6 Detail Coefficients')
subplot(7,1,3)
plot(cd5)
title('Level 5 Detail Coefficients')
subplot(7,1,4)
plot(cd4)
title('Level 4 Detail Coefficients')
subplot(7,1,5)
plot(cd3)
title('Level 3 Detail Coefficients')
subplot(7,1,6)
plot(cd2)
title('Level 2 Detail Coefficients')
subplot(7,1,7)
plot(cd1)
title('Level 1 Detail Coefficients')

% reconstruct!

X = waverec(c,l,LO_R,HI_R);

% reconstruct with specific bands
% replace the "cd" arrays with zeros that will not be used in the
% reconstruction with zeros
newC = [approx;cd6;zeros(size(cd5));zeros(size(cd4));zeros(size(cd3));zeros(size(cd2));zeros(size(cd1))];

newX = waverec(newC,l,LO_R,HI_R);
plot(X)
hold on
plot(newX)

%% MULTIRESOLUTION ANALYSIS

wt = modwt(clean.signal,'haar',4);
mra = modwtmra(wt,'haar');

figure(1)
subplot(5,1,1)
plot(wt(1,:))
hold all
for i=2:5
    subplot(5,1,i)
    plot(wt(i,:))
end

figure(2)
subplot(5,1,1)
plot(mra(1,:))
hold all
for i=2:5
    subplot(5,1,i)
    plot(mra(i,:))
end

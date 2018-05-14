%% Create dummy data. noise and some sinusoids in different levels and at different frequencies
data = zeros(32,30000);

noisy = awgn(data,20);

% plot(noisy(1,:));

N=1000;

fs = 500;

w = (1:N)*2*pi/fs;

x5 = .75 *sin(w*100);
x7 = .4  *sin(w*7);
x2 = .1  *sin(w*2);
x9 = .9  *sin(w*9);


D = zeros(32,1000);
for i=1:32
D(i,:) = rand*x5 + rand*x7 + rand*x2 + rand*x9;
end

D = awgn(0.4*D,20);

full = [noisy,D,noisy];
% plot(full(1,:))
% spectrogram(full(1,:),ceil(1*fs),ceil(0.9*fs),ceil(1*fs),fs,'yaxis')
% view(-45,65)

%% Apply svd

[U,S,V] = svdecon(full);

%% get one component

Stemp = zeros(size(S));

lambda = 1;
Stemp(lambda,lambda) = S(lambda,lambda);

full_partial_reconstr = U*Stemp*V';

plot(full_partial_reconstr(2,:))

%% PCA

[coeff,score,latent] = pca(full');
[coe2,score2] = pca(full');
% mu = mean(full');
mu2 = mean(full,2);
temp  = score(:,1:4)*coeff(:,1:4)' + mu;
el = 2;
figure
plot(full_partial_reconstr(1,:))
hold on
% figure
plot(score(:,1))

plot(temp(:,1))
hold on
plot(full(1,:))

isequal(temp(:,1),full(1,:)')


%%

[testV,~] = eig(full*full');
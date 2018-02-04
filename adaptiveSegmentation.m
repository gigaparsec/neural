%% Author: Jack Tchimino

% Adaptive segmentation function. This function will divide the signal into
% segments, using the spectral error measure criterion and the ACF criterion. 
% follows the algorithm in Ragnayyan's Biomedical Signal Analysis, page 489

%%
% threshold: the SEM threshold at which the signal will be segmented. P is
% the order of the AR model that willbe used. plotYN is a "Y" or "N" input,
% which makes the function generate the plots or not. method can be ACF or
% SEM. N: number of samples in window

function segments = adaptiveSegmentation(signal,fs,threshold,N,P,plotYN,method)

    n0 = N+P;

if strcmp(method,'SEM')==1
    initWindow = signal(n0-N:n0+N);
    ACF = autocorr(initWindow,P);
    temp = iddata(double(initWindow),[],1/fs);
    ARmod = ar(temp,P);
    predError = pe(ARmod,(initWindow),N);
    return;
    
    
    
elseif strcmp(method,'ACF')==1
        
else
     disp("Error: Try again");
    return;
end   
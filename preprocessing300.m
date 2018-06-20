function output = preprocessing300(sig,fs,fn,directory)

% filter at 300, notch, subsample at 1kHz

% temp = removeDC(prepareSignal(sig,fs,1,directory),5);
temp=prepareSignal(sig,fs,1,directory);
if isempty(fn)
    tempOut = filterHF(temp,300);
else
    tempF = filterHF(temp,300);
    tempOut = customNotch(tempF,fn,250,'N');
    tempOut = customNotch(tempOut,[30,50,90,120,180,210,240,270,330,360,390,420],100,'N');
end

output = prepareSignal(tempOut.signal,fs,30,directory);
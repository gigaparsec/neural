function preprocessedSignal = preprocessSignal(data,HighPassCutoff,LowPassCutoff,fNotch,subsamplingFactor)

signal = data.signal;
% tV = data.TimeVector;
% fV = data.FrequencyVector;
fs = data.SamplingFrequency;
Directory = data.FilePath;
Q = 70;
if HighPassCutoff>=LowPassCutoff
    error('High Pass filter cutoff cannot be higher than Low Pass Cutoff');
end
if subsamplingFactor<1 || subsamplingFactor>30000
    error('Invalid subsampling factor');
end

WLP = LowPassCutoff/(0.5*fs);
[b,a] = butter(5,WLP,'low');
preprocessed = filtfilt(b,a,signal);

if HighPassCutoff>0
    WHP = HighPassCutoff/(0.5*fs);
    [b,a] = butter(4,WHP,'high');
    preprocessed = filtfilt(b,a,preprocessed);
end

if ~isempty(fNotch)
    nyquist = fs/2;
    w0 = fNotch/nyquist;
    bw = w0/Q;
    for i=1:length(fNotch)
        [b,a] = iirnotch(w0(i),bw(i));%,20);
        preprocessed = filtfilt(b,a,preprocessed);
    end
end



preprocessedSignal = prepareSignal(preprocessed,fs,subsamplingFactor,Directory);


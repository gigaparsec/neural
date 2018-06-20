folder = uigetdir('C:\Users\Jack Tchimino\Documents\Recordings');

for i=1:1:32
    recPath = [folder,'\100_CH',num2str(i),'.continuous'];
    [raw,~,info] = load_open_ephys_data_faster(recPath);
    fs = info.header.sampleRate;
    clear info
    numOfSamples = length(raw);%ceil(seconds*fs);
    shortSignal = raw(1:numOfSamples);
    Channel = prepareSignal(shortSignal,fs,1);
    Channel.name = ['Channel ',num2str(i)];
    cd 'C:\Users\Jack Tchimino\Documents\MATLAB\Thesis\signal Mat files'
    save(['Channel',num2str(i)],'Channel')
    clear raw numOfSamples shortSignal Channel
    disp(['channel ',num2str(i),' done'])
end
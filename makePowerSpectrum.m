function h = makePowerSpectrum(data)

signal = data.signal;
fV = data.FrequencyVector;
fs = data.SamplingFrequency;

h=semilogx(fV,mag2db(abs((fft(detrend(signal))))),'x');
xlim([10^0,fs/2]);
grid on
xlabel('Frequency [Hz]');ylabel('Magnitude [db]');
% title('Neural Signal Power Spectrum');

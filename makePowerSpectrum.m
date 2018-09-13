function h = makePowerSpectrum(data,ind)

signal = data.signal;
fV = data.FrequencyVector;
fs = data.SamplingFrequency;

h=semilogx(fV,mag2db(abs((fft(detrend(signal))))),ind);
% h=semilogx(fV,(abs((fft(detrend(signal))))),ind);
xlim([10^-1,fs/2]);
grid on
xlabel('Frequency [Hz]');ylabel('Magnitude [db]');
% title('Neural Signal Power Spectrum');

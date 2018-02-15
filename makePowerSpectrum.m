function makePowerSpectrum(data)

signal = data.signal;
fV = data.FrequencyVector;

semilogx(fV,mag2db(abs((fft(detrend(signal))))),'x')
% xlim([10^0,1e4]);
grid on
xlabel('Frequency [Hz]');ylabel('Magnitude [db]');
title('Neural Signal Power Spectrum');

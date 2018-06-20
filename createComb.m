function [b,a] = createComb(freq,Q,fs)

b = cell(1,length(freq));
a = cell(1,length(freq));

nyquist = fs/2;
w0 = freq/nyquist;
bw = w0/Q;

for i=1:length(freq)
    [b_temp,a_temp] = iirnotch(w0(i),bw(i));
    b{i} = b_temp;
    a{i} = a_temp;
end
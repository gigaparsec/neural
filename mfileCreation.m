%%
fileID = fopen('raw2_1kHz.dat','w');
fwrite(fileID, raw2.signal,'double');
fclose(fileID);

m = memmapfile('raw.dat','Format','double');
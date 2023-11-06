function [f,f_pos] = freq_fft(L,spacelength)
% Function to get frequencies for fft
% f: Map to negative and positive frequencies: [0,fs,..,N/2*fs,-(N/2)*fs,...,fs]
% fp: Map to positive frequencies: [0,fs,...,N/2*fs]
fs = L/spacelength;
f = (0:L-1)*fs/L;
if mod(L,2)==0
    f(L/2+1:end) = f(L/2+1:end)-fs;
else
    f((L+1)/2+1:end) = f((L+1)/2+1:end)-fs;
end
f_pos = abs(f(1:floor(L/2)+1));

end
function [Y,ft] = fft_time_z(signal,time,nd,q,tap)

if q<1
    %ndq = (nd -1)/q
    %ndq = (nd)/(1-q)-1;
    ndq = nd + (nd-1)*((1/(1-q))-1);
else
    ndq=nd;
end

[NT,NZ] = size(signal);
%exp1 = nextpow2(NT)-1;
%NT1 = 2^exp1;
%signal = signal_in(1:NT1,:);
N = NT/nd;

qm = 1-q;
%overlapping
for i=1:ndq
    ind = int64(qm*(i-1)*N+1:(qm*(i-1)+1)*N);
    Pn1{i} = signal(ind,:);
    size(Pn1{i});
end


ndq = fix(ndq);
%N = 2^nextpow2(N);
N = 2*N;
Y = zeros(ndq,N,NZ);

if tap
    tapF = repmat(hanning(N),[1 NZ]);
    tapCoef = sqrt(8/3);
else
    tapF = repmat(ones(N,1),[1 NZ]);
    tapCoef = 1;
end

%FFT in time
for i=1:ndq
    %func = Pn1{i};
    func = [Pn1{i};zeros(N,NZ)];
    Y(i,:,:) = tapCoef*fft(func.*tapF,[],1);
end

[ft,ftp] = freq_fft(N,time(N)-time(1));

end

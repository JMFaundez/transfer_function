function [Sxx,Syy,Sxy,ft,fz,dat]= ordinary_spectra(xdata,ydata,time,z0,nd,q,tap)

%NS = 31; %number of sensors in one row
NS = length(z0);

dt = time(2)-time(1);

nt = length(time);
xn = xdata;%-avgx;
yn = ydata;%-avgy;

NT = length(time);

Xz = fft(xn,[],2);
Yz = fft(yn,[],2);

Lz = z0(end) - z0(1) + (z0(2)-z0(1));
[fz,fzp] = freq_fft(NS,Lz);

[X, ft] = fft_time_z(Xz,time,nd,q,tap);
[Y, ft] = fft_time_z(Yz,time,nd,q,tap);

[ndq,N,nz] = size(X);
Sxx = zeros(ndq,N,nz);
Syy = zeros(ndq,N,nz);
Sxy = zeros(ndq,N,nz);
for i=1:ndq
	Sxx(i,:,:) = squeeze(conj(X(i,:,:))).*squeeze(X(i,:,:));
	Syy(i,:,:) = squeeze(conj(Y(i,:,:))).*squeeze(Y(i,:,:));
	Sxy(i,:,:) = squeeze(conj(X(i,:,:))).*squeeze(Y(i,:,:));
end

Sxx = squeeze(mean(Sxx,1));
Syy = squeeze(mean(Syy,1));
Sxy = squeeze(mean(Sxy,1));


dat.time = time;
dat.xn = xn;
dat.yn = yn;
dat.X = squeeze(mean(X,1));
dat.Y = squeeze(mean(Y,1));

end

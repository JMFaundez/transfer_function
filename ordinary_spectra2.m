function [Sxx,Syy,Sxy,ft,fz]= ordinary_spectra2(xdata,ydata,time,z0,nd,q,tap)

%NS = 31; %number of sensors in one row
nz = length(z0);

dt = time(2)-time(1);

nt = length(time);
xn = xdata;%-avgx;
yn = ydata;%-avgy;


% fft along span
Xz = fft(xn,[],2);
Yz = fft(yn,[],2);

Lz = z0(end) - z0(1) + (z0(2)-z0(1));
[fz,fzp] = freq_fft(nz,Lz);

N = 2^nextpow2(nt/nd);

test = cpsd(Xz(:,1),Xz(:,1),fix(nt/nd),fix(q*nt/nd));
[N,~] = size(test);

Sxx = cpsd(Xz,Xz,fix(nt/nd),fix(q*nt/nd));
Syy = cpsd(Yz,Yz,fix(nt/nd),fix(q*nt/nd));
Sxy = cpsd(Xz,Yz,fix(nt/nd),fix(q*nt/nd));
[N,~] = size(Sxx)
%Sxx = zeros(N,nz);
%Syy = zeros(N,nz);
%Sxy = zeros(N,nz);

%for i=1:nz 
%	Sxx(:,i) = cpsd(Xz(:,i),Xz(:,i),fix(nt/nd),fix(q*nt/nd));
%	Syy(:,i) = cpsd(Yz(:,i),Yz(:,i),fix(nt/nd),fix(q*nt/nd));
%	Sxy(:,i) = cpsd(Xz(:,i),Yz(:,i),fix(nt/nd),fix(q*nt/nd));
%end


[ft,ftp] = freq_fft(N,time(N)-time(1));

%[X, ft] = fft_time_z(Xz,time,nd,q,tap);
%[Y, ft] = fft_time_z(Yz,time,nd,q,tap);
%
%[ndq,N,nz] = size(X);
%Sxx = zeros(ndq,N,nz);
%Syy = zeros(ndq,N,nz);
%Sxy = zeros(ndq,N,nz);
%for i=1:ndq
%	Sxx(i,:,:) = squeeze(conj(X(i,:,:))).*squeeze(X(i,:,:));
%	Syy(i,:,:) = squeeze(conj(Y(i,:,:))).*squeeze(Y(i,:,:));
%	Sxy(i,:,:) = squeeze(conj(X(i,:,:))).*squeeze(Y(i,:,:));
%end
%
%Sxx = squeeze(mean(Sxx,1));
%Syy = squeeze(mean(Syy,1));
%Sxy = squeeze(mean(Sxy,1));
%
%
%dat.time = time;
%dat.xn = xn;
%dat.yn = yn;
%dat.X = squeeze(mean(X,1));
%dat.Y = squeeze(mean(Y,1));

end

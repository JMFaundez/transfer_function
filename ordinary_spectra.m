function [Sxx,Syy,Sxy,ft,fz,dat]= ordinary_spectra(xdata,ydata,time,z0,nd,q,tap)

%NS = 31; %number of sensors in one row
NS = length(z0);
%linex = 2;
%liney = 4;

%nd = 2;
%q =0.2;
%tap = 1;

%% Load sensors measurements
% L = load('time_sensor.mat');
% A = L.A;
% 
% z0 = linspace(-0.0145,0.0145,31);
% timeR = A(2:end,1) - A(2,1);
% 
% time = linspace(0,timeR(end),floor(length(timeR)))';

%% zero padding

dt = time(2)-time(1);
% 
% nt = length(time);
% nt_n = pow2(nextpow2(nt));
% time = [time ; time(nt) + (time(nt) - time(nt-1))*[1:(nt_n-nt)]'];
% avgx = mean(xdata(:,1:end) ,1) ;
% xn = xdata;%-avgx;
% xn = [xn ;zeros(nt_n -nt,NS)];
% avgy = mean(ydata(:,1:end),1);
% yn = ydata;%-avgy;
% yn = [yn ;zeros(nt_n -nt,NS)];

nt = length(time);
%nt_n = pow2(nextpow2(nt)-1);
%time = time(nt-nt_n+1:end);
%time = time(1:nt_n);


%xdata = xdata(nt-nt_n+1:end,:);
%xdata = xdata(1:nt_n,:);
%avgx = mean(xdata(:,1:end) ,1) ;
xn = xdata;%-avgx;

%ydata = ydata(nt-nt_n+1:end,:);
%ydata = ydata(1:nt_n,:);
%avgy = mean(ydata(:,1:end),1);
yn = ydata;%-avgy;

NT = length(time);

Xz = fft(xn,[],2);
Yz = fft(yn,[],2);

[fz,fzp] = freq_fft(NS,z0(end)-z0(1));

[X, ft] = fft_time_z(Xz,time,nd,q,tap);
[Y, ft] = fft_time_z(Yz,time,nd,q,tap);

[nd,N,nz] = size(X);
Sxx = zeros(nd,N,nz);
Syy = zeros(nd,N,nz);
Sxy = zeros(nd,N,nz);
for i=1:nd
	Sxx(i,:,:) = squeeze(conj(X(i,:,:))).*squeeze(X(i,:,:));
	Syy(i,:,:) = squeeze(conj(Y(i,:,:))).*squeeze(Y(i,:,:));
	Sxy(i,:,:) = squeeze(conj(X(i,:,:))).*squeeze(Y(i,:,:));
end

Sxx = squeeze(mean(Sxx,1));
Syy = squeeze(mean(Syy,1));
Sxy = squeeze(mean(Sxy,1));


%C1 = fftshift(Sxx,1);
%C1 = fftshift(C1,2)';
%C2 = fftshift(Syy,1);
%C2 = fftshift(C2,2)';
%
%[FT,FZ] = meshgrid(sort(ft),sort(fz));
%
%figure()
%subplot(121)
%surf(FT,FZ,C1)
%view(2)
%shading interp
%colorbar()
%xlim([-20,20])
%ylim([-300,300])
%
%
%subplot(122)
%surf(FT,FZ,C2)
%view(2)
%shading interp
%colorbar()
%xlim([-20,20])
%ylim([-300,300])

dat.time = time;
dat.xn = xn;
dat.yn = yn;
dat.X = squeeze(mean(X,1));
dat.Y = squeeze(mean(Y,1));

end

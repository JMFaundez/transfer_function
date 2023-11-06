function res = ordinary_TF(y,z,zspan,time,cond)

%NS = 31; %number of sensors in one row
%liney = 2;
%linez = 3;
% if cond==1
% 	nd = 4;
% 	q =0.2;
% 	tap = 1;
% elseif cond==2
% 	nd = 2;
% 	q =0.2;
% 	tap = 1;
% elseif cond==3
% 	nd = 2;
% 	q =0.2;
% 	tap = 1;
% end
q = cond.q; %overlap
nd = cond.nd; % number of bins
tap = cond.tap; % 1 or 0

NS = length(zspan);
%zspan = linspace(-0.0145,0.0145,NS);

[Syy,Szz,Syz,ft,fz,dat]= ordinary_spectra(y,z,time,zspan,nd,q,tap);
[N,nz] = size(Syy);

TF = Syz./Syy;
gyz = ifft(TF,[],1);
gyz = ifft(gyz,[],2);

sum(gyz(1,:));
coherence = conj(Syz).*Syz./(Szz.*Syy);


[FT,FZ] = meshgrid(sort(ft),sort(fz));


%C = fftshift(coherence,1);
%C = fftshift(C,2);
%figure()
%surf(FT,FZ,C')
%view(2)
%shading interp
%colorbar()
%xlim([-30,30])
%ylim([-300,300])


%% Original Data: Convolution
fs = 2*max(ft);

timeN = time;linspace(0,(N-1)/fs,N);
%timeN = time+0.2;
y_int = y;y(1:N,:);
z_int = z;z(1:N,:);
size(y_int);
size(gyz);
N = length(time);
% y_int = interp_sens(liney,time);
% y_int = y_int - mean(y_int,1);
% z_int = interp_sens(linez,time);
% z_int = z_int - mean(z_int,1);

%%
%zconvT = conv2(real(gyz),y_int);

zconv = conv_jose(y_int,gyz,NS,N);


res.gyz = gyz;
res.coh = coherence;
res.ft = ft;
res.fz = fz;
res.zspan = zspan;
res.timeN = timeN;
res.zconv = zconv;
res.zdns = z_int;
res.Syy  = Syy;
res.Szz = Szz;
res.Syz = Syz;
res.Y = dat.X;
res.Z = dat.Y;

%name = ['gyz_l',num2str(liney),'_l',num2str(linez),'_c',num2str(cond),'.mat'];
%save(name,'res');


end


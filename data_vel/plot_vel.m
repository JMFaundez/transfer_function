clear all

L = load('Linear_vel');
iin = 1;
iout = 13;

uin = squeeze(L.u1(:,iin,:));
zin = squeeze(L.zz(iin,:));
uout = squeeze(L.u1(:,iout,:));
t = L.time;


[T,Z] = meshgrid(t,zin);

figure('Position',[300 300 1000 300])
surf(T,Z,uin')
shading interp
xlim([min(t),max(t)])
ylim([min(zin),max(zin)])
view(2)
colorbar()


figure('Position',[300 300 1000 300])
surf(T,Z,uout')
shading interp
xlim([min(t),max(t)])
ylim([min(zin),max(zin)])
view(2)
colorbar()

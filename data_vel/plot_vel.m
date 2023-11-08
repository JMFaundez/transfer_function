clear all

L = load('Linear_vel');
iin = 1;
iout = 13;

uin = squeeze(L.u1(:,iin,:));
zin = squeeze(L.z(iin,:));
t = L.time;


[T,Z] = meshgrid(t,zin);

figure()
surf(T,Z,uin')
shading interp
view(2)
colorbar()



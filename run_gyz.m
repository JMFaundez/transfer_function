clear all

% Welch parameters
cond.q = 0.75; %overlap
cond.nd = 5;  % Number of bins 
cond.tap = 1; % Use hanning window or not



%% Data Diego
tic
LD = load("../CF/Linear_Diego/wallStresses.mat",'Q','t');
toc
gridDiego = load("../CF/Linear_Diego/grid");


[Zd,Xd] = meshgrid(gridDiego.Z,gridDiego.X(:,1)); 
[Ns,Nz] = size(Xd);
time_diego = LD.t;
[Np,Nt] = size(LD.Q);


x0=0.10;
xf=0.20;
t0 = 0.6;
ix0 = find(Xd(:,1)>=x0,1,"first");
ixf = find(Xd(:,1)>=xf,1,"first");
it0 = find(time_diego>=t0,1,"first");
time = time_diego(it0:end);
nnt = floor(length(time)/cond.nd) * cond.nd;
dift = length(time) - nnt;
it0 = it0 + dift;
time = time(dift+1:end);

cf = (LD.Q(Np/3+1:2*Np/3,:));
cf = reshape(cf,Ns,Nz,Nt);
I = squeeze(cf(ix0,1:end-1,it0:end))';
O = squeeze(cf(ixf,1:end-1,it0:end))';

z = gridDiego.Z(1:end-1);

tic
out = gyz_1row(I,O,cond,time,z);
toc


[Z,T] = meshgrid(z,time);
val=0.5e-3;
cmap=parula;
figure('Position',[500 500 1000 400])
ax1 =subplot(211);
pcolor(T,Z,O)
shading interp
clim([-val,val])
colormap(ax1,cmap)
colorbar()
%axis equal
ylim([min(Z(:)),max(Z(:))])
%xlim([0.02,0.5])
title('Output')
ylabel('$z$','FontSize',18,'Interpreter','Latex')
xlabel('time','FontSize',18,'Interpreter','Latex')

ax2 =subplot(212);
pcolor(T,Z,out.est)
shading interp
clim([-val,val])
colormap(ax1,cmap)
colorbar()
%axis equal
ylim([min(Z(:)),max(Z(:))])
%xlim([0.02,0.5])
title('Prediction')
ylabel('$z$','FontSize',18,'Interpreter','Latex')
xlabel('time','FontSize',18,'Interpreter','Latex')

iz = 50;
figure()
hold on 
plot(time,O(:,iz),'DisplayName','Sim')
plot(time,out.est(:,iz),'DisplayName','Prediction')
ylabel('output','FontSize',18,'Interpreter','Latex')
xlabel('time','FontSize',18,'Interpreter','Latex')
legend('Location','best')
box on 



[FT,FZ] = meshgrid(fftshift(out.ft)*2*pi,fftshift(out.fz)*2*pi);
figure()
surf(FT,FZ,fftshift(out.coherence'))
shading interp
view(2)
clim([0,1])
xlim([0,400])
ylim([-1000,1000])
xlabel('$\omega$','FontSize',18,'Interpreter','Latex')
ylabel('$\beta$','FontSize',18,'Interpreter','Latex')
title('Coherence','FontSize',18,'Interpreter','Latex')
colorbar()

[~,ntg] = size(FT);
[Tg,Zg] = meshgrid(time(1:ntg)-time(1),z);

figure()
surf(Tg,Zg,real(fftshift(out.gyz))')
shading interp
view(2)
%xlim([0,400])
%ylim([-1000,1000])
xlabel('time','FontSize',18,'Interpreter','Latex')
ylabel('span','FontSize',18,'Interpreter','Latex')
title('TF','FontSize',18,'Interpreter','Latex')
colorbar()



%val=0.5e-3;
%cmap=parula;
%figure('Position',[500 500 1000 400])
%ax1 =subplot(211);
%pcolor(T,Z,I)
%shading interp
%clim([-val,val])
%colormap(ax1,cmap)
%colorbar()
%%axis equal
%ylim([min(Z(:)),max(Z(:))])
%%xlim([0.02,0.5])
%title('Input Linear')
%ylabel('$z$','FontSize',18,'Interpreter','Latex')
%xlabel('time','FontSize',18,'Interpreter','Latex')
%
%ax2 =subplot(212);
%pcolor(T,Z,O)
%shading interp
%clim([-val,val])
%colormap(ax1,cmap)
%colorbar()
%%axis equal
%ylim([min(Z(:)),max(Z(:))])
%%xlim([0.02,0.5])
%title('Output Linear')
%ylabel('$z$','FontSize',18,'Interpreter','Latex')
%xlabel('time','FontSize',18,'Interpreter','Latex')
%
%cfLD = squeeze(cf(:,:,600));
%figure('Position',[500 500 1000 400])
%ax1 =subplot(211);
%pcolor(Xd,Zd,cfLD)
%shading interp
%clim([-val,val])
%colormap(ax1,cmap)
%colorbar()
%axis equal
%ylim([min(Z(:)),max(Z(:))])
%xlim([0.02,0.5])
%ylabel('$z$','FontSize',18,'Interpreter','Latex')
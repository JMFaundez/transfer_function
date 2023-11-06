clear all

% Welch parameters
cond.q = 0.75; %overlap
cond.nd = 10;  % Number of bins 
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


x0=0.05;
xf=0.1;
ix0 = find(Xd(:,1)>=x0,1,"first");
ixf = find(Xd(:,1)>=xf,1,"first");
it0 = find(time_diego>=0.3,1,"first");
time = time_diego(it0:end);
nnt = floor(length(time)/cond.nd) * cond.nd;
it0 = it0 + length(time)-nnt;


cf = (LD.Q(Np/3+1:2*Np/3,:));
cf = reshape(cf,Ns,Nz,Nt);
I = squeeze(cf(ix0,1:end-1,it0:end))';
O = squeeze(cf(ixf,1:end-1,it0:end))';

z = gridDiego.Z(1:end-1);

out = gyz_1row(I,O,cond,time,z);


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
title('Input Linear')
ylabel('$z$','FontSize',18,'Interpreter','Latex')
xlabel('time','FontSize',18,'Interpreter','Latex')

ax2 =subplot(212);
pcolor(T,Z,out.z_s)
shading interp
clim([-val,val])
colormap(ax1,cmap)
colorbar()
%axis equal
ylim([min(Z(:)),max(Z(:))])
%xlim([0.02,0.5])
title('Output Linear')
ylabel('$z$','FontSize',18,'Interpreter','Latex')
xlabel('time','FontSize',18,'Interpreter','Latex')

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
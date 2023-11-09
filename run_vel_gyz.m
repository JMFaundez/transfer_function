clear all

% Welch parameters
cond.q = 0.50; %overlap
cond.nd = 5;  % Number of bins 
cond.tap = 1; % Use hanning window or not



%% Data Diego
tic
LD = load("data_vel/Linear_vel.mat");
toc

%[Ns,Nz] = size(Xd);
time = LD.time;
[Nt,Nprobes,nz] = size(LD.u1);



ix0 = [8];
ixf = [13];

t0 = 0.6;
it0 = find(time>=t0,1,"first");
time = time(it0:end);
nnt = floor(length(time)/cond.nd) * cond.nd;
dift = length(time) - nnt;
it0 = it0 + dift;
time = time(dift+1:end);

O = squeeze(LD.u1(it0:end,ixf,:));

[n1,n2] = size(O);
nvel=3;
I = zeros(n1,n2,length(ix0)*nvel);
for i=1:length(ix0)
    if nvel==1
        Ii = squeeze(LD.u1(it0:end,ix0(i),:));
        I(:,:,i) = Ii;
    elseif nvel==2 
        Ii = squeeze(LD.u1(it0:end,ix0(i),:));
        Iii = squeeze(LD.u2(it0:end,ix0(i),:));
        I(:,:,(i-1)*nvel+1) = Ii;
        I(:,:,(i-1)*nvel+2) = Iii;
    elseif nvel==3 
        Ii = squeeze(LD.u1(it0:end,ix0(i),:));
        Iii = squeeze(LD.u2(it0:end,ix0(i),:));
        Iiii = squeeze(LD.u3(it0:end,ix0(i),:));
        I(:,:,(i-1)*nvel+1) = Ii;
        I(:,:,(i-1)*nvel+2) = Iii;
        I(:,:,(i-1)*nvel+3) = Iiii;
    end
end

z = LD.zz(ixf,:);

tic
%out = gyz_1row(I,O,cond,time,z);
out = gyz_nrows(I,O,cond,time,z);
toc
%return
%disp("Error="+num2str(out.error))
figure()
plot(out.error)

iz = 50;
[Z,T] = meshgrid(z,time);
val=0.08;
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
xlim([min(T(:)),max(T(:))])
title('Output')
ylabel('$z$','FontSize',18,'Interpreter','Latex')
xlabel('time','FontSize',18,'Interpreter','Latex')

ax2 =subplot(212);
hold on
pcolor(T,Z,out.est)
shading interp
clim([-val,val])
colormap(ax1,cmap)
colorbar()
yline(z(iz))
%axis equal
ylim([min(Z(:)),max(Z(:))])
xlim([min(T(:)),max(T(:))])
title('Prediction')
ylabel('$z$','FontSize',18,'Interpreter','Latex')
xlabel('time','FontSize',18,'Interpreter','Latex')

figure()
hold on 
plot(time,O(:,iz),'DisplayName','Sim')
plot(time,out.est(:,iz),'DisplayName','Prediction')
%for i=1:length(ix0)
%    plot(time,out.z_s(:,iz,i),'DisplayName',"$x_{"+num2str(i)+"}$")
%end
ylabel('output','FontSize',18,'Interpreter','Latex')
xlabel('time','FontSize',18,'Interpreter','Latex')
legend('Location','best','Interpreter','latex')
box on 



%[FT,FZ] = meshgrid(fftshift(out.ft)*2*pi,fftshift(out.fz)*2*pi);
ft_shift = fftshift(out.ft);
fz_shift  = fftshift(out.fz);
dft = ft_shift(2)-ft_shift(1);
dfz = fz_shift(2)-fz_shift(1);
xx = (ft_shift - dft/2)*2*pi;
yy = (fz_shift - dfz/2)*2*pi;
[FT,FZ] = meshgrid(xx,yy);

figure()
hold on
imagesc(ft_shift*2*pi,fz_shift*2*pi,fftshift(out.coherence'))
%shading interp
view(2)
clim([0,1])
%mesh(FT,FZ,FZ*0,'FaceAlpha',0,'LineWidth',1,'EdgeColor','k')
xlim([0,floor(400/dft)*dft])
ylim([0,floor(1000/dfz)*dfz])
xlabel('$\omega$','FontSize',18,'Interpreter','Latex')
ylabel('$\beta$','FontSize',18,'Interpreter','Latex')
title('Coherence','FontSize',18,'Interpreter','Latex')
colorbar()

[~,ntg] = size(FT);
[Tg,Zg] = meshgrid(time(1:ntg)-time(1),z);

figure()
surf(Tg,Zg,real((out.gyz(:,:,1)))')
shading interp
view(2)
%xlim([0,400])
%ylim([-1000,1000])
xlabel('time','FontSize',18,'Interpreter','Latex')
ylabel('span','FontSize',18,'Interpreter','Latex')
title('TF','FontSize',18,'Interpreter','Latex')
colorbar()
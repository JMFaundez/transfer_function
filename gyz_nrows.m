function out = gyz_nrows(I,O,cond,t,z0)
%% Parameters
%dtdns = 1e-5;
%z0 = linspace(-0.035,0.035,21);
%z0 = z0(1:end-1);
NS = length(z0);
[nt,nz,nin] = size(I);


%nd = 8;  q=0.75;  tap =1;
%c1 = 1e-4;
%tresh = 0.1e-3;
w1 = 100;
%noise_opt = 2; %1: my old version, 2: With weigths, Jose step function. 3:With weigths, Diego step function.
%liney1 = 2; % x=0.025,0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2
%liney2 = 3;
%linez = 6;
%% Load uncontrolled data

% Ls = load('timesignals_unc.mat');
% A = Ls.A;
% tsensor = Ls.t;
% B = load('ref_val.mat').Bref;
% yunc1 = A(1:end,NS*(liney1-1)+1:NS*liney1)-B(1,NS*(liney1-1)+1:NS*liney1);%A(1,NS*(liney-1)+1:NS*liney);
% yunc2 = A(1:end,NS*(liney2-1)+1:NS*liney2)-B(1,NS*(liney2-1)+1:NS*liney2);
% zunc = A(1:end,NS*(linez-1)+1:NS*linez)-B(1,NS*(linez-1)+1:NS*linez);%A(1,NS*(linez-1)+1:NS*linez);
% 
% %%
% % downsampling and choose a nice N
% iin = 3;
% step = 100;
% ys1 = yunc1(iin:step:end,:); ys2 = yunc2(iin:step:end,:); 
% zs = zunc(iin:step:end,:);
% 
% times = tsensor(iin:step:end);
% 
% cut = mod(length(times),nd);
% ys1 = ys1(cut+1:end,:);
% ys2 = ys2(cut+1:end,:);
% zs = zs(cut+1:end,:);
% times = times(cut+1:end);

%cond.q = q; cond.nd = nd; cond.tap = tap;


%% Compute Gyz

res_in = {};
res_in_out = {};
for i=1:nin
    in1 = squeeze(I(:,:,i));
    res_in_out{i} = ordinary_TF(in1,O,z0,t,cond);
    for j=1:nin
        in2 = squeeze(I(:,:,j));
        res_in{i,j} = ordinary_TF(in1,in2,z0,t,cond);
    end
end


[ntf,nz] = size(res_in{1,1}.Syy);

%A = zeros(ntf,nz,nin,nin);

Gyz = zeros(ntf,nz,nin);
A = zeros(ntf,nz,nin,nin);
b = zeros(ntf,nz,nin);

for i=1:nin
    b(:,:,i) = res_in_out{i}.Syz;
    for j=1:nin
        A(:,:,i,j) = res_in{i,j}.Syz;
    end
end

for i=1:ntf
    for j=1:nz
        Ai = squeeze(A(i,j,:,:));
        bi = squeeze(b(i,j,:));
        Gyz(i,j,:) = Ai\bi;
    end
end

gyz = ifft(Gyz,[],1);
gyz = ifft(gyz,[],2);

z_s = zeros(length(t),nz,nin);

for i=1:nin
    gyz_i = real(squeeze(gyz(:,:,i))); 
    in1 = squeeze(I(:,:,i));
    z_s(:,:,i) = conv_jose(in1,gyz_i,NS,length(t));
end

out.gyz = real(gyz);
out.z_s = z_s;
out.ft = res_in_out{1}.ft;
out.fz = res_in_out{1}.fz;
out.coherence = res_in_out{1}.coh;



%%

%Sy1y1 = res3.Syy;
%Sy2y2 = res3.Szz;
%Sy1y2 = res3.Syz;
%Sy2y1 = res4.Syz;
%
%Sy1z = res1.Syz;
%Sy2z = res2.Syz;
%
%% Gyz =  | Gyz_1   Gyz_2 |
%Sx = step_2(res1.ft*2*pi,w1,200,1e12,1e4);
%Sx = ifftshift(Sx);
%R2 = min(Sx);
%SxM = repmat(Sx',[1 20]);
%SxM = SxM/R2;
%SxM_inv = 1./SxM;
%
%
%
%%Gyz_1 = (Sy1z./Sy1y1 + Sy2z./Sy1y2).*SxM_inv;
%%Gyz_2 = (Sy1z./Sy2y1 + Sy2z./Sy2y2).*SxM_inv;
%
%Gyz_s1 = Sy1z./Sy1y1.*SxM_inv;
%Gyz_s2 = Sy2z./Sy2y2.*SxM_inv;
%
%
%gyz_s1 = ifft(Gyz_s1,[],1);
%gyz_s1 = ifft(gyz_s1,[],2);
%
%gyz_s2 = ifft(Gyz_s2,[],1);
%gyz_s2 = ifft(gyz_s2,[],2);
%
%det = Sy1y1.*Sy2y2 - Sy1y2.*Sy2y1;
%Gyz_1 = (Sy2y2.*Sy1z - Sy1y2.*Sy2z)./det.*SxM_inv; 
%Gyz_2 = (-Sy2y1.*Sy1z + Sy1y1.*Sy2z)./det.*SxM_inv; 
%
%
%gyz_1 = ifft(Gyz_1,[],1);
%gyz_1 = ifft(gyz_1,[],2);
%
%gyz_2 = ifft(Gyz_2,[],1);
%gyz_2 = ifft(gyz_2,[],2);
%
%
%z1 = conv_jose(ys1,real(gyz_1),NS,length(times));
%z2 = conv_jose(ys2,real(gyz_2),NS,length(times));
%
%zest = z1+z2;
%
%
%z_s1 = conv_jose(ys1,real(gyz_s1),NS,length(times));
%z_s2 = conv_jose(ys2,real(gyz_s2),NS,length(times));
%
%%%
%
%
%% zi = 11;
%% figure('Position',[500 500 1200 900])
%% 
%% subplot(311)
%% hold on
%% plot(times, z_s1(:,zi),'g','DisplayName','$z_{y_1}$')
%% plot(times, zs(:,zi),'k','DisplayName','$z_{meas}$')
%% box on
%% legend('Interpreter','latex','FontSize',18)
%% xlabel('time','FontSize',18,'Interpreter','latex')
%% ylabel('$z$','FontSize',18,'Interpreter','latex')
%% xlim([0,max(times)])
%% 
%% subplot(312)
%% hold on
%% plot(times, z_s2(:,zi),'b','DisplayName','$z_{y_2}$')
%% plot(times, zs(:,zi),'k','DisplayName','$z_{meas}$')
%% box on
%% legend('Interpreter','latex','FontSize',18)
%% xlabel('time','FontSize',18,'Interpreter','latex')
%% ylabel('$z$','FontSize',18,'Interpreter','latex')
%% xlim([0,max(times)])
%% 
%% subplot(313)
%% hold on
%% plot(times, zest(:,zi),'r','DisplayName','$z_{y_1+y_2}$')
%% plot(times, zs(:,zi),'k','DisplayName','$z_{meas}$')
%% box on
%% legend('Interpreter','latex','FontSize',18)
%% xlabel('time','FontSize',18,'Interpreter','latex')
%% ylabel('$z$','FontSize',18,'Interpreter','latex')
%% xlim([0,max(times)])
%
%
%
%error1 = mse(z_s1-zs)/mse(zs);
%error2 = mse(z_s2-zs)/mse(zs);
%error12 = mse(zest-zs)/mse(zs);
%
%out.error = [error1,error2,error12];
%out.gyz_1 = gyz_1;
%out.gyz_2 = gyz_2;
%out.est = zest;
%

end
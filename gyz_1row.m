function out = gyz_1row(I,O,cond,t,z0)
%% Parameters
%dtdns = 1e-5;
%z0 = linspace(-0.035,0.035,21);
%z0 = z0(1:end-1);
NS = length(z0);

w1 = 500;

%% Compute Gyz
res1 = ordinary_TF(I,O,z0,t,cond);

%%
Syy = res1.Syy;

Syz = res1.Syz;

% Gyz =  | Gyz_1   Gyz_2 |
Sx = step_2(res1.ft*2*pi,w1,600,1e12,1e4);
Sx = ifftshift(Sx);
R2 = min(Sx);
SxM = repmat(Sx',[1 NS]);
SxM = SxM/R2;
SxM_inv = 1./SxM;
figure()
plot(SxM_inv)

%Gyz_1 = (Sy1z./Sy1y1 + Sy2z./Sy1y2).*SxM_inv;
%Gyz_2 = (Sy1z./Sy2y1 + Sy2z./Sy2y2).*SxM_inv;

%Gyz = Syz./Syy.*SxM_inv;
Gyz = Syz./Syy;

gyz = ifft(Gyz,[],1);
gyz = ifft(gyz,[],2);


z_s = conv_jose(I,real(gyz),NS,length(t));

error = mse(z_s-O)/mse(O);

out.error = error;
out.est = z_s;
out.gyz = gyz;
out.Gyz = Gyz;
out.coherence = res1.coh;
out.ft = res1.ft;
out.fz = res1.fz;

end

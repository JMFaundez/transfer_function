function out = gyz_nrows(I,O,cond,t,z0)
%% Parameters
%dtdns = 1e-5;
%z0 = linspace(-0.035,0.035,21);
%z0 = z0(1:end-1);
NS = length(z0);
[nt,nz,nin] = size(I);

nd = cond.nd;
q = cond.q;
tap = cond.tap;
w1 = 100;

%% Compute Gyz

SII = {};
SIO = {};
for i=1:nin
    in1 = squeeze(I(:,:,i));
    %res_in_out{i} = ordinary_TF(in1,O,z0,t,cond);
    [~,SOO,SIO{i},ft,fz,~] = ordinary_spectra(in1,O,t,z0,nd,q,tap);
    %[~,SOO,SIO{i},ft,fz] = ordinary_spectra2(in1,O,t,z0,nd,q);
    for j=1:nin
        in2 = squeeze(I(:,:,j));
        %res_in{i,j} = ordinary_TF(in1,in2,z0,t,cond);
        [~,~,SII{i,j},ft,fz,~] = ordinary_spectra(in1,in2,t,z0,nd,q,tap);
        %[~,~,SII{i,j},ft,fz] = ordinary_spectra2(in1,in2,t,z0,nd,q);
    end
end


[ntf,nz] = size(SII{1,1});

%A = zeros(ntf,nz,nin,nin);

Gyz = zeros(ntf,nz,nin);
A = zeros(ntf,nz,nin,nin);
b = zeros(ntf,nz,nin);

Sx = step_2(ft*2*pi,300,500,1e12,1e4);
Sx = ifftshift(Sx);
R2 = min(Sx);
SxM = repmat(Sx',[1 NS]);
SxM = SxM/R2;
SxM_inv = 1./SxM;

for i=1:nin
    b(:,:,i) = SIO{i};
    for j=1:nin
        A(:,:,i,j) = SII{i,j};
    end
end

for i=1:ntf
    for j=1:nz
        Ai = squeeze(A(i,j,:,:));
        bi = squeeze(b(i,j,:));
        Gyz_ij = Ai\bi;
        Gyz(i,j,:) = Gyz_ij;

    end
end

%for i=1:nin
%    Gyz(:,:,i) = Gyz(:,:,i).*SxM_inv;
%end

%size(Ai)
%size(bi)

gyz = ifft(Gyz,[],1);
gyz = ifft(gyz,[],2);

z_s = zeros(length(t),nz,nin);

for i=1:nin
    gyz_i = real(squeeze(gyz(:,:,i))); 
    in1 = squeeze(I(:,:,i));
    z_s(:,:,i) = conv_jose(in1,gyz_i,NS,length(t));
end
estimation = squeeze(sum(z_s,3));
error = rms(estimation-O,1)./rms(O,1);
coherence1 = conj(SIO{1}).*SIO{1}./(SOO.*SII{1,1});

out.gyz = real(gyz);
out.z_s = z_s;
out.est = estimation; 
out.ft = ft;
out.fz = fz;
out.coherence = coherence1;
out.error = error;

end
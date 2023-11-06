function S = step_2(omega,omax,odm,v1,v2)
%omega : vector of frequencies 0,...,-f1
%omax = 150; start of damping
%odm = 300; end of damping
%v1 = 1e14; value after odm
%v2 = 1e3; % value before omax
omegas = fftshift(omega);
ostar = 1-(omegas-odm)/(omax-odm);
S = (v1-v2)./(1+exp(1./(ostar-1) +1./ostar))+v2;
S(ostar<=0) =v2;
S(ostar>=1)=v1;
io0 = find(omegas>=0,1,'first');
nom = length(omega);
if mod(nom,2)==0
    S(2:io0-1) = flip(S(io0+1:end));
    S(1) = S(2);
else
    S(1:io0-1) = flip(S(io0+1:end));
    %Sx(1) = Sx(2);
end

end
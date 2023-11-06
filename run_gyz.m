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

cf = (LD.Q(Np/3+1:2*Np/3,:));
cf = reshape(cf,Ns,Nz,Nt);


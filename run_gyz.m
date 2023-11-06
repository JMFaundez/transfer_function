clear all

% Welch parameters
cond.q = 0.75; %overlap
cond.nd = 10;  % Number of bins 
cond.tap = 1; % Use hanning window or not



%% Data Diego
tic
LD = load("../Linear_Diego/wallStresses.mat",'Q','t');
toc
gridDiego = load("../CF/Linear_Diego/grid");


[Zd,Xd] = meshgrid(gridDiego.Z,gridDiego.X(:,1)); 
[Ns,Nz] = size(Xd);
time_diego = LD.t;
[Np,Nt] = size(LD.Q);
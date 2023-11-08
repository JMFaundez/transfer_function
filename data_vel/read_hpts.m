clear all

casei=1;
switch casei
case 1
	fnameopen = "FST.his";
	filename = "Linear_vel";
	time_steps = 480+480+480+480+480;
	iinit=0;
	nxLE=3;
	nyLE=4;
	nxBL=4;
	nyBL=3;
	nz = 150;
end


disp('================== STARTING PROGRAM ======================')
disp('Opening file...')
fileID = fopen(fnameopen,'r');
A = fscanf(fileID,'%f');
disp('file open!')
disp(fnameopen)
n_fields = 4;
npoints = A(1);
coord =  A(2:npoints*3+1);
xyz = zeros(npoints,3);
xyz(:,1) = coord(1:3:end);
xyz(:,2) = coord(2:3:end);
xyz(:,3) = coord(3:3:end);

NLE = nxLE*nyLE;
NBL = nxBL*nyBL;
nprobes=NLE+NBL;
%nprobes = nx*ny
%nz = npoints/nprobes
%nplane = nz*ny


%time_steps = 161;33; (length(A)-3*npoints-1)/(npoints*(n_fields+1));

fields = A(3*npoints+2:end);
clear A
time = fields(1:(n_fields+1):end);
vx = fields(2:(n_fields+1):end);
vy = fields(3:(n_fields+1):end);
vz = fields(4:(n_fields+1):end);
pr = fields(5:(n_fields+1):end);
nt = time_steps;

u1 = zeros(nt,nprobes,nz); u2 = zeros(nt,nprobes,nz); u3 = zeros(nt,nprobes,nz);
xx = zeros(nprobes,nz); yy = zeros(nprobes,nz); zz = zeros(nprobes,nz);
p = zeros(nt,nprobes,nz);

%u1 = zeros(nt,nx,ny,nz); u2 = zeros(nt,nx,ny,nz); u3 = zeros(nt,nx,ny,nz);
%xx = zeros(nx,ny,nz); yy = zeros(nx,ny,nz); zz = zeros(nx,ny,nz);
%p = zeros(nt,nx,ny,nz);
for t=1:nt
	t
for i=1:nprobes
	%index = (t-1)*nplane*nx+nplane*(i-1)+1:(t-1)*nplane*nx+nplane*i;
	index = (t-1)*nz*nprobes+nprobes*(i-1)+1:(t-1)*nz*nprobes+nprobes*i;
	if t==1
		inn=[index(1),index(end),index(end)-index(1)];
	end
	u1(t,i,:) = reshape(vx(index),nz,ny)';	
	u2(t,i,:) = reshape(vy(index),nz,ny)';	
	u3(t,i,:) = reshape(vz(index),nz,ny)';
  	p(t,i,:) = reshape(pr(index),nz,ny)';

%	v1 = reshape(vx(index),nz,ny)';	
%	v2 = reshape(vy(index),nz,ny)';	
%	v3 = reshape(vz(index),nz,ny)';
%  	p = reshape(pr(index),nz,ny)';
	%fnamei = filename+num2str(i)+"_it_"+num2str(t+iinit)
	%save(fnamei,'v1','v2','v3','p','-v7.3')
	%disp("saving file "+fnamei)
end
end

%disp('Saving files...')

%xx = zeros(nprobes,nz);
%yy = zeros(nprobes,nz);
%zz = zeros(nprobes,nz);

for i=1:nprobes
	index = nz*(i-1)+1:nz*i;
	xx(i,:) = reshape(xyz(index,1),nz,ny)';	
	yy(i,:) = reshape(xyz(index,2),nz,ny)';	
	zz(i,:) = reshape(xyz(index,3),nz,ny)';	
%	save(filename+num2str(i)+"_mesh",'xx','yy','zz')
%	for it=1:nt
%	v1 = squeeze(u1(it,i,:,:));
%	v2 = squeeze(u2(it,i,:,:));
%	v3 = squeeze(u3(it,i,:,:));
%	pr = squeeze(p(it,i,:,:));
%	fnamei = filename+num2str(i)+"_it_"+num2str(it+iinit);
%	save(fnamei,'v1','v2','v3','pr','-v7.3')
%	disp("saving file "+fnamei)
% 	end
end

%%

save(filename,'u1','v2','u3','p','xx','yy','zz','-v7.3')
disp('files SHOULD be saved...')
disp('================== END PROGRAM ======================')

%con = mean(u1,3;
%X = xx(:,:,20);
%Y = yy(:,:,20);

%figure(100)
%contourf(X,Y,con)
%axis('equal')
%colorbar()

% t,vx,vy,vz,p
%to_file = zeros(time_steps,n_fields+1,nz);
%tic
%for i=1:nprobes
%    for k=1:nz
%        to_file(:,1,k) = time(k+(i-1)*nz:npoints:end);
%        to_file(:,2,k) = vx(k+(i-1)*nz:npoints:end);
%        to_file(:,3,k) = vy(k+(i-1)*nz:npoints:end);
%        to_file(:,4,k) = vz(k+(i-1)*nz:npoints:end);
%        to_file(:,5,k) = pr(k+(i-1)*nz:npoints:end);
%    end
%    filename = ['./probes/probe_',num2str(i),'.mat'];
%    save(filename,'to_file')
%    XYZ = xyz((i-1)*nz+1:i*nz,:);
%    filename = ['./probes/xyz_',num2str(i),'.mat'];
%    save(filename,'XYZ')
%end
%toc

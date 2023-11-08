function read_hpts(casei)
%clear all

switch casei
case 1
	fnameopen = "FST_BF.his";
	filename = "mesh1_BF/FST25_BF_";
	time_steps = 1;
	iinit=0;
	nx=58;
	ny=50;
case 2
	fnameopen = "FST_unc.his";
	filename = "mesh1_unc/FST25_unc_";
	time_steps = 161+160;
	iinit=0;
	nx=58;
	ny=50;
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


nprobes = nx*ny
nz = npoints/nprobes
nplane = nz*ny


%time_steps = 161;33; (length(A)-3*npoints-1)/(npoints*(n_fields+1));

fields = A(3*npoints+2:end);
clear A
time = fields(1:(n_fields+1):end);
vx = fields(2:(n_fields+1):end);
vy = fields(3:(n_fields+1):end);
vz = fields(4:(n_fields+1):end);
pr = fields(5:(n_fields+1):end);
nt = time_steps


%u1 = zeros(nt,nx,ny,nz); u2 = zeros(nt,nx,ny,nz); u3 = zeros(nt,nx,ny,nz);
%xx = zeros(nx,ny,nz); yy = zeros(nx,ny,nz); zz = zeros(nx,ny,nz);
%p = zeros(nt,nx,ny,nz);
for t=1:nt
	t
for i=1:nx
	index = (t-1)*nplane*nx+nplane*(i-1)+1:(t-1)*nplane*nx+nplane*i;
	if t==1
		inn=[index(1),index(end),index(end)-index(1)];
	end
%	u1(t,i,:,:) = reshape(vx(index),nz,ny)';	
%	u2(t,i,:,:) = reshape(vy(index),nz,ny)';	
%	u3(t,i,:,:) = reshape(vz(index),nz,ny)';
%  p(t,i,:,:) = reshape(pr(index),nz,ny)';

	v1 = reshape(vx(index),nz,ny)';	
	v2 = reshape(vy(index),nz,ny)';	
	v3 = reshape(vz(index),nz,ny)';
  p = reshape(pr(index),nz,ny)';
	fnamei = filename+num2str(i)+"_it_"+num2str(t+iinit)
	save(fnamei,'v1','v2','v3','p','-v7.3')
	disp("saving file "+fnamei)
end
end

%disp('Saving files...')

for i=1:nx
	index = nplane*(i-1)+1:nplane*i;
	xx = reshape(xyz(index,1),nz,ny)';	
	yy = reshape(xyz(index,2),nz,ny)';	
	zz = reshape(xyz(index,3),nz,ny)';	
	save(filename+num2str(i)+"_mesh",'xx','yy','zz')
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
end

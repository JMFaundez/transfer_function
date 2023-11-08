function zconv = conv_jose(y_int,gyz,NS,N)
    %zconv = zconvT(1:N,1:nz);
    %z2side = zconvT(1:N,nz+1:end);
    %zconv(:,1:nz-1) = zconv(:,1:nz-1) +  z2side;
    nz = NS;
    ny = nz;
    [ngyz,~] = size(gyz);
    nint = floor(ngyz/1);
    zconvT = zeros(N,nz,ny);
    for k=1:nz
        for m=1:ny
            ind = m + k - 1;
            if ind>ny
                ind2 = m + k - 1 -ny;
                zconvi = conv(y_int(:,ind2),real(gyz(1:nint,m)));
            else
                zconvi = conv(y_int(:,ind),real(gyz(1:nint,m)));
            end
            zconvT(:,k,m) = zconvi(1:N);
        end
    end
    zconv = sum(zconvT,3);
    
    %dif = abs(zconv-z_int);
    %error = mean(dif.^2,'all');
end
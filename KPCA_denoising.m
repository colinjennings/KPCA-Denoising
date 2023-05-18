function out = KPCA_denoising(dwi,mask,kernel,flag_shrinkage)
%[out,noisemap] = KPCA_denoising(dwi,mask,kernel,flag_shrinkage)

% Returns out (denoised image) and noisemap (estimated noise map).

% flag_shrinkage : 1, apply shrinkage, 0 normal KPCA
% get biased estimates. 

%G.Ramos Llorden 1st August 2021

kernel = kernel + (mod(kernel, 2)-1);   % needs to be odd.
k = (kernel-1)/2; kx = k(1); ky = k(2); kz = k(3);
center = sub2ind([kernel(1),kernel(2),kernel(3)],kx+1,ky+1,kz+1);
N = prod(kernel);

[nx,ny,nz,M]=size(dwi);
out = dwi;

[cSure,hSure] = optimal_sigma_rank_SURE(dwi,mask,kernel,0); %estimation of optimal c and h 


rangey = (ky+1):(ny-ky);
rangez = (kz+1):(nz-kz);
rangex = (kx+1):(nx-kx);

%denoising begins there, we use parfor, remove_if_errors.

parfor nxx=rangex
    for nzz=rangez
              for nyy=rangey
                if mask(nxx,nyy,nzz)~=0
                    X = dwi(nxx-kx:nxx+kx, nyy-ky:nyy+ky, nzz-kz:nzz+kz,:);
                    X = reshape(X, N, M); X = X';       
                    c= cSure(nxx,nyy,nzz),
                    h = hSure(nxx,nyy,nzz);
                    [ Signal ] = KPCA_per_block(X, h, c, center,flag_shrinkage);
                    out(nxx,nyy,nzz,:) = Signal;
                else
                    out(nxx,nyy,nzz,:)=0;
                end
             end
     end
end

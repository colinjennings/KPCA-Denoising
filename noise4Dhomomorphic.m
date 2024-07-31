function Noisemap = noise4Dhomomorphic(Im,SNR,LPF,Mode)

addpath(genpath('./utils'))
addpath(genpath('./utils/noise_map_homomorphic'))

[nx,ny,nz,M]=size(Im);
h2D = fspecial('gaussian',[nx,ny],LPF);
h = fspecial('gaussian',[nz,1],LPF);

for nzz=1:nz
    h3D(:,:,nzz) = h2D*h(nzz);
end 
   
h3D=h3D./max(h3D(:));
save('LPF.mat','h3D');


for ndir = 1:M
    [MapGaussian]=rice_homomorf_est3D(squeeze(Im(:,:,:,ndir)),SNR,LPF,Mode);
    Noisemap(:,:,:,ndir) = MapGaussian;
end


system('rm LPF.mat')

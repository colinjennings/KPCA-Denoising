function [cSure,hSure,SURE,Noisemap] = optimal_sigma_rank_SURE(Im,mask,kernel_estim,flag_shrinkage)
 %[cSure,hSure] = optimal_sigma_rank_SURE(Im,mask,kernel_estim,flag_shrinkage)
 % GRamos Llorden et al , MGH Martinos center
addpath(genpath('./utils'))

Noisemap = noise4Dhomomorphic(Im,0,2,2);
%noise = min(Noisemap,[],4);
noise = mean(Noisemap,4);
kernel_estim = kernel_estim + (mod(kernel_estim, 2)-1);   % needs to be odd.


% First, we divide image into non-overlapping blocks

[Imcell,~,~,~] = divide4Dimageinsquares(Im,kernel_estim,mask);

[Noisecell,icenter,jcenter,kcenter] = divide4Dimageinsquares(noise,kernel_estim,mask);
Tot = numel(icenter(:));
[nx,ny,nz,M]=size(Im);

%parforsetup(50); %use 50 % of cores

c = logspace(log10(0.6),log10(6),5);
h = unique(round(logspace(log10(4),log10(40),10)));

%We process each block separately
Hopt_all=cell(Tot,1);
Copt_all = Hopt_all;

parfor nn=1:Tot
    data_cell = Imcell{nn};
    noise_cell = mean(Noisecell{nn}(:));
    [h_cell, c_cell,sure] = sureKPCA_shrinkage( data_cell, h, c, noise_cell, flag_shrinkage);  %actual optimization happens here
    Hopt_all{nn} = h_cell;
    Copt_all{nn} = c_cell;
    SURE{nn}=sure;
    %fprintf('Progress %d %% \n:', round(100*nn/Tot))
end

% Final maps
cSure = unpatch_blocks4D(Copt_all,icenter,jcenter,kcenter, kernel_estim,[nx,ny,nz]);
hSure = unpatch_blocks4D(Hopt_all,icenter,jcenter,kcenter, kernel_estim,[nx,ny,nz]);
save('maps.mat','cSure','hSure','-v7.3')

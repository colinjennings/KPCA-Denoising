function [hopt,copt,SURE ] = sureKPCA_shrinkage( data_cell, h, c, sigma_noise, flag_shrinkage )

%[1] S. Ramani, T. Blu, M. Unser, "Monte-Carlo Sure: A Black-Box Optimization of Regularization Parameters
%for General Denoising Algorithms", IEEE TIP 2008

hvector = h; %grid of h
cvector = c; %grid of c
P = numel(data_cell(:));
P_max = size(data_cell,2);
SURE = zeros(numel(cvector),numel(hvector));
eps_scale = 1e-2;


%Brute force search to find the minimum of the Monte-Carlo SURE 

for cc=1:numel(cvector)
    for hh=1:numel(hvector)
            c = cvector(cc);
            h = hvector(hh);
            h=min(h,P_max);
            denoised_cell = vect_KPCA_shrinkage( data_cell, h, c, flag_shrinkage); %output of KPCA denoising of signal.   f_lambda(y) Eq.9 of [1]
            denoised_cell = denoised_cell(:);
            %% Sure_random 
            b =randn(size(data_cell));
            eps_vector_pos = data_cell + eps_scale*b;
            denoised_cell_pos  = vect_KPCA_shrinkage( eps_vector_pos, h, c, flag_shrinkage); %Outpuf of KPCA denoising of positive increment of signal f_lambda( y+epsilon*  ) Eq.14 of [1]
            denoised_cell_pos = denoised_cell_pos(:);
            
            b=b(:);
            divSure_random = (1/eps_scale)*b' *(denoised_cell_pos - denoised_cell);  %Eq.17 of [1]  Limit of expectaction operator that approaches the divergence Eq.14
            First_term1 = (1/P) *(norm(data_cell(:)-denoised_cell)).^2;
            SURE(cc,hh) = First_term1 - sigma_noise^2 + ((2*sigma_noise^2)/P) *divSure_random; % Eq.5 of [1]
    end
end
%SURE
[~,Jindex]=min( SURE(:) + 0*sigma_noise^2); % We need only to minimize terms in * which depends on cc,hh
[indx_cc,indx_hh] = ind2sub([numel(cvector),numel(hvector)],Jindex);
copt = cvector(indx_cc);
hopt = hvector(indx_hh);



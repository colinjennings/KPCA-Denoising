function [ Signal ] = vect_KPCA_shrinkage( Data, h, c, flag_shrinkage)

N = size(Data,2);
onesv = ones(N,1);
onesvT=onesv';
%% Kernel matrix definition K(xi,xj)
X2 = diag(Data' * Data);
Distance_l2 =  real(sqrt ( X2 + X2'  - 2* (Data')*Data));
Distance_l2_dummy = Distance_l2 + diag(ones(N,1)*inf);

disp_class = sqrt ( (1/N)*sum( min(Distance_l2_dummy.^2,[],2) ) );
sigma = c  *  disp_class; 
variance = sigma^2;
h_scale = 2*variance;
K = exp(-(Distance_l2.^2)./(h_scale)); 
H = eye(N,N) - (1/N) * (onesv * (onesvT)); %Centering matrix
Khat = H*K*H; % Centered Kernel matrix
Khat(isinf(Khat))=0;
Khat(isnan(Khat))=0;
%[A,Lambda,~] = svd(Khat); %Singular value decomposition %is symetric so U=V.
%
[A,E] = eig(Khat*Khat'); % surprisingly, it turns out that this is generally faster than MATLAB's svd, svds, or eigs commands
[~,idx] = sort(abs(diag(E)),'descend');
A = A(:,idx);
new=diag(E);
Lambda= sqrt(diag(new(idx)));

%
Normalization = 1./sqrt(Lambda);
Normalization(isinf(Normalization)) = 0;
A = A*Normalization;
%% Given a point x and a dimension of projected space, H,

%Definition of Ktilde Eq.2
%kx = s
%ky = 
sumkx = sum(K,2);
sumky = sumkx';
Ktilde = K - (1/N)*repmat(sumkx,1,N) - (1/N)*repmat(sumky,N,1) + (1/N^2)*onesvT*(K)*onesv;  %  Xeval  X Xj

% Beta definition Eq.1
Beta =  ( A(:,1:h)' * Ktilde');   %h  X Xeval

%% OptShrink

if flag_shrinkage
    Lambda_KPCA = diag(Lambda)/N;
    Singular_value_KPCA = sqrt(Lambda_KPCA*(N));
    rank_estim = h;
    Noise_only_SV = Singular_value_KPCA( (rank_estim+1):end);
    Signal_only_SV = Singular_value_KPCA(1:rank_estim);
    Cov_noise = diag(Noise_only_SV);
    dim_noise = N-(rank_estim);

    matrix_core = zeros(rank_estim,dim_noise,dim_noise);
    z_matrix_core = matrix_core;
    matrix_core_squared = matrix_core;
    z2_matrix_core_squared = matrix_core;
    Dtransform = zeros(rank_estim,1);
    Dtransform_deriv = zeros(rank_estim,1);

    for n=1:rank_estim
        matrix_core(n,:,:) = inv(   (Signal_only_SV(n).^2)*eye(dim_noise) - Cov_noise*Cov_noise'  );
        z_matrix_core(n,:,:) = Signal_only_SV(n) * matrix_core(n,:,:);
        Dtransform(n) = (1/dim_noise) * trace(squeeze(z_matrix_core(n,:,:))) .* (1/dim_noise)*trace(squeeze(z_matrix_core(n,:,:)));
        matrix_core_squared(n,:,:) = squeeze(matrix_core(n,:,:)) * squeeze(matrix_core(n,:,:));
        z2_matrix_core_squared(n,:,:) = -2*(Signal_only_SV(n).^2)*matrix_core_squared(n,:,:);

        Dtransform_deriv(n) = (1/dim_noise) * trace( squeeze(z_matrix_core(n,:,:)) ) * (1/dim_noise)*trace( squeeze(z2_matrix_core_squared(n,:,:)) + squeeze(matrix_core(n,:,:)) );
        Dtransform_deriv(n) = 2*Dtransform_deriv(n);  
    end

    wopt = -2*Dtransform./Dtransform_deriv;
    Attenuation = wopt./Signal_only_SV;
    Beta = Beta(1:h,:).*Attenuation;
end
%%

% Gamma definition Eq.7 
gamma = A(:,1:h) * Beta(1:h,:); %Xj  x Xeval
gammatilde = gamma + (1/N) * (1- sum(gamma,1)); %do the centering  %Xj  x Xeval

% Calculation of Eq.4
M=0;
for hh=1:h
    M = M + A(:,hh) * A(:,hh)';
end
%%till here is the same as the for-loop method
HMatrix = H' * M * H;
VX = HMatrix * ( K' - (1/N)*K*onesv) ;  %Xj x X eval on the right 
Term2 =  ( ((1/N)*K*onesv)' - 2*K)* VX;  %Xi x Xeval
Term1 = diag(K*VX); % Xeval
Term12 = repmat(Term1',N,1) + Term2;
Termi =  (1/N^2)*onesvT*(K)*onesv + diag(K) - (2/N)*sumkx;
Dist2 = Term12 +  Termi;  % Xi   X  Xeval
Dist = 1 - 0.5*Dist2;
Weights = gammatilde.*Dist;  % Xi x Xeval
Xtilde = Data * Weights;
Constant = sum(Weights,1);

Signal = Xtilde./Constant;

end


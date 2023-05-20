# DiffusionMRIDenoisingKPCA
"SNR-enhanced diffusion MRI with structure-preserving low-rank denoising in reproducing kernel Hilbert spaces"

We exploit nonlinear redundancy of the dMRI data by means of kernel principal component analysis (KPCA), a nonlinear generalization of PCA to reproducing kernel Hilbert spaces. By mapping the signal to a high-dimensional space, a higher level of redundant information is exploited, thereby enabling better denoising than linear PCA. We implement KPCA with a Gaussian kernel, with parameters automatically selected from knowledge of the noise statistic.

Use:

function out = KPCA_denoising(dwi,mask,kernel,flag_shrinkage)

dwi: diffusion-weighted image (no b0-image)
mask:
kernel:
flag_shrinkage: 

# DiffusionMRIDenoisingKPCA
"SNR-enhanced diffusion MRI with structure-preserving low-rank denoising in reproducing kernel Hilbert spaces"

We exploit nonlinear redundancy of the dMRI data by means of kernel principal component analysis (KPCA), a nonlinear generalization of PCA to reproducing kernel Hilbert spaces. By mapping the signal to a high-dimensional space, a higher level of redundant information is exploited, thereby enabling better denoising than linear PCA. We implement KPCA with a Gaussian kernel, with parameters automatically selected from knowledge of the noise statistic.

Use:

function out = KPCA_denoising(dwi,mask,kernel,flag_shrinkage)

1) dwi: diffusion-weighted image (no b0-image) (nx x ny x nz x ndir)
2) mask: (nx x ny x nz)
3) kernel: kernel for denoising (odd size)
4) flag_shrinkage:  1--> Apply shrinkage ,   0--> Normal KPCA


If you find this code useful, please cite the following papers,

Ramos-Llordén G, Vegas-Sánchez-Ferrero G, Liao C, Westin CF, Setsompop K, Rathi Y. SNR-enhanced diffusion MRI with structure-preserving low-rank denoising in reproducing kernel Hilbert spaces. Magn Reson Med. 2021 Sep;86(3):1614-1632. doi: 10.1002/mrm.28752.

The algorithm makes use fo the noise estimator proposed in 

Aja-Fernández S, Pieciak T, Vegas-Sánchez-Ferrero G. Spatially variant noise estimation in MRI: a homomorphic approach. Med Image Anal. 2015 Feb;20(1):184-97. doi: 10.1016/j.media.2014.11.005. Epub 2014 Nov 24. PMID: 25499191.




![Imagen 1](https://github.com/gabrll/DiffusionMRIDenoisingKPCA/assets/49204215/15ac0991-c3e0-4e3e-8882-432adbe7fd91)
![Uploading color_FA.gif…]()

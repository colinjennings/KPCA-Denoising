
% Helpful script to run and save denoised file 
% - requires freesurfer matlab scripts

% To run: matlab -nodisplay -nosplash -nodesktop -r "cd('/path/to/KPCA-Denoising'); denoise('/path/to/dwi_scan.nii.gz', '/path/to/b0_mask.nii.gz', '/path/to/output.nii.gz'); exit;"
function denoise(subject_file, mask_file, output_file)
	addpath /data/pnl/soft/pnlpipe3/fs7.1.0/matlab
	nifti_file = MRIread(subject_file);
	data = nifti_file.vol;
	mask = MRIread(mask_file).vol;

	% default size
	kernel = [5 5 5];

    temp_lpf = regexprep(output_file, '\.nii\.gz$', '_lpf.mat');
    % 0 is no default no shrinkage
	denoised = KPCA_denoising(data,mask,kernel,0,temp_lpf);

	nifti_file.vol = denoised;

	MRIwrite(nifti_file, output_file);
end
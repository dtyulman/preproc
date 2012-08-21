function [normfiles, normsource, parameters] = spm2_norm(template, spm_params, spm2_dir)
% BOLD Normalization using SPM2
 
disp('Converting to Analyze format...')
!unset LD_LIBRARY_PATH && mri_convert bold.nii -ot spm tmp_spm-
!unset LD_LIBRARY_PATH && mri_convert refvol.nii -ot spm refvol-

addpath(spm_params);
addpath(spm2_dir);

spm_defaults;

P = spm_get('Files', '' ,'tmp_spm-*.img');
meanf = spm_get('Files', '', 'refvol-001.img');

% Get header information for images
V = spm_vol(P); 
Vm = spm_vol(meanf);

% Normalization
defaults.normalise.write.vox = [2 2 2];
defaults.normalise.write.interp = 7;
defaults.normalise.write.wrap = [0 1 0];
defaults.normalise.write.bb = [[-90 -126 -72];[ 90 90 108]];

% Estimate unwarping parameters
disp('Estimating normalization parameters...')
matname = [spm_str_manip(Vm.fname, 'sd') '_sn.mat'];
VG = template;
params = spm_normalise(VG, Vm, matname, '','', defaults.normalise.estimate);
snMask = spm_write_sn(V,params,defaults.normalise.write,'mask');
spm_write_sn(Vm,params,defaults.normalise.write, snMask);

% Write normalized
disp('Normalizing...')
for ii=1:length(V),
    spm_write_sn(V(ii),params,defaults.normalise.write,snMask);
end

disp('Converting BOLD scan back to NIFTI format...')
!unset LD_LIBRARY_PATH && fslmerge -t wtmp_spm.nii.gz wtmp_spm-???.img
disp('Restoring correct BOLD orientation (flip across x axis)...')
!unset LD_LIBRARY_PATH && fslswapdim wtmp_spm.nii.gz -x y z bold_atl.nii.gz
disp('Restoring correct reference volume orientation (flip across x axis) and converting back to NIFTI...')
!unset LD_LIBRARY_PATH && fslswapdim wrefvol-001.img -x y z wrefvol.nii.gz

quit




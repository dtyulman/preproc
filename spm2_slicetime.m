function [timecorrected_file] = spm2_slicetime(order, refslice, TR, nslices, spm_params, spm2_dir)

addpath(spm_params);
addpath(spm2_dir);

spm_defaults

disp('Converting to Analyze format...')
!unset LD_LIBRARY_PATH && mri_convert bold.nii -ot spm tmp_spm-

P = spm_get('Files','','tmp_spm-*.img');

timing = [TR/nslices, (TR-TR/nslices)/(nslices-1)];

if strcmp(order, 'ascending')
    sliceorder = [1:1:nslices];
else if strcmp(order, 'descending')
        sliceorder = [nslices:-1:1];
    else if strcmp(order, 'interleaved (middle-top)')
            for k = 1:nslices,
                sliceorder(k) = round((nslices-k)/2 + (rem((nslices-k),2) * (nslices - 1)/2)) + 1;
            end
        else if strcmp(order, 'interleaved (bottom -> up)')
                sliceorder = [1:2:nslices 2:2:nslices];
            else if strcmp(order, 'interleaved (top -> down)')
                    sliceorder = [nslices:-2:1, nslices-1:-2:1];
                else
                    sprintf('%s is not a valid entry for slice order. Defaulting to interleaved (bottom -> up).', order)
                    sliceorder = [1:2:nslices 2:2:nslices];
                end
            end
        end
    end
end

disp('Performing slice timing correction using SPM2')
spm_slice_timing_noui(P, sliceorder, refslice, timing);

disp('Converting BOLD scan back to NIFTI format...')
!unset LD_LIBRARY_PATH && fslmerge -t bold_faln atmp_spm-*.img

quit





% add tools
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

% directory containing DMNFC maps
indir = '/scratch/groups/leanew1/xcpd_outP50_36p_bp/bv_rs';
files = dir(fullfile(indir, '*_DMNFC_map.dscalar.nii'));

% initialize accumulator
for i = 1:length(files)
    filepath = fullfile(files(i).folder, files(i).name);
    cifti = read_cifti(filepath);

    if i == 1
        sum_data = cifti.cdata;
    else
        sum_data = sum_data + cifti.cdata;
    end
end

% compute mean by divide by n
mean_data = sum_data / length(files);

% Save group average
group_cifti = cifti;
group_cifti.cdata = mean_data;

outpath = '~/GroupAvg_DMNFC_map.dscalar.nii';
ciftisave(group_cifti, outpath);

niftiout='~/GroupAvg_DMNFC_map.nii.gz';
% also save it as a nifti
command = ['wb_command -cifti-separate ' outpath ' COLUMN -volume-all ' niftiout];
system(command)

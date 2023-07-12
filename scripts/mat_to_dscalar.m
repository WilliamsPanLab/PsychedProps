function mat_to_dscalar(subj)

% add matlab paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));

% directory
OutputFolder = ['/oak/stanford/groups/leanew1/users/apines/data/SingleParcel_1by1/' subj '-concat.dtseries.nii/IndividualParcel_Final_sbj1_comp4_alphaS21_1_alphaL300_vxInfo1_ard0_eta0/'];

% resampled V of interest (group or individ - CHANGE AS FIT TO MATCH NAME)
finalUVFile=[OutputFolder 'final_UV.mat'];
Loading_Mat=load(finalUVFile);
V=Loading_Mat.V;

% extract from struct
V=V{:};

% cifti to replace cdata in
HP=read_cifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients.dscalar.nii');
HP.cdata(1:59412,1:4)=V(1:59412,1:4);
outputfile=[OutputFolder 'SoftParcel.dscalar.nii'];
write_cifti(HP,outputfile)

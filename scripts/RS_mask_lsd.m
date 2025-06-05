function RS_mask_lsd(subj,sesh)
% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% load in boolean scan boundaries tmask to glean RS
Boundmaskfp=['/scratch/users/apines/LSD_ICL/rest_proc/' subj '_' sesh '_sr_bpss_tmask.txt'];
Boundmask=load(Boundmaskfp);

% load in aggregate cifti 
CiftiFp=['/scratch/users/apines/LSD_ICL/rest_proc/' subj '_' sesh '_sr_bpss.dtseries.nii'];;
Cifti_file=read_cifti(CiftiFp);

% ensure length of boolean RS mask matches length of aggregate cifti
if size(Cifti_file.cdata, 2) ~= length(Boundmask)
    error('Error: Dimensions of boolean RS mask do not match the length of the aggregate CIFTI file.');
end

% ensure 660 length
if size(Cifti_file.cdata, 2) ~= 660
    error('Error: Non 660 length scan');
end

% make directory
subjDir=['/scratch/users/apines/LSD_ICL/rest_proc/' subj ];
system(['mkdir ' subjDir]);

% split scans into rs1 music and rs2
Ciftifp_rs1=[subjDir '/' subj '_' sesh '_rs1.dtseries.nii'];
Ciftifp_rs2=[subjDir '/' subj '_' sesh '_rs2.dtseries.nii'];
Ciftifp_mus=[subjDir '/' subj '_' sesh '_mus.dtseries.nii'];

% segment out indy scans
Cifti_rs1=Cifti_file.cdata(:,1:220);
Cifti_rs2=Cifti_file.cdata(:,441:660);
Cifti_mus=Cifti_file.cdata(:,221:440);

% reset cifti metadata
Cifti_file.diminfo{2}.length=220;
% save each out
Cifti_file.cdata=Cifti_rs1;
write_cifti(Cifti_file,Ciftifp_rs1);
Cifti_file.cdata=Cifti_rs2;
write_cifti(Cifti_file,Ciftifp_rs2);
Cifti_file.cdata=Cifti_mus;
write_cifti(Cifti_file,Ciftifp_mus);

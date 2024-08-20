ToolFolder = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder))

ProjectFolder = '/oak/stanford/groups/leanew1/users/apines/data/SingleParcellation'
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf/'

% for surface data
surfL = [SubjectsFolder 'lh.pial'];
surfR = [SubjectsFolder 'rh.pial'];

surfML = [SubjectsFolder 'lh.Mask_SNR.label'];
surfMR = [SubjectsFolder 'rh.Mask_SNR.label'];

[surfStru, surfMask] = getFsSurf(surfL, surfR, surfML, surfMR);

% uncomment to implement without SNR mask.
%surfMask.l=ones(length(surfMask.l),1);
%surfMask.r=ones(length(surfMask.r),1);

gNb = createPrepData('surface', surfStru, 1, surfMask);

% save gNb into file for later use
prepDataName = [ProjectFolder '/CreatePrepData.mat'];
save(prepDataName, 'gNb');

function run_FlowPortait_mice(subj,run)

% add path to flow portrait
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/FLOWPortrait/'));
% thank you Nate Linden!!! https://github.com/natejlinden/FLOWPortrait/tree/master
% seriously, thank you thank you thank you

% base fp for LSD
basefp='/scratch/users/apines/p50_mice/proc/20200228/'

% load in specified scan
if run==1
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_preLSD0p3mgkg_1/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
%% add if/else
% post 1
if run==2
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_0/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 2
if run==3
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_5/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 3
if run==4
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_10/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 4
if run==5
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_15/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 5
if run==6
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_20/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end

% convert to double
data=data;
size(data)
% convert mask
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet=networks.components_matrix';
Mask=networks.mask;
% need to re-boolean from the interpolation in the downsample
Mask=Mask==1;
% need to transpose mask to match python-matlab difference: to be confirmed via visualization
Mask=Mask';

% initializing with 20
integrationLen=20;
flowPortrait(data,integrationLen,'white_background',true,'mask',Mask,'filename',['~/' subj '_' num2str(run) '.png'],'save_ftle',true);


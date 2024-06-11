function Viz_FTLE(subj,sesh)
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts')
% load in ftle (i'm aware that the .png_ftle.mat is dumb, will fix)
ftlefp=['~/' subj '_' num2str(sesh) '.png_ftle.mat'];
ftle=load(ftlefp).ftle;
% load in time series (interpolated)
%basefp='/scratch/users/apines/p50_mice/proc/20200228/'
%ofInterpFp=[basefp '/' subj '_run-' num2str(sesh) '_interp.mat'];
%data=load(ofInterpFp);
% get length of time series
%lenOpFl=size(data.InterpData);
%lenOpFl=lenOpFl(3)
% for some timepoints
for t=300:450
    filename=['/scratch/users/apines/' subj '_' num2str(sesh) '_t' num2str(t) '_']; 
    % save out png (imagesc of signal in inferno)
    figure;
    imagesc(ftle.f(:,:,t));  % Adjust to your specific data structure
    colormap(inferno);
    caxis([0.0; .015]);
    print([filename 'f.png'],'-dpng','-r600')
    figure;
    imagesc(ftle.b(:,:,t));  % Adjust to your specific data structure
    colormap(inferno);
    caxis([0.0; .015]);
    print([filename 'b.png'],'-dpng','-r600')
end


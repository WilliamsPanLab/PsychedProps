% load in inferno colormap through addpath
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts');
subjects={'m2000','m7502','m7507','m7520','m7589','m7594'};

sessions = {1,2,3,4,5,6};
% for each session
for j = 1:length(sessions)
	sesh = sessions{j};
	% initialize mean ftle
	meanftle_f=zeros(67,70,6);
	meanftle_b=zeros(67,70,6);
	% obs counter
	obs=0;
	% for each subject
	for i=1:length(subjects)
	        subj = subjects{i}
                % set filepath that things would be saved to
                outFP=['/scratch/users/apines/data/mouse'];
                % if file exists, then continue
                filename=[outFP '/' subj '_' num2str(sesh) '_Prop_Feats_gro.csv'];
                if isfile(filename)
			% load in ftle (i'm aware that the .png_ftle.mat is dumb, will fix)
			ftlefp=['~/' subj '_' num2str(sesh) '.png_ftle.mat'];
			ftle=load(ftlefp).ftle;
			% save mean f and b from each mice as continuous map
			meanftle_f(:,:,i)=mean(ftle.f,3);
			meanftle_b(:,:,i)=mean(ftle.b,3);
			obs=obs+1;
		else
		end
	end
	% average by number of observations
	meanftle_f=mean(meanftle_f,3);
	meanftle_b=mean(meanftle_b,3);
	filename1=['/home/users/apines/meanFTLE_f_' num2str(sesh) '.png'];
	filename2=['/home/users/apines/meanFTLE_b_' num2str(sesh) '.png'];
	figure;
	% something about this keeps crashing matlab... lets try rounding the numbers
	%meanftle_f=round(meanftle_f);
	% it's not round...
	imagesc(meanftle_f);  % Adjust to your specific data structure
	colormap(inferno);
	colorbar
	caxis([0.0; .003]);	
	print(filename1,'-dpng','-r600')
	figure;
	imagesc(meanftle_b);  % Adjust to your specific data structure
	colormap(inferno);
	colorbar
	caxis([0.0; .003]);	
	print(filename2,'-dpng','-r600')
end

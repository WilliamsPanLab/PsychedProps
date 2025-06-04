function Extract_DMNMag_mice_Spun(subj,sesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take optical flow results, get a bottom-up and top-down resultant vector in x,y coords for each pixel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% Load in flatmouse opflow calc
% note this is for LSD only! Adapt recording date for other drugs if needed
childfp='/scratch/users/apines/p50_mice/proc/20200228/'
datafp=[childfp subj '_vf_out_' num2str(sesh) '.mat']
% adding in if it exists: ends all the way at the end of the script
if exist(datafp)
data=load(datafp)
% vector fields
vf=data.outputs.velocityFields;
% and actual signal
signalGrid=data.outputs.filteredSignal;
% load in mask to repopulate consensus onto 133x139 grid
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Mask=networks.mask;
% need to re-boolean from the interpolation in the downsample
Mask=Mask==1;
% need to transpose mask to match python-matlab difference: to be confirmed via visualization
Mask=Mask';
% get size of mask and indices to include
MaskSize=sum(sum(Mask));
MaskInds=find(Mask);
% save out mask as csv

% extract size of time series
NumTRs=size(vf);
NumTRs=NumTRs(3);
lenOpFl=NumTRs;

% and an "over time" version of the mask to apply to the time series
fullMask = repmat(Mask, [1, 1, lenOpFl]);

% initialize out dataframes
Propvec=[];
stringVec={};

% and azez and els for opflow vectors
OpF_x=zeros(size(vf,1),size(vf,2),lenOpFl);
OpF_y=zeros(size(vf,1),size(vf,2),lenOpFl);
for tp=1:lenOpFl
	OpF_x(:,:,tp)=real(vf(:,:,tp));
	OpF_y(:,:,tp)=imag(vf(:,:,tp));
end

% mask out angles
OpF_x(~fullMask)=0;
OpF_y(~fullMask)=0;

%%%%%%%%%%%%%%%%%%%% Code to plot stuff currently loaded in to make sure it is lined up
% print out OpFl Vectors from tp 1, mask, masked vectors to confirm alignment
%fig = figure;
% plots signal from one timepoint
%tp=1;
%imagesc(signalGrid(:, :, tp));  % Adjust if 'signalGrid' has different dimensions
%hold on;  % Hold on to plot quiver on top
% pull this timepoint of vectors
%currentOpF_x=OpF_x(:,:,1);
%currentOpF_y=OpF_y(:,:,1);
% Overlay the vector field using quiver
%[x, y] = meshgrid(1:size(OpF_x, 2), 1:size(OpF_x, 1));
%quiver(x(Mask), y(Mask), currentOpF_x(Mask), currentOpF_y(Mask), 'y');
% Add labels and titles for clarity
%title('Optical Flow Vectors on Masked Signal Grid');
%xlabel('X Coordinate');
%ylabel('Y Coordinate');
%axis tight;  % Fit the axes tightly to the data
%hold off;  % Release the hold
% Save the figure to a file
%saveas(fig, '~/flow_vectors.png');
% Close the figure programmatically
%close(fig);
% Beautiful. Commenting out for now.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize output vector
SpunMags=zeros(1,2000);

% load in spins
SpunDMNS=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_permutations_10k.mat');
for k=1:2000
	Dnet=SpunDMNS(:,:,k);
	% create pixel-wise network mask
	DMN_bool=Dnet;
	DMN_bool(DMN_bool>.6)=1;
	DMN_bool(DMN_bool<.6)=0;
	% adding mask into this step 7/1/24
	DMN_bool(Mask==0)=0;
	DMN_bool=logical(DMN_bool);
	
	% initialize matrix for each pixel to saveout to scratch
	faceMatrix=zeros(sum(sum(DMN_bool)));

	%%%%%%%%%%%%%%%%%%%%%%%%%%% temporary visualization code to triple check stuff
	%%%%% DMN gradient + vectors
	%figure;
	%imagesc(net);
	%colormap('jet');
	%hold on;
	% Create a grid for the quiver plot
	%[x, y] = meshgrid(1:size(net, 2), 1:size(net, 1));
	% Hold on to the current image and overlay the quiver plot
	%hold on;
	%quiver(x, y, nGx, nGy);
	%hold off;
	%print('~/DMUnder.png','-dpng','-r600');
	%%%%% Ca2+ signal+vectors, DMN gradient vectors
       	%for fr=100:200
	%	fig=figure;
	%	imagesc(signalGrid(:,:,fr));
	%	colormap('jet');
	%	caxis([min(min(signalGrid(:,:,30)))-0.002,max(max(signalGrid(:,:,30)))+0.002]);
	%	hold on;
		% Create a grid for the quiver plot
        %	[x, y] = meshgrid(1:size(net, 2), 1:size(net, 1));
        	% Hold on to the current image and overlay the quiver plot
        %	hold on;
        %	quiver(x, y, nGx, nGy);
        %	hold on;
		% quiver for of vectors
	%	quiver(x,y,OpF_x(:,:,fr),OpF_y(:,:,fr),'w')
	%	hold off;
	%	print(fig, ['~/DMNGrad_SignalUnder_' num2str(fr)], '-dpng', '-r600')
	%end
	%%%%%% Boolean DMN mask, Ca2+ signal and DMN gradient vectors
	%fig=figure;
	%DMN_bool=Dnet;
        %DMN_bool(DMN_bool>.6)=1;
        %DMN_bool(DMN_bool<.6)=0;
        %imagesc(DMN_bool);
        %colormap('jet');
        %hold on;
        % Create a grid for the quiver plot
        %[x, y] = meshgrid(1:size(net, 2), 1:size(net, 1));
        % Hold on to the current image and overlay the quiver plot
        %hold on;
        %quiver(x, y, nGx, nGy);
        %hold on;
        % quiver for of vectors
        %quiver(x,y,OpF_x(:,:,1),OpF_y(:,:,1),'w')
        %hold off;
        %print(fig, '~/DMNGrad_BooleanUnder', '-dpng', '-r600')	
	% reconversion to logical
	%DMN_bool=logical(DMN_bool);
	% looks good, commenting out for now
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % initialize angular distance vector for each network (l and r) above
        Mags=zeros(sum(sum(DMN_bool)),lenOpFl);
	% get angular distance for each face for each timepoint
        for F=1:sum(sum(DMN_bool));
                % loop over each tp
                for fr=1:lenOpFl
			% pull out optical flow vectors for this pixel (denoted by F, because it represents faces in human workflow)
			curOpF_x=OpF_x(:,:,fr);
			curOpF_y=OpF_y(:,:,fr);
			curOpF_x_DMN=curOpF_x(DMN_bool);
			curOpF_y_DMN=curOpF_y(DMN_bool);
                        % get optical flow vector
                        OpFlVec=[curOpF_x_DMN(F) curOpF_y_DMN(F)];
			
			% calculate magnitude!
                        Mags(F,fr)=sqrt(sum(OpFlVec .^2));
			%%%% quadruple-check figure: plot one DMN grad angle, one opfl angle, set title to a (angular distance)
			%fig=figure('Visible','off');
			% Plot the reference vector (assuming nVec is already defined in your workspace)
			%quiver(0, 0, nVec(1), nVec(2), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'k', 'DisplayName', 'Reference Vector');
			%hold on;
			% optical flow vectors
			%quiver(0, 0, OpFlVec(1), OpFlVec(2), 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Color', 'r', 'DisplayName', 'Optical Flow Vector');
			%legend show;
			%title(sprintf('Angular Distance: %.2f degrees', a));
			%axis equal;
			%xlim([-max(abs([nVec, OpFlVec]))*1.5, max(abs([nVec, OpFlVec]))*1.5]);
			%ylim([-max(abs([nVec, OpFlVec]))*1.5, max(abs([nVec, OpFlVec]))*1.5]);
			% Plot the optical flow vector
			%saveas(fig, ['~/vector_plot' num2str(fr) '.png']);
			%%%% Seems to work, commenting out

                % end tp loop
                end
	% end each face loop
        end
	% calc outFP
	outFP=['/scratch/users/apines/data/mouse/'];
	% time series population
	OutTs=Mags;
	writematrix(OutTs,[outFP subj '_' num2str(sesh) '_Mag_TS_dmn.csv'])
	% average angular distances across hemispheres
        avgD=mean(mean(Mags));
        SpunMags(k)=avgD;
end
% add label
stringVec = compose("Spin%d", 1:2000);
%save out as csv
T=table(stringVec','RowNames',stringVec);
% write out
writetable(T,[outFP subj '_' num2str(sesh) '_Prop_Feats_Mag_Spun.csv'],'WriteRowNames',true)
else
	disp('file not found')
end

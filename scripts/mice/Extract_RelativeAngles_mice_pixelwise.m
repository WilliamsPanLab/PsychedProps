function Extract_RelativeAngles_mice(subj,sesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take optical flow results, get a bottom-up and top-down resultant vector in x,y coords for each pixel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% Load in flatmouse opflow calc
% note this is for LSD only! Adapt recording date for ketamine if needed
childfp='/scratch/users/apines/p50_mice/proc/20200228/'
datafp=[childfp subj '_vf_out_' num2str(sesh) '.mat']
% adding in if it exists: ends all the way at the end of the script
if exist(datafp)
data=load(datafp)
% load in mask
%fn='/scratch/users/apines/p50_mice/proc/20200228/thy1gc6s_0p3mgkg_m2000_preLSD0p3mgkg_1/masked_dff_DS_BP_Smoothed.h5'
%formask=h5read(fn, '/mask');
%mask = cellfun(@(x) strcmpi(x, 'TRUE'), formask);
% might need a surf
%surf=mask;
% vector fields
vf=data.outputs.velocityFields;
% and actual signal
signalGrid=data.outputs.filteredSignal;
% get incenters of pixels
% TR_L = TriRep(F_L,V_L);

% USED TO get INBETWEEN SPOTS IN VIS SCRIPT
%x=linspace(1, size(signalGrid,2), size(vf,2))
%y=linspace(1, size(signalGrid,1), size(vf,1))

%P_L = TR_L.incenters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make medial wall vector
%g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
% load in mask to repopulate consensus onto 133x139 grid
% smoothed dmn version
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet=networks.components_matrix';
Mask=networks.mask;
% need to re-boolean from the interpolation in the downsample
Mask=Mask==1;
% need to transpose mask to match python-matlab difference: to be confirmed via visualization
Mask=Mask';
% get size of mask and indices to include
MaskSize=sum(sum(Mask));
MaskInds=find(Mask);
% extract size of time series
NumTRs=size(vf);
NumTRs=NumTRs(3);
lenOpFl=NumTRs;

% and an "over time" version of the mask to apply to the time series
% fullMask = repmat(Mask, [1, 1, lenOpFl]);
% NOT FOR PIXELWISE

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
% translate xyz spherical coordinates to az/el/r
%[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
% convert from radians to degrees
%azd_L=rad2deg(az_L);
%eld_L=rad2deg(el_L);

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

% add lukas lang functions
%addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))

% AP - 5/13/24: just testing full-brain at this point can break down by networks later if needed
%for k=1:4
for k=1
	nets=Dnet;
	% create pixel-wise network mask
	DMN_bool=Dnet;
	DMN_bool(DMN_bool>.6)=1;
	DMN_bool(DMN_bool<.6)=0;
	DMN_bool=logical(DMN_bool);
	% initialize matrix for each pixel to saveout to scratch
	faceMatrix=zeros(sum(sum(DMN_bool)));
        % network of interest
        net=Dnet;
	% calculate gradient of DMN
        [nGx, nGy] =imgradientxy(net);

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
       	%fig=figure;
	%imagesc(signalGrid(:,:,1));
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
	%print(fig, '~/DMNGrad_SignalUnder', '-dpng', '-r600')
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
        NangDs=zeros(size(net,1),size(net,2),lenOpFl);
	% get angular distance for each face for each timepoint
	% for each x and y (will be redundant because of symmetry of x y iterations)
        for X=1:size(net,1);
		for Y=1:size(net,2);
                	% get vector for each face (network vector)
                	nVec=[nGx(X,Y) nGy(X,Y)];
                	% loop over each tp
                	for fr=1:lenOpFl
				% pull out optical flow vectors for this pixel (denoted by F, because it represents faces in human workflow)
				curOpF_x=OpF_x(X,Y,fr);
				curOpF_y=OpF_y(X,Y,fr);
                	        % get optical flow vector
                	        OpFlVec=[curOpF_x curOpF_y];
				% if in valid pixel
				if sum(OpFlVec)~=0;
					% get angular distance at that timepoint (degrees)
                	        	a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
                        
					% populate vector
                        		NangDs(X,Y,fr)=a;
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
				else
				end

                	% end tp loop
                	end
		% end each y
		end
	% end each x
	end
        % average values over time and plop into facematrix for this participant
        faceMatrix=NangDs;
	% calc outFP
	outFP=['/scratch/users/apines/data/mouse/'];
	% save out faceMatrix with subject ID as csv to /scratch/users/apines/gp/PropFeatsTemp
	writematrix(faceMatrix,['/scratch/users/apines/gp/PropFeats/' subj '_' num2str(sesh) '_faceMatrix_gro_pixelwise.csv'])
end
else
	disp('file not found')
end

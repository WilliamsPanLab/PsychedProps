function Extract_RelativeAngles_CaudR(subj, sesh, task)
% Measure % bottom-up optical flow direction relative to DMN gradient in 3D
ToolFolder = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% Set filepaths
atlas_fp = '/oak/stanford/groups/leanew1/users/apines/maps/Tian_Subcortex_S1_3T_32k.dscalar.nii';
dmn_fp   = '~/GroupAvg_DMNFC_map.nii.gz';
childfp  = ['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj '/' sesh '/func'];
flow_fp  = [childfp '/' subj '_CaudR_vf_out_' sesh '_' task '.mat'];

% load it in
atlas = read_cifti(atlas_fp);
dmn   = niftiread(dmn_fp);
load(flow_fp);
% segment out this subcort struct
atlasInds = atlas.cdata;
subcortVoxels = [];
for s = 3:length(atlas.diminfo{1}.models)
    subcortVoxels = [subcortVoxels; atlas.diminfo{1}.models{s}.voxlist];
end
caud_R = find(atlasInds == 8) - 59412;
caud_R_coords = subcortVoxels(caud_R, :) + 1;
% set caud R mask
caud_R_mask = false(size(dmn,1), size(dmn,2), size(dmn,3));
% now set voxels in caud R as true
for i = 1:size(caud_R_coords,1)
    x = caud_R_coords(i,1);
    y = caud_R_coords(i,2);
    z = caud_R_coords(i,3);
    caud_R_mask(x, y, z) = true;
end
% apply to dmn
dmn(~caud_R_mask) = 0;
% extract boxes that these voxels live in (we're going to restrict box by 1 to live less on the edge)
% caud R
caud_R_min = min(caud_R_coords,[],1)+1;
caud_R_max = max(caud_R_coords,[],1)-1;
% same adjustment for caudate as in optical flow script
caud_R_min(2)=caud_R_min(2)+1;
% now zoom in on the bounding box
dmn = dmn(caud_R_min(1):caud_R_max(1), ...
                caud_R_min(2):caud_R_max(2), ...
                caud_R_min(3):caud_R_max(3), :);

% initialize counts for EZ tracking across x y and z planes
BUP_count = 0;
total_count = 0;

% get dmn slices to iterate over
x_slices_to_check = find(squeeze(any(any(dmn > 0, 2), 3)));  % x-dimension slices with any DMN loadings

% for each X slice
for x = x_slices_to_check'
	% get this dmn slice
	dmn_slice = squeeze(dmn(x,:,:));		
	% get gradient of it
	[Gy, Gz] = imgradientxy(dmn_slice);
	% threshold DMN map at .1, has to be more permissive than cortical thresholding
	Gy(dmn_slice < 0.1) = 0;
	Gz(dmn_slice < 0.1) = 0;
	% get angular distance for each timepoint
	for t = 1:size(outputs.velocityFields_x,4)
		% now get valued x coords in this slice
		column_to_check=find(any(abs(Gy)>0));
		OF_y = squeeze(real(outputs.velocityFields_x(x,:,:,t)));
		OF_z = squeeze(imag(outputs.velocityFields_x(x,:,:,t)));
		% for each voxel column
		for column = column_to_check;
			% find row coords to check
			row_to_check=find(abs(Gy(:,column))>0);
			% for each voxel row
			for row = row_to_check'
				DMNy=Gy(row,column);
				DMNz=Gz(row,column);
				OF_yval=OF_y(row,column);
				OF_zval=OF_z(row,column);
				% combine
				nVec=[DMNy DMNz];
				OpFlVec=[OF_yval OF_zval];
				% calculate angular distance, same as it ever was
			 	a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
				% add to bup if it's there
				BUP_count = BUP_count + (a < 90);
				total_count = total_count+1;
			end
		end
		% optional visualization module
		%%if ismember(t, [5, 15, 30]) && x == 4  % adjust timepoints as desired
		%%t
	    	%%fig = figure('Visible', 'off'); 
	    	%%tiledlayout(2,3);
	    	% first tile: full slice with vectors
	    	%%nexttile;
	    	%%imagesc(dmn_slice); axis image; hold on;
		%%colorbar;
		% keep same clims
		%%clims = caxis;
	    	%%[Y, Z] = meshgrid(1:size(dmn_slice,2), 1:size(dmn_slice,1));
	    	%quiver(Y, Z, Gy, Gz, 'c'); % DMN gradient in cyan
	    	%%quiver(Y, Z, OF_y, OF_z, 'g'); % OpFlow in green
	    	%%title(sprintf('DMN slice x=%d t=%d', x, t));
	    	%%axis off;
	    	% pick 5 example voxels
	    	%%[r_all, c_all] = find(dmn_slice > 0.1);
	    	%%samp_idx = randsample(length(r_all), min(5, length(r_all)));
	    	%%for k = 1:length(samp_idx)
	        %%	rr = r_all(samp_idx(k));
	        %%	cc = c_all(samp_idx(k));
	        %%	DMNy = Gy(rr,cc); DMNz = Gz(rr,cc);
	        %%	OF_yval = OF_y(rr,cc); OF_zval = OF_z(rr,cc);
	        %%	nVec = [DMNy DMNz]; OpFlVec = [OF_yval OF_zval];
	        %%	if norm(nVec)==0 || norm(OpFlVec)==0
	        %%   		ang = NaN;
	        %%	else
	        %%    		ang = acosd(min(1,max(-1, dot(nVec, OpFlVec) / norm(nVec) / norm(OpFlVec))));
	        %%	end
	        	% plot
	        %%	nexttile;
			% define patch and vector spot
		%%	rmin = max(rr-1, 1); cmin = max(cc-1, 1);
		%%	patch = dmn_slice(rmin:min(rr+1,end), cmin:min(cc+1,end));
		%%	plot_row = rr - rmin + 1;
		%%	plot_col = cc - cmin + 1;
		%%	imagesc(patch);
		%%	caxis(clims);
		%%	axis image; hold on;
    		    	% normalize for plotting
		%%	DMN_norm = sqrt(DMNy^2 + DMNz^2 + eps);
		%%	quiver(plot_col, plot_row, DMNy / DMN_norm, DMNz / DMN_norm, 'c', 'LineWidth', 2);
    		%%    	OF_norm = sqrt(OF_yval^2 + OF_zval^2 + eps);
		%%	quiver(plot_col, plot_row, OF_yval / OF_norm, OF_zval / OF_norm, 'g', 'LineWidth', 2);
		%%	title(sprintf('Vox (%d,%d), a=%.1f°', rr, cc, ang));
        	%%	axis off;
    		%%end
    		% save
    		%%outdir = fullfile('/scratch/users/apines/DMN_OpFlow_diag');
    		%%if ~exist(outdir, 'dir'), mkdir(outdir); end
    		%%	outfile = fullfile(outdir, sprintf('x%d_t%d_diag.png', x, t));
    		%%	saveas(fig, outfile);
    		%%	close(fig);
		%%end
	end
end

% get dmn slices to iterate over
y_slices_to_check = find(squeeze(any(any(dmn > 0, 1), 3)));

% for each Y slice
for y = y_slices_to_check
        % get this dmn slice
        dmn_slice = squeeze(dmn(:,y,:));
        % get gradient of it
        [Gx, Gz] = imgradientxy(dmn_slice);
        % threshold DMN map at .1
        Gx(dmn_slice < 0.1) = 0;
        Gz(dmn_slice < 0.1) = 0;
        % get angular distance for each timepoint
        for t = 1:size(outputs.velocityFields_y,4)
                % now get valued x coords in this slice
                column_to_check=find(any(abs(Gx)>0));
                OF_x = squeeze(real(outputs.velocityFields_y(:,y,:,t)));
                OF_z = squeeze(imag(outputs.velocityFields_y(:,y,:,t)));
                % for each voxel column
                for column = column_to_check;
                        % find rows coords to check
                        row_to_check=find(abs(Gx(:,column))>0);
                        % for each voxel row
                        for row = row_to_check'
                                DMNx=Gx(row,column);
                                DMNz=Gz(row,column);
                                OF_xval=OF_x(row,column);
                                OF_zval=OF_z(row,column);
                                % combine
                                nVec=[DMNx DMNz];
                                OpFlVec=[OF_xval OF_zval];
                                % calculate angular distance, same as it ever was
                                a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
                                % add to bup if it's there
                                BUP_count = BUP_count + (a < 90);
                                total_count = total_count+1;
                        end
                end
	end
end

% for each z slice
% get dmn slices to iterate over
z_slices_to_check = find(squeeze(any(any(dmn > 0, 1), 2)));

% for each z slice
for z = z_slices_to_check'
        % get this dmn slice
        dmn_slice = squeeze(dmn(:,:,z));
        % get gradient of it
        [Gx, Gy] = imgradientxy(dmn_slice);
        % threshold DMN map at .1
        Gx(dmn_slice < 0.1) = 0;
        Gy(dmn_slice < 0.1) = 0;
        % get angular distance for each timepoint
        for t = 1:size(outputs.velocityFields_z,4)
                % now get valued x coords in this slice
                column_to_check=find(any(abs(Gx)>0));
                OF_x = squeeze(real(outputs.velocityFields_z(:,:,z,t)));
                OF_y = squeeze(imag(outputs.velocityFields_z(:,:,z,t)));
                % for each voxel column
                for column = column_to_check;
                        % find row coords to check
                        row_to_check=find(abs(Gx(:,column))>0);
                        % for each voxel row
                        for row = row_to_check'
                                DMNx=Gx(row,column);
                                DMNy=Gy(row,column);
                                OF_xval=OF_x(row,column);
                                OF_yval=OF_y(row,column);
                                % combine
                                nVec=[DMNx DMNy];
                                OpFlVec=[OF_xval OF_yval];
                                % calculate angular distance, same as it ever was
                                a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
                                % add to bup if it's there
                                BUP_count = BUP_count + (a < 90);
                                total_count = total_count+1;
                        end
                end
        end
end

% compute % bottom-up
pct_bottom_up = BUP_count / total_count;
T=table(pct_bottom_up,'RowNames',string('AngD_DMN'));
outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh]
writetable(T,[outFP '/' subj '_' sesh '_' task '_CaudR_Prop_Feats_gro.csv'],'WriteRowNames',true)
end


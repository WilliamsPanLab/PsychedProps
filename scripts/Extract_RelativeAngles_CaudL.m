function Extract_RelativeAngles_CaudL(subj, sesh, task)
% Measure % bottom-up optical flow direction relative to DMN gradient in 3D
ToolFolder = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% Set filepaths
atlas_fp = '/oak/stanford/groups/leanew1/users/apines/maps/Tian_Subcortex_S1_3T_32k.dscalar.nii';
dmn_fp   = '~/GroupAvg_DMNFC_map.nii.gz';
childfp  = ['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj '/' sesh '/func'];
flow_fp  = [childfp '/' subj '_CaudL_vf_out_' sesh '_' task '.mat'];

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
caud_L = find(atlasInds == 16) - 59412;
caud_L_coords = subcortVoxels(caud_L, :) + 1;
% set caud L mask
caud_L_mask = false(size(dmn,1), size(dmn,2), size(dmn,3));
% now set voxels in caud L as true
for i = 1:size(caud_L_coords,1)
    x = caud_L_coords(i,1);
    y = caud_L_coords(i,2);
    z = caud_L_coords(i,3);
    caud_L_mask(x, y, z) = true;
end
% apply to dmn
dmn(~caud_L_mask) = 0;
% extract boxes that these voxels live in (we're going to restrict box by 1 to live less on the edge)
% caud L
caud_L_min = min(caud_L_coords,[],1)+1;
caud_L_max = max(caud_L_coords,[],1)-1;
% same adjustment for caudate as in optical flow script
caud_L_min(2)=caud_L_min(2)+1;
% now zoom in on the bounding box
dmn = dmn(caud_L_min(1):caud_L_max(1), ...
                caud_L_min(2):caud_L_max(2), ...
                caud_L_min(3):caud_L_max(3), :);

% initialize counts for EZ tracking across x y and z planes
BUP_count = 0;
total_count = 0;

% get dmn slices to iterate over
x_slices_to_check = find(squeeze(any(any(dmn > 0, 2), 3)));  % x-dimension slices with any DMN loadings

% for each X slice
for x = x_slices_to_check'
	% get this dmn slice
	dmn_slice = squeeze(dmn(x,:,:));		
	% get gradient of it: note Gx and Gy are RELATIVE, not absolute
	[Gx, Gy] = imgradientxy(dmn_slice);
	% threshold DMN map at .1, has to be more permissive than cortical thresholding, still only returns 58 voxels
	Gx(dmn_slice < 0.1) = 0;
	Gy(dmn_slice < 0.1) = 0;
	% get angular distance for each timepoint
	for t = 1:size(outputs.velocityFields_x,4)
		% now get valued x coords in this slice
		column_to_check=find(any(abs(Gx)>0));
		vx = squeeze(real(outputs.velocityFields_x(x,:,:,t)));
		vy = squeeze(imag(outputs.velocityFields_x(x,:,:,t)));
		% for each x voxel column (relative x!)
		for column = column_to_check;
			% find y coords to check
			row_to_check=find(abs(Gx(:,column))>0);
			% for each y voxel row
			for row = row_to_check'
				DMNx=Gx(row,column);
				DMNy=Gy(row,column);
				OF_x=vx(row,column);
				OF_y=vy(row,column);
				% combine
				nVec=[DMNx DMNy];
				OpFlVec=[OF_x OF_y];
				% calculate angular distance, same as it ever was
			 	a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
				% add to bup if it's there
				BUP_count = BUP_count + (a < 90);
				total_count = total_count+1;
			end
		end
	end
end

dmn_slice = squeeze(dmn_vol(x,:,:));
    if nnz(dmn_slice) == 0, continue; end
    [Gx, Gy] = gradient(dmn_slice);

    for t = 1:size(outputs.velocityFields_x,4)
        vx = squeeze(real(outputs.velocityFields_y(x,:,:,t)));  % note Y
        vy = squeeze(imag(outputs.velocityFields_y(x,:,:,t)));  % note Y

        magG = sqrt(Gx.^2 + Gy.^2);
        magV = sqrt(vx.^2 + vy.^2);
        dotProd = Gx.*vx + Gy.*vy;
        cosTheta = dotProd ./ (magG .* magV + eps);
        angles = acosd(max(-1, min(1, cosTheta)));

        mask = dmn_slice > 0 & ~isnan(angles);
        good_counts = good_counts + nnz(angles(mask) < 90);
        total_counts = total_counts + nnz(mask);
    end
end 







% Loop through dimensions: x, y, z
for dim = ["x", "y", "z"]
    vf = VF.(['velocityFields_' dim]);

    dims = size(vf);
    slice_dim = find(dim == "xyz");

    % Iterate over slices in the given orientation
    for s = 1:dims(slice_dim)
        % Extract DMN slice
        switch dim
            case "x"
                dmn_slice = squeeze(dmn_vol(s,:,:));
                vf_slice = squeeze(vf(s,:,:,:));
            case "y"
                dmn_slice = squeeze(dmn_vol(:,s,:));
                vf_slice = squeeze(vf(:,s,:,:));
            case "z"
                dmn_slice = squeeze(dmn_vol(:,:,s));
                vf_slice = squeeze(vf(:,:,s,:));
        end

        % Threshold DMN for gradient masking
        dmn_thresh = dmn_slice;
        dmn_thresh(dmn_thresh < 0.1) = 0;

        % Compute gradient of DMN
        [Gx, Gy] = imgradientxy(dmn_thresh);

        for t = 1:size(vf_slice, 3)
            Vx = real(vf_slice(:,:,t));
            Vy = imag(vf_slice(:,:,t));

            for i = 1:size(Gx,1)
                for j = 1:size(Gx,2)
                    ref = [Gx(i,j), Gy(i,j)];
                    vec = [Vx(i,j), Vy(i,j)];

                    % Only evaluate if gradient is nonzero (i.e., meaningful direction)
                    if norm(ref) > 0.1 && norm(vec) > 0
                        total_valid = total_valid + 1;
                        angle = acosd(max(min(dot(ref, vec)/(norm(ref)*norm(vec)), 1), -1));
                        if angle < 90
                            total_bottom_up = total_bottom_up + 1;
                        end
                    end
                end
            end
        end
    end
end

% Compute % bottom-up
pct_bottom_up = total_bottom_up / total_valid;

% Save result
T = table(pct_bottom_up, 'VariableNames', {'PctBottomUp'});
writetable(T, outfp_csv);
fprintf('Saved %%BottomUp: %.2f%% to %s\n', pct_bottom_up*100, outfp_csv);
end


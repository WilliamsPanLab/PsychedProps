function [gfl, gfl_ind] = constructW_flat(neiR)
% generate sparse graph for 2D image
% indexMat: a 3D matrix containing the regions of intesest (ROI), the value of non-zero voxel refer to cifti index
% neiR: the radius of spatial neighborhood


% AP - mask load in
% AP load in data: pre LSD
fn='/scratch/users/apines/p50_mice/proc/20200228/thy1gc6s_0p3mgkg_m2000_preLSD0p3mgkg_1/masked_dff_DS_BP_Smoothed.h5'
formask=h5read(fn, '/mask');
mask = cellfun(@(x) strcmpi(x, 'TRUE'), formask);
maskMat = int32(mask);
% we are testing this by visualizing:
imwrite(mask, '~/testMaskMice.png');
% set neighborhood radius
% 5-8-24 now taken as input argument

% get size of matrix
flSz = size(mask);
if length(flSz)~=2
    error('  Only support 2D image for now');
end

% if length of the flat isnt = length of neighbor input (which should be 1)
if length(flSz)~=length(neiR)
    neiRvec = repmat(neiR,1,length(flSz));
else
    neiRvec = neiR;
end

% number of elements - nonzero
eleNum = sum(mask(:)~=0);
% initialize cell structure the size of nonzero elements
gfl = cell(eleNum,1);
gfl_ind = zeros(eleNum,1);
%assign a label (1 to length(non-zero voxel)) for each non-zero voxel in the mask
labMaskMat = maskMat;
nonZero = maskMat~=0;
eleNum = sum(nonZero(:));
labMaskMat(maskMat~=0) = 1:eleNum;


% initialize a counter 
gi = 0;
% for each y element of flat
for yi=1:flSz(2)
   % for each x element of flat
   for xi=1:flSz(1)
	% label of this particular center to check out the neighborhood of
        cenLab = labMaskMat(xi,yi);
	% if cenLab is a real element (nonmask)
        if cenLab~=0
	    % update iterator
	    gi = gi + 1;
       	    % return global label of cenLab in order of masked pixels
	    gfl_ind(gi) = cenLab;
	    % find nstart: looks like the [1,1,] is just there to catch the edge case of the starting iteraitons (max of negative numbers)
            nstart = max([xi,yi]-neiRvec,[1,1]);
	    % same for end: ending two digits (for x and y)
            nend = min([xi,yi]+neiRvec,flSz);
	    for yni=nstart(2):nend(2)
		for xni=nstart(1):nend(1)
			neiLab = labMaskMat(xni,yni);
				if neiLab~=0 && neiLab~=cenLab
					gfl{gi}(end+1) = neiLab;
					% set this guy in the neighbor's vector too. Indent appears off-standard in original code
					gfl{neiLab}(end+1) = cenLab;
				end
			end
		end
        end
    end
end

% save out a few example neighborhoods
targetPixelLabel = gfl_ind(500); % This should be the label of the 500th 'valid' pixel

% Create an RGB image where each pixel starts off as black [0 0 0]
imageOutput = zeros([size(mask), 3]);

% Make in-mask pixels white
for x = 1:size(mask, 1)
    for y = 1:size(mask, 2)
        if mask(x, y)  % Check if the pixel is in the mask
            imageOutput(x, y, :) = [1 1 1]; % White
        end
    end
end

% Color the target pixel red
[targetX, targetY] = find(labMaskMat == targetPixelLabel); % Find the coordinates of the target pixel
imageOutput(targetX, targetY, :) = [1 0 0]; % Red

% Color the neighbors blue
neighbors = gfl{500}; % Get the list of neighbors from the 500th row
for i = 1:length(neighbors)
    [nx, ny] = find(labMaskMat == neighbors(i)); % Find the coordinates of each neighbor
    imageOutput(nx, ny, :) = [0 0 1]; % Blue
end

% Save the image
outputFilename = '~/testMaskMice_Neighbors500.png'; % Specify the path and filename for output
imwrite(imageOutput, outputFilename);

% save out a few example neighborhoods
targetPixelLabel = gfl_ind(5000); % This should be the label of the 500th 'valid' pixel

% Create an RGB image where each pixel starts off as black [0 0 0]
imageOutput = zeros([size(mask), 3]);

% Make in-mask pixels white
for x = 1:size(mask, 1)
    for y = 1:size(mask, 2)
        if mask(x, y)  % Check if the pixel is in the mask
            imageOutput(x, y, :) = [1 1 1]; % White
        end
    end
end

% Color the target pixel red
[targetX, targetY] = find(labMaskMat == targetPixelLabel); % Find the coordinates of the target pixel
imageOutput(targetX, targetY, :) = [1 0 0]; % Red

% Color the neighbors blue
neighbors = gfl{5000}; % Get the list of neighbors from the 500th row
for i = 1:length(neighbors)
    [nx, ny] = find(labMaskMat == neighbors(i)); % Find the coordinates of each neighbor
    imageOutput(nx, ny, :) = [0 0 1]; % Blue
end

% Save the image
outputFilename = '~/testMaskMice_Neighbors5000.png'; % Specify the path and filename for output
imwrite(imageOutput, outputFilename);

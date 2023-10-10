function Vert_ASL(subj,sesh)
% compute the vertex-wise anterior/superior/left vector from SigStreams

% load in sigStreams
% load in surface for euclidean locations
% initialize vector score matrix (2562 x 3 or subset of 2562, x 3)
% for each vertex
for v=vertices
	% get euclidean location of this vertex (x,y,z)
	% get row of interest from sigstreams
	% initialize vector score
	% find each vertex with a nonzerovalue in this row
	nonzvs=find(
	for v2=nonzvs
		% get current cell
		% get euclidean location of cell
		% get vector from v to v2
		% multiply by counts, can later ammend to sqrt of counts
		% add derived vector to vector score
	end
	% populate vector score matrix
end
% saveout vector score matrix

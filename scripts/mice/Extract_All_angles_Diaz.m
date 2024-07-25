% addpaths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% mouse list
mList={'m7531','m7576','m7520','m7522','m7590'};
% for each mouse extract Relative angles
for i=1:length(mList)
	m=mList{i}
	% for each run
	for r=1:6
		Extract_RelativeAngles_mice_Diaz(m,r)
	end
end
disp('donezo')

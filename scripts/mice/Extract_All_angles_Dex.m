% addpaths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% mouse list
mList={'m1','m824','m961','m962','m964'};
% for each mouse extract Relative angles
for i=1:length(mList)
	m=mList{i}
	% for each run
	for r=1:6
		Extract_RelativeAngles_mice_Dex(m,r)
	end
end
disp('donezo')

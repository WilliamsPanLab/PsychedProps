% addpaths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% set IS as FB
FB='IS'
% mouse list
mList={'m2000','m7507','m7520','m7522','m7589','m7594'};
%mList={'m7520','m7522','m7589','m7594'};
% for each mouse extract Relative angles
for i=1:length(mList)
	m=mList{i}
	% for each run
	for r=1:6
		Extract_RelativeAngles_mice(m,r,FB)
	end
end
disp('donezo')

function sbjData = prepareFuncData_cifti_func(sbjListFile,tNum)

% get sbj list
disp('Get subject list...');
fid = fopen(sbjListFile,'r');
sbjList = textscan(fid,'%s');
sbjList = sbjList{1};
fclose(fid);

disp('read flat data...');
sbjData = cell(length(sbjList),1);
for si=1:length(sbjList)
    fileName = sbjList{si};
    disp([num2str(si),'. ',fileName]); 
    if exist(fileName,'file')
	viddata=h5read(fileName, '/processed_data');
	% load in mask to recover indices of each valid points
	%formask=h5read(fileName, '/atlas');
	% NOTE MASK HAS TO BE DERIVED FROM STATIC MOUSE: MOUSE STRUCTURE HAS DIFFERENT SUM OF ONES FROM MOUSE-TO-MOUSE: AMENDED 4/15/24 AP in DS_data.ipynb
	formask=h5read('/scratch/users/apines/p50_mice/proc/20200228/thy1gc6s_0p3mgkg_m2000_preLSD0p3mgkg_1/masked_dff_DS_BP_Smoothed.h5','/mask');
	mask = cellfun(@(x) strcmpi(x, 'TRUE'), formask);
	% Find indices of the valid mask points
	[rows, cols] = find(mask);
	% intialize subj data as number of valid points by time	
	resultData = zeros(length(rows), tNum);
	for i = 1:length(rows)
	    resultData(i, :) = squeeze(viddata(rows(i), cols(i), 1:tNum)).';
	end
	size(resultData);
	% changed out 4/24/24 AP for flat mice
	sbjData{si}=resultData';
	%sbjData{si} = (data.(1:91282,:))';  % change from 12/21/20 to run on surface only - (1:59412 added) % changed back 10/6/22

%	sbjData{si} = (cii.cdata(1:59412,:))';
    else
        disp('  this file does not exist !');
        continue;
    end
end

nonEmpty = ones(length(sbjData),1);
for si=1:length(sbjData)
    if isempty(sbjData{si})
        nonEmpty(si,1) = 0;
    end
end
sbjData = sbjData(logical(nonEmpty),1);

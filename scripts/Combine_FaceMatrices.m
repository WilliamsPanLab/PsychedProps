% read in subject dosage correspondence, has to be before addpath for some silly reason
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');

% combine angular distance from reference streams across subjects
subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
sessions = {'ses-00'};
tasks = {'rs1', 'rs2','gambling','emotion','wm'};
% initialize a histc counts matrix for each stream: 17 angular bins (0-10 degrees, 10-20 degrees, etc.)
% last dimension (3) is for the 3 different sessions
DMNcounts=zeros(18,14,3);
Dorsalcounts=zeros(18,14,3);
Ventralcounts=zeros(18,14,3);
Insularcounts=zeros(18,14,3);
MedPostcounts=zeros(18,14,3);
MedAntcounts=zeros(18,14,3);

iterator=0
% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% aggregate each subject, session, and task for each streams
for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17]
	% add 1 to iterator
	iterator=iterator+1
	% get subject name
	subj=subjList{s}
	% get session info
	seshInfo=subSeshDose{s,2:5};
	% for placebo, 80mg, and 120mg
	seshArray={seshInfo{2} seshInfo{3} seshInfo{4}};
	% make an iterator for session
	seshIterator=0;
	for sessioncell=seshArray
		% update session iterator
		seshIterator=seshIterator+1;
		% initialize DMN, Dorsal Stream, Ventral Stream, Insular Stream, Medial Posterior, and Medial Anterior streams
		DMN = [];
		Dorsal = [];
		Ventral = [];
		Insular = [];
		MedPost = [];
		MedAnt = [];
		sesh=sessioncell{1};
		for taskcell=tasks
			task=taskcell{1};
			% load in data
			fp=['/scratch/users/apines/gp/PropFeats/' subj '_' sesh '_' task '_faceMatrix.mat'];
			% if file exists
			if exist(fp,'file')
			faceMatrix=load(fp).faceMatrix;
			% stream 1, DMN
			st1=faceMatrix(:,1,:);
			% dorsal stream
			st2=faceMatrix(:,2,:);
			% ventral stream
	                st3=faceMatrix(:,3,:);
			% insular stream
                        st4=faceMatrix(:,4,:);
                        % posterior medial stream
                        st5=faceMatrix(:,5,:);
                        % anterior medial stream
                        st6=faceMatrix(:,6,:);
			% insert 1st dimension into DMN vector, 2nd to Dorsal, etc.
			DMN = [DMN st1(:)'];
			% eliminate null/0 values
			DMN = DMN(DMN~=0);
			Dorsal = [Dorsal st2(:)'];
			Dorsal = Dorsal(Dorsal~=0);
			Ventral = [Ventral st3(:)'];
			Ventral = Ventral(Ventral~=0);
			Insular = [Insular st4(:)'];
			Insular = Insular(Insular~=0);
			MedPost = [MedPost st5(:)'];
			MedPost = MedPost(MedPost~=0);
			MedAnt = [MedAnt st6(:)'];
			MedAnt = MedAnt(MedAnt~=0);
			% if it doesnt exist
			else
				fp
				disp('no file found')
			end
		end
		% get counts
		DMNcount=histc(DMN,0:10:180);
		DMNcounts(:,iterator,seshIterator)=DMNcount(1:18);
		Dorsalcount=histc(Dorsal,0:10:180);
        	Dorsalcounts(:,iterator,seshIterator)=Dorsalcount(1:18);
		Ventralcount=histc(Ventral,0:10:180);
        	Ventralcounts(:,iterator,seshIterator)=Ventralcount(1:18);
		Insularcount=histc(Insular,0:10:180);
        	Insularcounts(:,iterator,seshIterator)=Insularcount(1:18);
		MedPostcount=histc(MedPost,0:10:180);
        	MedPostcounts(:,iterator,seshIterator)=MedPostcount(1:18);
		MedAntcount=histc(MedAnt,0:10:180);
        	MedAntcounts(:,iterator,seshIterator)=MedAntcount(1:18);	
	% end session
	end
% end subject
end
% saveout placebo
writetable(table(DMNcounts(:,:,1)),'~/str_pl_DMNcounts.csv');
writetable(table(Dorsalcounts(:,:,1)),'~/str_pl_DScounts.csv');
writetable(table(Ventralcounts(:,:,1)),'~/str_pl_VScounts.csv');
writetable(table(Insularcounts(:,:,1)),'~/str_pl_INScounts.csv');
writetable(table(MedPostcounts(:,:,1)),'~/str_pl_MPcounts.csv');
writetable(table(MedAntcounts(:,:,1)),'~/str_pl_MAcounts.csv');
% save out 80 mg
writetable(table(DMNcounts(:,:,2)),'~/str_m8_DMNcounts.csv');
writetable(table(Dorsalcounts(:,:,2)),'~/str_m8_DScounts.csv');
writetable(table(Ventralcounts(:,:,2)),'~/str_m8_VScounts.csv');
writetable(table(Insularcounts(:,:,2)),'~/str_m8_INScounts.csv');
writetable(table(MedPostcounts(:,:,2)),'~/str_m8_MPcounts.csv');
writetable(table(MedAntcounts(:,:,2)),'~/str_m8_MAcounts.csv');
% save out 120mg
writetable(table(DMNcounts(:,:,3)),'~/str_m12_DMNcounts.csv');
writetable(table(Dorsalcounts(:,:,3)),'~/str_m12_DScounts.csv');
writetable(table(Ventralcounts(:,:,3)),'~/str_m12_VScounts.csv');
writetable(table(Insularcounts(:,:,3)),'~/str_m12_INScounts.csv');
writetable(table(MedPostcounts(:,:,3)),'~/str_m12_MPcounts.csv');
writetable(table(MedAntcounts(:,:,3)),'~/str_m12_MAcounts.csv');

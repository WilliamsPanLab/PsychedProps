% combine angular distance from reference streams across subjects
subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
sessions = {'ses-00'};
tasks = {'rs1', 'rs2'};
% initialize a histc counts matrix for each stream: 17 angular bins (0-10 degrees, 10-20 degrees, etc.)
DMNcounts=zeros(18,14);
Dorsalcounts=zeros(18,14);
Ventralcounts=zeros(18,14);
Insularcounts=zeros(18,14);
MedPostcounts=zeros(18,14);
MedAntcounts=zeros(18,14);
iterator=0
% aggregate each subject, session, and task for each streams
for subcell=subjects
	iterator=iterator+1
	subj=subcell{1}
	% initialize DMN, Dorsal Stream, Ventral Stream, Insular Stream, Medial Posterior, and Medial Anterior streams
	DMN = [];
	Dorsal = [];
	Ventral = [];
	Insular = [];
	MedPost = [];
	MedAnt = [];
	for sessioncell=sessions
		sesh=sessioncell{1};
		for taskcell=tasks
			task=taskcell{1};
			% load in data
			fp=['/scratch/users/apines/gp/PropFeats/' subj '_' sesh '_' task '_faceMatrix.mat'];
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
		end
	end
	% get counts
	DMNcount=histc(DMN,0:10:180);
	DMNcounts(:,iterator)=DMNcount(1:18);
	Dorsalcount=histc(Dorsal,0:10:180);
        Dorsalcounts(:,iterator)=Dorsalcount(1:18);
	Ventralcount=histc(Ventral,0:10:180);
        Ventralcounts(:,iterator)=Ventralcount(1:18);
	Insularcount=histc(Insular,0:10:180);
        Insularcounts(:,iterator)=Insularcount(1:18);
	MedPostcount=histc(MedPost,0:10:180);
        MedPostcounts(:,iterator)=MedPostcount(1:18);
	MedAntcount=histc(MedAnt,0:10:180);
        MedAntcounts(:,iterator)=MedAntcount(1:18);	
end
% saveout
writetable(table(DMNcounts),'~/str_DMNcounts.csv');
writetable(table(Dorsalcounts),'~/str_DScounts.csv');
writetable(table(Ventralcounts),'~/str_VScounts.csv');
writetable(table(Insularcounts),'~/str_INScounts.csv');
writetable(table(MedPostcounts),'~/str_MPcounts.csv');
writetable(table(MedAntcounts),'~/str_MAcounts.csv');

save(['/scratch/users/apines/RefStreamCounts.'],'faceMatrix')

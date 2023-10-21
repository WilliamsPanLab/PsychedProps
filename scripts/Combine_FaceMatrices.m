% combine angular distance from reference streams across subjects
subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
sessions = {'ses-00'};
tasks = {'rs1', 'rs2'};
% initialize a histc counts matrix for each stream: 17 angular bins (0-10 degrees, 10-20 degrees, etc.)
DMNcounts=zeros(17,14);
Dorsalcounts=zeros(17,14);
Ventralcounts=zeros(17,14);
Insularcounts=zeros(17,14);
MedPostcounts=zeros(17,14);
MedAntcounts=zeros(17,14);
% aggregate each subject, session, and task for each streams
for subcell=subjects
	sub=subcell{1};
	% initialize DMN, Dorsal Stream, Ventral Stream, Insular Stream, Medial Posterior, and Medial Anterior streams
	DMN = [];
	Dorsal = [];
	Ventral = [];
	Insular = [];
	MedPost = [];
	MedAnt = [];
	for sessioncell=sessions
		session=sessioncell{1};
		for taskcell=tasks
			task=taskcell{1};
			% load in data
			fp=['/scratch/users/apines/gp/PropFeats/' subj '_' sesh '_' task '_faceMatrix.csv'];
			faceMatrix=csvread(fp);
			% insert 1st dimension into DMN vector, 2nd to Dorsal, etc.
			DMN = [DMN faceMatrix(:,1,:)];
			% eliminate null/0 values
			DMN = DMN(DMN~=0);
			Dorsal = [Dorsal faceMatrix(:,2,:)];
			Dorsal = Dorsal(Dorsal~=0);
			Ventral = [Ventral faceMatrix(:,3,:)];
			Ventral = Ventral(Ventral~=0);
			Insular = [Insular faceMatrix(:,4,:)];
			Insular = Insular(Insular~=0);
			MedPost = [MedPost faceMatrix(:,5,:)];
			MedPost = MedPost(MedPost~=0);
			MedAnt = [MedAnt faceMatrix(:,6,:)];
			MedAnt = MedAnt(MedAnt~=0);
		end
	end
end

% read in subject dosage correspondence, has to be before addpath for some silly reason
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');

% combine angular distance from reference streams across subjects
subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
tasks = {'wm'};
% initialize a histc counts matrix for each stream: 17 angular bins (0-10 degrees, 10-20 degrees, etc.)
% last dimension (3) is for the 3 different sessions
VMmedcounts=zeros(18,14,3);
Dcounts=zeros(18,14,3);
OFCcounts=zeros(18,14,3);
PostMedMotcounts=zeros(18,14,3);
PostMedPCCcounts=zeros(18,14,3);
AntLatMotcounts=zeros(18,14,3);
SupInscounts=zeros(18,14,3);
InfInscounts=zeros(18,14,3);
AudCounts=zeros(18,14,3);
MidCingcounts=zeros(18,14,3);
VMlatcounts=zeros(18,14,3);

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
	% for baseline
	seshArray={seshInfo{1}};
	% make an iterator for session
	seshIterator=0;
	for sessioncell=seshArray
		% update session iterator
		seshIterator=seshIterator+1;
		% initialize streams
		VMmed = [];
		D = [];
		OFC = [];
		PostMedMot = [];
		PostMedPCC = [];
		AntLatMot = [];
		SupIns = [];
		InfIns = [];
		Aud = [];
		MidCing = [];
		VMlat = [];
		sesh=sessioncell{1};
		for taskcell=tasks
			task=taskcell{1};
			% load in data
			fp=['/scratch/users/apines/gp/PropFeats/' subj '_' sesh '_' task '_faceMatrix.mat'];
			% if file exists
			if exist(fp,'file')
			faceMatrix=load(fp).faceMatrix;
			% first non-zero stream is ventral (medial)
			st1=faceMatrix(:,2,:);
			% dorsal stream
			st2=faceMatrix(:,3,:);
			% OFC stream
	                st3=faceMatrix(:,4,:);
			% Premotor
                        st4=faceMatrix(:,5,:);
                        % posterior medial stream
                        st5=faceMatrix(:,6,:);
                        % premotor antlateral
                        st6=faceMatrix(:,7,:);
			% superior insular
			st7=faceMatrix(:,8,:);
			% inferior insular
			st8=faceMatrix(:,9,:);
			% audtiory
			st9=faceMatrix(:,10,:);
			% mid cingulate
			st10=faceMatrix(:,11,:);
			% ventral (lateral)
			st11=faceMatrix(:,12,:);
			% insert 1st dimension into DMN vector, 2nd to Dorsal, etc.
			VMmed = [VMmed st1(:)'];
			D = [D st2(:)'];
			OFC = [OFC st3(:)'];
			PostMedMot = [PostMedMot st4(:)'];
			PostMedPCC = [PostMedPCC st5(:)'];
			AntLatMot = [AntLatMot st6(:)'];
			SupIns = [SupIns st7(:)'];
			InfIns = [InfIns st8(:)'];
			Aud = [Aud st9(:)'];
			MidCing = [MidCing st10(:)'];
			VMlat = [VMlat st11(:)'];
			% eliminate null/0 values
			VMmed = VMmed(VMmed~=0);
			D = D(D~=0);
			OFC = OFC(OFC~=0);
			PostMedMot = PostMedMot(PostMedMot~=0);
			PostMedPCC = PostMedPCC(PostMedPCC~=0);
			AntLatMot = AntLatMot(AntLatMot~=0);
			SupIns = SupIns(SupIns~=0);
			InfIns = InfIns(InfIns~=0);
			Aud = Aud(Aud~=0);
			MidCing = MidCing(MidCing~=0);
			VMlat = VMlat(VMlat~=0);
			% if it doesnt exist
			else
				fp
				disp('no file found')
			end
		end
		% get counts
		VMmedcount=histc(VMmed,0:10:180);
		VMmedcounts(:,iterator,seshIterator)=VMmedcount(1:18);
		Dcount=histc(D,0:10:180);
        	Dcounts(:,iterator,seshIterator)=Dcount(1:18);
		OFCcount=histc(OFC,0:10:180);
        	OFCcounts(:,iterator,seshIterator)=OFCcount(1:18);
		PostMedMotcount=histc(PostMedMot,0:10:180);
        	PostMedMotcounts(:,iterator,seshIterator)=PostMedMotcount(1:18);
		PostMedPCCcount=histc(PostMedPCC,0:10:180);
        	PostMedPCCcounts(:,iterator,seshIterator)=PostMedPCCcount(1:18);
		AntLatMotcount=histc(AntLatMot,0:10:180);
        	AntLatMotcounts(:,iterator,seshIterator)=AntLatMotcount(1:18);	
		SupInscount=histc(SupIns,0:10:180);
		SupInscounts(:,iterator,seshIterator)=SupInscount(1:18);
		InfInscount=histc(InfIns,0:10:180);
                InfInscounts(:,iterator,seshIterator)=InfInscount(1:18);
		Audcount=histc(Aud,0:10:180);
                Audcounts(:,iterator,seshIterator)=Audcount(1:18);
		MidCingcount=histc(MidCing,0:10:180);
                MidCingcounts(:,iterator,seshIterator)=MidCingcount(1:18);
		VMlatcount=histc(VMlat,0:10:180);
                VMlatcounts(:,iterator,seshIterator)=VMlatcount(1:18);
	% end session iteration
	end
% end subject
end
% saveout working memory
writetable(table(VMmedcounts(:,:,1)),'~/str_wm_VMcounts.csv');
writetable(table(Dcounts(:,:,1)),'~/str_wm_Dcounts.csv');
writetable(table(OFCcounts(:,:,1)),'~/str_wm_OFCcounts.csv');
writetable(table(PostMedMotcounts(:,:,1)),'~/str_wm_PMMOTcounts.csv');
writetable(table(PostMedPCCcounts(:,:,1)),'~/str_wm_PMPCCcounts.csv');
writetable(table(AntLatMotcounts(:,:,1)),'~/str_wm_ALMcounts.csv');
writetable(table(SupInscounts(:,:,1)),'~/str_wm_SIcounts.csv');
writetable(table(InfInscounts(:,:,1)),'~/str_wm_IIcounts.csv');
writetable(table(Audcounts(:,:,1)),'~/str_wm_Acounts.csv');
writetable(table(MidCingcounts(:,:,1)),'~/str_wm_MCcounts.csv');
writetable(table(VMlatcounts(:,:,1)),'~/str_wm_VLcounts.csv');

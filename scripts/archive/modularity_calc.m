% add bct path
addpath('~/2019_03_03_BCT/')
% baseline for now
sesh='ses-00'
Modularity=zeros(1,17)
% subj name list
% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')
% for each subj
for s=[1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17]
	% load fc matrix
	subj=char(subjList(s));
	fcPath=['/scratch/groups/leanew1/xcpd_outMDMA_36p_despike_bp/xcp_d/' subj '/' sesh '/func/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_atlas-Schaefer417_den-91k_measure-pearsoncorrelation_conmat.txt'];
	fc=load(fcPath);
	% caclulate modularity
	[ci,q]=modularity_und(fc);
	Modularity(s)=q;
end
% save out
csvwrite('~/Modularity_MDMAbv.csv',Modularity)

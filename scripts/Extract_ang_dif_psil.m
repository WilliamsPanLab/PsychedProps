% extract network props from each session

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% get subj list
subjPrefix=repmat('PS',11,1);
subjSuffix=["02","03","16","18","19","21","24","93","96","98","99"];
subjList=strcat(subjPrefix,subjSuffix')
% for collecting each task
for task=["rs1" "rs2"]
	% initialize output table: maximum of 8 observations per condition per subj (11 subjs x 8 possible, although much lower in most instances) 
	outDF=zeros(88,20);
	% extra column just for subj name
	SubjNameCol=cell(88,1);
	% set common fp
	commonFP=['/scratch/users/apines/data/psil/'];
	% for each subj except 2
	for s=[2 3 4 5 6 7 8 9 10 11];
		disp(s)
		% for up to 8 iterations
		for i=1:8
			% populate subject name column
			SubjNameCol((s+(i*11))-11)=cellstr(subjList(s));
	        	% grab baseline
			bvFP=[commonFP subjList(s) '/Baseline' num2str(i) '/' subjList(s) '_Baseline' num2str(i) '_' task '_k4_Prop_Feats_gro.csv'];
	        	bvFP=strjoin(bvFP,'');
			% p = between
			pFP=[commonFP subjList(s) '/Between' num2str(i)  '/' subjList(s) '_Between' num2str(i) '_' task '_k4_Prop_Feats_gro.csv'];
	        	pFP=strjoin(pFP,'');
			% m1 = after
	        	m1FP=[commonFP subjList(s) '/After' num2str(i) '/' subjList(s) '_After' num2str(i) '_' task '_k4_Prop_Feats_gro.csv'];
	        	m1FP=strjoin(m1FP,'');
			% m2 = drug (to be decoded post-hoc in r)
	        	m2FP=[commonFP subjList(s) '/Drug' num2str(i) '/' subjList(s) '_Drug' num2str(i) '_' task '_k4_Prop_Feats_gro.csv'];
	        	m2FP=strjoin(m2FP,'');
			% load in baseline csvs
			if exist(bvFP,'file')
				bv=readmatrix(bvFP);
				% s+i*11-11 to insert into subj-sesh row 
				outDF(((s+(i*11))-11),1)=bv(1,2);
				outDF(((s+(i*11))-11),2)=bv(2,2);
				outDF(((s+(i*11))-11),3)=bv(3,2);
				outDF(((s+(i*11))-11),4)=bv(4,2);
				% load in remaining frames
                		bvCSIfp=[commonFP subjList(s) '/Baseline' num2str(i) '/' subjList(s) '_Baseline' num2str(i) '_task-' task '_ValidSegments_Trunc.txt'];
                		bvCSIfp=strjoin(bvCSIfp,'');
                		CSI = importdata(bvCSIfp);
                		% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
                		numTRsVS=sum(CSI(:,2));
                		outDF(((s+(i*11))-11),17)=numTRsVS;	
			else
				disp([bvFP ' not found'])
			end
			% load in Between csvs
			if exist(pFP,'file')
				p=readmatrix(pFP);
				% extract in group consensus atlas: placebo
				outDF(((s+(i*11))-11),5)=p(1,2);
				outDF(((s+(i*11))-11),6)=p(2,2);
				outDF(((s+(i*11))-11),7)=p(3,2);
				outDF(((s+(i*11))-11),8)=p(4,2);
				pCSIfp=[commonFP subjList(s) '/Between' num2str(i) '/' subjList(s) '_Between' num2str(i) '_task-' task '_ValidSegments_Trunc.txt'];
        		        pCSIfp=strjoin(pCSIfp,'');
       				CSI = importdata(pCSIfp);
                		numTRsVS=sum(CSI(:,2));
                		outDF(((s+(i*11))-11),18)=numTRsVS;
			else
                                disp([pFP ' not found'])
                        end
			% load After
			if exist(m1FP,'file')
				m1=readmatrix(m1FP);
				% extact in group consensus atlas: mdma 1
				outDF(((s+(i*11))-11),9)=m1(1,2);
				outDF(((s+(i*11))-11),10)=m1(2,2);
				outDF(((s+(i*11))-11),11)=m1(3,2);
				outDF(((s+(i*11))-11),12)=m1(4,2);
        		        m1CSIfp=[commonFP subjList(s) '/After' num2str(i) '/' subjList(s) '_After' num2str(i) '_task-' task '_ValidSegments_Trunc.txt'];
                		m1CSIfp=strjoin(m1CSIfp,'');
                		CSI = importdata(m1CSIfp);
               			numTRsVS=sum(CSI(:,2));
                		outDF(((s+(i*11))-11),19)=numTRsVS;
			else
                                disp([m1FP ' not found'])
			end
			% load Drug
			if exist(m2FP,'file')
				m2=readmatrix(m2FP);
				% extract in group consensus atlas: mdma 2
				outDF(((s+(i*11))-11),13)=m2(1,2);
				outDF(((s+(i*11))-11),14)=m2(2,2);
				outDF(((s+(i*11))-11),15)=m2(3,2);
				outDF(((s+(i*11))-11),16)=m2(4,2);
        		        m2CSIfp=[commonFP subjList(s) '/Drug' num2str(i) '/' subjList(s) '_Drug' num2str(i) '_task-' task '_ValidSegments_Trunc.txt'];
                		m2CSIfp=strjoin(m2CSIfp,'');
              		 	CSI = importdata(m2CSIfp);
                		numTRsVS=sum(CSI(:,2));
                		outDF(((s+(i*11))-11),20)=numTRsVS;
			else
                                disp([m2FP ' not found'])
                        end
		end
	end
	% save out matrix
	writematrix(outDF,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_Psil_propsMerged.csv'],''))
	% and subject ID column
	writetable(table(SubjNameCol),strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_Psil_propsMerged_subjOrder.csv'],''))
end

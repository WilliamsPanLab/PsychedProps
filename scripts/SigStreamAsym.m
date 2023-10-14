test=load('/scratch/users/apines/SimStreams/sub-MDMA001_ses-00_sigStreams.mat')
% initialize vertexwise values
vertvals=zeros(1,length(test.sig_Streams));
% for each vertex value
for i = 1:length(test.sig_Streams)
	vertvals(i)=sum(test.sig_Streams(i,:)-test.sig_Streams(:,i)');
end

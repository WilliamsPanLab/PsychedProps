% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1:numSpins
	% new networks for revisions
        Extract_ang_dif_Spun_VIS(n);
        Extract_VISMag_dif_Spun(n);
        Extract_VISSeg_dif_Spun(n);
        Extract_TAutocor_dif_Spun_VIS(n);
end

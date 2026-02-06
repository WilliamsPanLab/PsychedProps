% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1000:2000
	% new networks for revisions
        Extract_ang_dif_lsd_Spun_VIS(n);
        Extract_VISMag_dif_lsd_Spun(n);
        Extract_VISSeg_dif_lsd_Spun(n);
        Extract_TA_dif_lsd_Spun_VIS(n);
	n
end

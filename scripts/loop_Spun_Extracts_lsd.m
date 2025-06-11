% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1:numSpins
	Extract_ang_dif_lsd_Spun(n);
	Extract_DMNMag_dif_lsd_Spun(n);
	Extract_DMNSeg_dif_lsd_Spun(n);
	Extract_TA_dif_lsd_Spun(n);
	n
end

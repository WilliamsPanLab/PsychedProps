% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1684:2000
	Extract_ang_dif_lsd_Spun(n);
	Extract_DMNMag_dif_lsd_Spun(n);
	Extract_DMNSeg_dif_lsd_Spun(n);
	Extract_TA_dif_lsd_Spun(n);
	n
end
for n=1:1683
	Extract_DMNSeg_dif_lsd_Spun(n);
	n
end

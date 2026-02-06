% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1000:2000
        Extract_ang_dif_lsd_Spun_MOT(n);
        Extract_MOTMag_dif_lsd_Spun(n);
        Extract_MOTSeg_dif_lsd_Spun(n);
        Extract_TA_dif_lsd_Spun_MOT(n);
	n
end

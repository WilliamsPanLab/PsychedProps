% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1:numSpins
        Extract_ang_dif_psil_Spun_MOT(n);
        Extract_MOTMag_dif_psil_Spun(n);
        Extract_MOTSeg_dif_psil_Spun(n);
        Extract_TA_dif_psil_Spun_MOT(n);
	n
end

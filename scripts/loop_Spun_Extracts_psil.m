% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1000:numSpins
	Extract_ang_dif_psil_Spun(n);
	Extract_DMNMag_dif_psil_Spun(n);
	Extract_DMNSeg_dif_psil_Spun(n);
	Extract_TA_dif_psil_Spun(n);
end

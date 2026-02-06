% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1000:numSpins
	Extract_ang_dif_Spun(n);
	Extract_DMNMag_dif_Spun(n);
	Extract_DMNSeg_dif_Spun(n);
	Extract_TAutocor_dif_Spun(n);
end

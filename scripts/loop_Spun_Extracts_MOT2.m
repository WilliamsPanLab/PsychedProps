% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1:numSpins
        Extract_ang_dif_Spun_MOT(n);
        Extract_MOTMag_dif_Spun(n);
        Extract_MOTSeg_dif_Spun(n);
        Extract_TAutocor_dif_Spun_MOT(n);
end

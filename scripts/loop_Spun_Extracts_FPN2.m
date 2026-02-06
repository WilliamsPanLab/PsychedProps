% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1:numSpins
        Extract_ang_dif_Spun_FPN(n);
        Extract_FPNMag_dif_Spun(n);
        Extract_FPNSeg_dif_Spun(n);
        Extract_TAutocor_dif_Spun_FPN(n);
end

% loop over spins to produce equivalent outputs as true measures
numSpins=2000;
% this runs for mdma
for n=1:numSpins
        Extract_ang_dif_psil_Spun_FPN(n);
        Extract_FPNMag_dif_psil_Spun(n);
        Extract_FPNSeg_dif_psil_Spun(n);
        Extract_TA_dif_psil_Spun_FPN(n);
end

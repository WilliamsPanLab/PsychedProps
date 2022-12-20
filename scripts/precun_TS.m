function precun_TS(subj,sesh,infn,outfn)
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti'));
% load ROIs
roi = read_cifti(['/oak/stanford/groups/leanew1/users/apines/maps/Schaefer2018_100Parcels_7Networks_order.dscalar.nii']);
precunL = find(roi.cdata == 50);
precunR = find(roi.cdata == 100);
% load the time series
timeSeries = read_cifti(infn);
% extract cortex
timeSeriesCort=%%%% LEFT OFF HERE, MOVING ON FROM OPAQUE FT_ LIBRARY
tmpL = mean(timeSeries.cdata(precunL,:));
tmpR = mean(timeSeries.cdata(precunR,:));
% saveout as csv
csvwrite(outfn, table(tmpL,tmpR));


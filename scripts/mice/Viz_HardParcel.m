function Viz_HardParcel(initFP,k)
% load in mask to repopulate consensus onto 133x139 grid
fn='/scratch/users/apines/p50_mice/proc/20200228/thy1gc6s_0p3mgkg_m2000_preLSD0p3mgkg_1/masked_dff_DS_BP_Smoothed.h5'
formask=h5read(fn, '/mask');
mask = cellfun(@(x) strcmpi(x, 'TRUE'), formask);

% load in init
res=load(initFP);

% extract init v
[ ~ , HP]=max(res.initV,[],2);

% viz with imshow and good ol' jet
outgrid=zeros(200,208);
outgrid(mask)=HP;
cmap = jet(max(HP));  % Create a colormap with as many colors as the max value in HP
RGB = ind2rgb(outgrid, cmap);
outputfilename=['~/HP_' num2str(k) '.png'];
% save image out
imwrite(RGB, outputfilename)


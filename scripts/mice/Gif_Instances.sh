subj=$1
# specifiy frequency band as argument 2
FB=$2
count=0

# for 100 pngs in this mouse dir
for f in `ls -v /scratch/users/apines/mouseViz/${subj}/vecField_${FB}_*.png | head -n +350`
do
printf -v counts "%02d" $count
# symbolic link to compile into one gif? takes up less space
ln -s $f /scratch/users/apines/mouseViz/${subj}/frame_${counts}.png
count=`expr $count + 1`
done

ffmpeg \
  -framerate 15 \
  -i  /scratch/users/apines/mouseViz/${subj}/frame_%02d.png \
  -vf scale=600:-1 \
  /scratch/users/apines/mouseViz/${subj}_${FB}.gif \
;


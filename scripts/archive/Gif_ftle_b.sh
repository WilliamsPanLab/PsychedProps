subj=$1
# specify run as argument 2
sesh=$2
count=0

# for 100 pngs in this mouse dir
for f in `ls -v /scratch/users/apines/${subj}_${sesh}_t*_b.png | head -n +100`
do
printf -v counts "%02d" $count
# symbolic link to compile into one gif? takes up less space
ln -s $f /scratch/users/apines/mouseViz/${subj}_${sesh}_frame_${counts}_b.png
count=`expr $count + 1`
done

ffmpeg \
  -framerate 15 \
  -i  /scratch/users/apines/mouseViz/${subj}_${sesh}_frame_%02d_b.png \
  -vf scale=600:-1 \
  /scratch/users/apines/mouseViz/${subj}_${sesh}_b.gif \
;


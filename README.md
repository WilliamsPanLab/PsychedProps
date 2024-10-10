# Guide to the code behind "Psychedelics disrupt activity propagations in the default mode network of humans and mice"

This document outlines the steps and methods used in the project. Below is a structured guide for image processing, derivations, and analyses. All image processing was run in a Linux environment using a SLURM cluster for high-compute jobs. In this context, sbatch refers to submitting a job to the SLURM job scheduler. Note that fmriprep and xcpd calls utilize their singularity images, which need to be installed locally. In addition to fmriprep and xcpd, a TON of this code leverages tools built by other people. A summary table is available at the end of this markdown. 

I'll occasionaly refer to study 1, study 2, and study 3. Study 1 is our MDMA sample, 2 is psilocybin, and 3 is LSD/mice.

This is non-comprehensive, but you might find [sbatch_OpFl.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_OpFl.sh) (a parent script that sequentially launches individual processing scripts on scan sessions) useful for further orienting yourself to order-of-operations. If you are interested in cortical propagations broadly, you might also find [the replication guide](https://github.com/PennLINC/DevProps) from our previous [paper](https://www.sciencedirect.com/science/article/pii/S0896627323000387?via%3Dihub) useful. 

## 1. Preprocessing

### 1A. fMRI Preprocessing
  This is the [fmriprep](https://fmriprep.org/en/stable/) call used. These derivatives become the inputs of TSNR (temporal signal to noise ratio) scripts as well as xcpd (eXtensible connectivity pipeline) scripts. It is written to run on a SLURM cluster (Simple Linux Utility for Resource Management): - [full_fmriprep.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/full_fmriprep.sh)
  
  After fmriprep has ran, this script is the first in deriving a TSNR mask. It simply masks out the brain from a nifti. This provides an out-of-brain comparison for the magnitude of signal fluctuations: [TSNR_mask_0_antimask.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/TSNR_mask_0_antimask.sh)
  
  This uses derivatives to calculate SNR for each individual subject: [TSNR_mask_1_ExtractSNR_subjectwise.py](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/TSNR_mask_1_ExtractSNR_subjectwise.py)
  
  This script then combines average SNR maps across participants for a global average map. [TSNR_mask_2_Combine_SNRmaps.py](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/TSNR_mask_2_Combine_SNRmaps.py)

  Next is a quick downsampling of this mask. We'll need it in fsaverage5 space for regularized non negative matrix factorization. Here's that [script, DS_surf_TSNR_fs5.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DS_surf_TSNR_fs5.sh)

  Next, this script combines the averaged left and right hemisphere TSNR maps. The bottom 15 percent of voxels (lowest TSNR) are masked out for all subsequent analyses. Note that psilocybin TSNR maps indicated that psilocybin study data all met the 15th percentile calculated from the MDMA study, likely due to use of multiecho. The same mask is applied to that dataset for equivalence: [TSNR_mask_3_avSNR-mw_to_label](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/TSNR_mask_3_avSNR-mw_to_label.m)

  Now all we need to do is downsample the mask to appropriate resolution for optical flow, fsaverage4. Here's that [script](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DS_surf_TSNR.sh). Shoutout to [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) (wb_command) for making this vanilla and simple.

  That's it for SNR. Now, let's implement the [xcpd](https://xcp-d.readthedocs.io/en/latest/) sbatching script (useful for head-motion-prone scan protocols, runs on fmriprep outputs): [sbatch_xcpd.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_xcpd.sh)

  The WashU crew was kind enough to send their data over already processed. More information on their dataset is available [here](https://www.nature.com/articles/s41586-024-07624-5).For full transparency, the code used to download the images and their associated files are linked below:

  *Download psilocybin data, text files*: [DL_psilo_txts.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DL_psilo_txts.sh)
  
  *Download psilocybin data, framewise displacement files* - [DL_psilo_fd.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DL_psilo_fd.sh)
  
  *Download psilocybin data: neuroimages* - [DL_psilo.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DL_psilo.sh)

  OK, now we have all our processed data. To make it suitable for NMF and optical flow, we are going to motion-mask it. This goes a step further than normal temporal censoring: we are only interested in sequences of neuroimages uninterrupted by head motion. So rather than simply dropping high-motion frames, we'll only retain sequences of at least 8 TRs that are uninterrupted-by-motion. We'll save out the indices of these sequences for optical flow, as we're only interested in activity movement within clean segments (i.e., if frames 9 and 10 have high motion and we retain frames 1-8 and 11-20, we don't want to try and estimate activity movement between 8 and 11, just within the clean segments.) The script that takes care of this is [MotMask.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/MotMask.m). It also converts files to .mgh internally, which are easier for matlab to work with. As will be the case for many scripts below, an [equivalent script](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/MotMask_psil.m) exists for the psilocybin study. In general, the only differences are the input filepath and output file path. For this script in particular, some additional lines of code exist. As we get further in the pipeline, the file derivatives will be more equivalent, requiring fewer differences between MDMA study and Psilocybin Study scripts.

  Nice. Done with fMRI preprocessing!
  
### 1B. Ca2+ Preprocessing
  The first thing to do is to list out sessions to process. This will enable us to loop over sessions to-be-processed within python: [sesh_lister.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/sesh_lister.sh). Note that several different sesh listers exist for different drugs and pre vs. post conditions, but all use the same framework. More information on this dataset is available [here](https://www.nature.com/articles/s41586-020-2731-9).

  Next, we'll obtain a group-level mask so that no mouse has pixels included that other mice don't. We'll combine this with downsampling to saveout a mask that is applicable to processing in downsampled space. As for human data, there are two different downsampled resolutions, with the greater resolution being provided to NMF and the lesser resolution being provided to optical flow for computational tractability.

  *Obtain group-level mask at factor of 2 downsample (NMF)* - [Group_Mask_and_DS_oneHalf.py](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Group_Mask_and_DS_oneHalf.py)

  *Obtain group-level mask at factor of 6 downsample (OpFl)* - [Group_Mask_and_DS_oneSixth.py](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Group_Mask_and_DS_oneSixth.py)

  Just for cleanliness, these scripts are ran independently of downsampling the data for use in NMF and optical flow (even though internal calculations are matched). Here are those scripts. Note they also include gaussian smoothing [for SNR](https://support.brainvoyager.com/brainvoyager/functional-analysis-preparation/29-pre-processing/86-spatial-smoothing) purposes.

  *Downsample by factor of 2 for NMF* - [BP_Smooth_oneHalf.py](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/BP_Smooth_oneHalf.py)

  *Downsample by factor of 6 for Optical flow* - [DS_Smooth_oneSixth_Drug.py](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/DS_Smooth_oneSixth_Drug.py) (adapted for each drug with file path and extensions)

   By the end of this section, you might notice that the mouse data is generally a little less cumbersome to process. That pattern will continue. 

### 1C. DMN Derivation - humans
Before evaluating how DMN function is altered, we have to derive our DMN definition. This is done with regularized non negative matrix factorization. See [this link](https://www.sciencedirect.com/science/article/pii/S1053811917303944?via%3Dihub) for the original paper. In case it's useful context, we've used this approach fairly regularly across half a dozen publications/preprints now: [Cui et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32078800/) [Pines et al. 2022](https://www.nature.com/articles/s41467-022-30244-4) [Shanmugan et al., 2022](https://www.pnas.org/doi/10.1073/pnas.2110416119) [Keller et al., 2023](https://www.nature.com/articles/s41467-023-44087-0) [Hermosillo et al., 2024](https://www.nature.com/articles/s41593-024-01596-5), etc. etc.

*Downsampling to fs5*
  First, the ciftis are downsampled to fsaverage5 to play nicely with NMF. That downsample script, DS_surf_ts_mdma_fs5.sh, is available [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DS_surf_ts_mdma_fs5.sh)
  
*NMF steps*
  Alright, now that we have data in fsaverage5, we'll run it through the NMF pipeline. Basically the steps are to 1) prepare the data for internal NMF processing, 2) create a few candidate solutons using different combinations of brain scans 3) determine the optimal solution of the candidate solutions 4) convert the solution from a matrix to a proper brainmap for further use. Note that both human and mouse NMF utilized extra scans to for greater data volume and reliability. Also note you'll need the helper scripts that come with this code. If you are internal to our lab, you can find the helper scripts at:
> /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts

If you are external to our lab, check out Hongming Li's [repository](https://github.com/hmlicas/Collaborative_Brain_Decomposition).

  NMF script 1: this script will pull some information about the cortical surface you are using (fsaverage5 here) as well as the SNR mask generated in 1A. Here's the link to [Step_1_CreatePrepData.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/Step_1_CreatePrepData.m).
  
  NMF script 2: this script will launch a bunch of parallel matrix decomposition jobs (into a SLURM queue). The internal workflow is to concatenate the functional timeseries data into a single giant matrix, and to decompose the matrix into spatiotemporal components with NMF. Regularization parameters were set a priori, in accordance with previous literature. k=4 components was set in accordance with a nice coarse network solution we previously derived in an independent sample. The reason for using a coarse solution is outlined in the manuscript, but essentially it boils down to inclusive masking to capture propagations entering and exiting the core DMN. Here's the link to [Step_2nd_ParcellationInitialize.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/Step_2nd_ParcellationInitialize.m)
  
  NMF script 3: this script picks out the optimal NMF solution generated by the previous script. Optimal is defined as consistency with other solutions and operationalized via [normalized cuts](https://ieeexplore.ieee.org/abstract/document/868688). This step runs rather quickly and doesn't require a bunch of job submissions and factorization of giant concatenated matrices. Here's the link to [Step_3rd_SelRobustInit.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/Step_3rd_SelRobustInit.m).

  NMF script 4: this script just takes the derived solution and converts it to connectome workbench style file formats (giftis). Simple as. Here's the link to [Step_4_Extract_Group_Atlas.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/Step_4_Extract_Group_Atlas.m).
 
*Downsampling to fs4*
Downsampling the NMF solution to fsaverage4 for optical flow: this step is likely familiar by this point. To downsample fsaverage5 resolution networks to fsaverage4, use [DS_surf_Networks_fs5tofs4.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/DS_surf_Networks_fs5tofs4.sh).

*func.gii to .mat*
The last step is to convert this brainmap back into a plain ol' matrix. This makes matlab less fussy down the road when loading it back in. You can use [Netgiis_2_mat.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/Netgiis_2_mat.m) is pretty straightforward though.

### 1C. DMN Derivation - mice
  Thankfully, we can use the same code to derive functional networks in mice with some slight adaptations. The adaptations are pretty in-the-weeds but pretty simple at a conceptual level. Basically NMF uses an abutment matrix for spatial regularization, where all points in space are arranged in an adjacency matrix with all other points in space. Neighboring voxels are assigned a "1" to reflect sharing a physical edge; this adjacency matrix is leveraged by regularization ultimately making them more likely to be assigned to the same network. The main shift from the original code is from a volumetric image to a 2D brain surface captured via cranial windows. So essentially we just have to remove a dimension. If you look at the equivalent [step 1](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/Step_1st_CreatePrepData_AP_flat.m) script, which we will return to below, it calls [constructW_flat.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/constructW_flat.m) instead of the original script (the original script has to perform across a z-dimension, which is obviated here). This script also has some visualization checks to make sure everything is running as intended. Everything else is extremely similar to the human-brain-adapted NMF we'vea already done by this point. For compelteness, here's an equivalent set of instructions to the human section above.

  NMF script 1: this script will pull some information about the cortical surface you are using (cranial window field-of-view here). You can find it [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/Step_1st_CreatePrepData_AP_flat.m).
  
  NMF script 2: this script will launch a bunch of parallel matrix decomposition jobs (into a SLURM queue). The internal workflow is to concatenate the functional timeseries data into a single giant matrix, and to decompose the matrix into spatiotemporal components with NMF. Regularization parameters did need some adaptation to Ca2+/mouse data from fMRI, but all optical flow and statistical analyses were ran after crystallizing these paramaters on the basis of correspondence with prior accounts of the mouse DMN (see [Whitesell et al](https://pubmed.ncbi.nlm.nih.gov/33290731/) for a good example). k=13 components was set because 1) the DMN takes up a smaller portion of cortical area in mice and 2) the medial wall is not readily visible from the cranial window, meaning that in addition to the DMN occupying less area of mouse cortex we also are not capturing a chunk of the area it does occupy on the medial wall. For those reasons, greater granularity is needed to pull out specific DMN components. The script can be found [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/Step_2nd_ParcellationInitialize_AP_flat.m).
  
  NMF script 3: Again, this script picks out the optimal NMF solution generated by the previous script and itself runs rather quickly. The script can be found [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/Step_3rd_SelRobustInit_AP_flat.m).

  NMF script 3.5: This one didn't get it's own arabic numeral because it's just a plotting script. Nice to triple check that everything is working as intended. Script can be found [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/NMF/Step_3andAHalf_Viz_Consensus_flat.m)

  Great. Now we just have to downsample the mouse DMN to optical flow resolution. It's in python to mirror the same downsampling used for mouse time series data. The script to use for this step is [DS_Smooth_NMFtoOF_resolution.py](https://github.com/WilliamsPanLab/PsychedProps/blob/9744f8f672e1a46085af4f98db34927fa04b9032/scripts/mice/DS_Smooth_NMFtoOF_resolution.py).
  
### 1D. DMN validation

  Ok, these are just some simple checks to make sure our derived DMN maps onto previous accounts of the DMN. Specifically, we'll use the spin test to see if our DMN is more aligned with previously published maps than 10,000 spun nulls. The [spatial permutation approach](https://github.com/spin-test/spin-test) allows us to better-account for autocorrelation in brain maps. There's a subfolder in the scripts directory for 1D, titled [DMN_validate](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DMN_validate/).

  We'll start with the ROIs from Goldstein-Piekarski et al., [Biological Psychiatry 2021](https://www.biologicalpsychiatryjournal.com/article/S0006-3223(21)01437-2/fulltext). More specifically, these ROIs are available on our [lab's github page](https://github.com/WilliamsPanLab/2021-masks). 

First, we have to convert them from volumetric to fs5 so we can use the spin test/have them in the space space as the DMN delineation we just conducted. The script to use for this procedure is [biol_psych_rois_to_fs5.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DMN_validate/biol_psych_rois_to_fs5.sh). 

In parallel, we can "spin" our DMN 10,000 times to get our null maps. That script explicitly uses the spin test repo linked. The script is named [spin_DMN.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DMN_validate/spin_DMN.m).

After the spins are complete (they might take a while, consider sbatching them with something like [sbatch_matlab](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_matlab.sh) (replacing the script ran in matlab with spin_DMN) so you don't have to sit and wait. The equivalent of sbatch is qsub if you're on an SGE (sun grid engine) system. You can use this script to formally evaluate above-chance localization of our current DMN boundaries to those derived in biol. psych. 2021, via t-test. There's a baby script in the main .rmd that will remain to generate figures/stats after running these, but this is the majority of the computation. If you want to spot check the tiny chunk of code in R, you can skip ahead to [Stats_n_Viz.md](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz.md): control-F for '###### spin tests' to find where these analyses were run.

The other previously established DMN map we evaluate is from the epic 35-figure [Yeo et al., 2011](https://journals.physiology.org/doi/full/10.1152/jn.00338.2011). If you don't already know, his lab runs a [great repo](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation), which is where I pulled the .annot files from. Because our NMF-DMN maps are already spun 10,000x from the previous comparison, we don't have to respin to generate our null distributions. Just loop over the spins and run the same t-tests in [yeo7_overlap.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DMN_validate/yeo7_overlap.m).

That's it for step 1! As noted above, data were disimilar upon receipt, but are gradually standardized into similar filenames, extensions, variable names, etc. as they go through more and more shared code. This means that the rest of the scripts should be more equivalent. 


## 2. Optical Flow Derivations

Optical flow is relatively similar to image registration tools we already commonly use in neuroimaging fields. Just like a lot of image registrations, optical flow seeks to find the deformation field that explains the displacement of signal between two images. In registration, we might use a similar algorithm to find the transformation matrix needed to map one brain image onto another. In optical flow, we instead use this optimization to estimate how signal moves over time between two temporally adjacent images in a timeseries. BOLD signal for humans, and calcium signal for mice. Just like image co-registrations, this is computationally intensive.

For mice/3D data (x and y over time), this is fairly tractable with some downsampling. Specifically We use [NeuroPattToolbox](https://github.com/BrainDynamicsUSYD/NeuroPattToolbox) from the Brain Dynamics Group to delineate activity movement trajectories in our mouse data. More info on this software in the table at the end of this markdown.

This gets a little more complicated for human data, because it's acquired in 4D (x y z over time). The workaround for this is to take advantage of [freesurfer](https://surfer.nmr.mgh.harvard.edu/) which maps BOLD, to the cortical surface, and as an intermediate step, to a sphere. Once the data is in spherical space, we can then leverage a nice [spherical optical flow codebase](https://github.com/lukaslang/ofd) that makes this all computationally tractable on the cortical surface (largely through spherical basis functions). This was originally designed to track the [movement of progenitor cells](https://www.semanticscholar.org/paper/Decomposition-of-optical-flow-on-the-sphere-Kirisits-Lang/c7db8ae08ce2f48513198fa5c724657cf8c0d330) on zebrafish gastrulae, but we have sucessfully adapted it for brain data in our [prior publication](https://pubmed.ncbi.nlm.nih.gov/36803653/). Citations for freesurfer and spherical optical flow can also be found in the table at the end of this markdown.

Once we have our resulant vector fields, which describe the movement of BOLD/Ca2+, we then compare it to the [gradient](https://en.wikipedia.org/wiki/Image_gradient) of the default mode network (NMF) to determine if the direction of activity movement is in the bottom-up or top-down direction. 

### 2A. Running optical flow
  OpFl MDMA: To run optical flow on data from study 1 (MDMA), you can use this [OpFl_mdma.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/OpFl_mdma.m). This will take the .mgh output from the motion mask script, ensure it's aligned with the fsaverag4 spherical cortical surface, and run optical flow between each temporally adjacent image WITHIN low-motion segments. The latter is ensure by loading in the _ValidSegments_Trunc.txt file generated during the motion masking step. The real meat of this script comes between lines 103 and 114, where optical flow is ran via the *of* function. This function is directly pulled from the Kirisits et al. [paper](https://www.semanticscholar.org/paper/Decomposition-of-optical-flow-on-the-sphere-Kirisits-Lang/c7db8ae08ce2f48513198fa5c724657cf8c0d330) and [repository](https://github.com/lukaslang/ofd). We'll populate a matlab struct (us) with the fields vf_left and vf_right to store vector fields describing the motion of BOLD signal between timepoints. Parameters are the same as those provided by default in the code and used in our Neuron paper. These will undergo further processing to get our metrics of interest. There's two scripts to aid in submitting these scripts to the cluster. The first is [sbatch_OpFl.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_OpFl.sh), which runs the whole pipeline for optical flow data for an individual subject session. It's useful to comment certain chunks out if you don't want to run the whole pipeline. The second helper script is [sbatch_all.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_all.sh), which submits sbatch_OpFl.sh for each subject for each session as a slurm job. Optional code is in there to NOT submit a job if a certain file (i.e., and output file indicating that the pipeline ran successfully already) exists.

  OpFl Psilocybin: This [script](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/OpFl_psil.m) is identical, but works with file paths/extensions from the psilocybin processing stream. The equivalent helper scripts are [sbatch_OpFl_psil.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_OpFl_psil.sh) for individual subject sessions, and [sbatch_psil_all.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_psil_all.sh) for launching all subjects/sessions.

  OpFl Mice: Within the "mice" folder, you'll find [Mouse_OpFl.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Mouse_OpFl.m). This script is also equivalent, but not exactly the same. First, it has to deal with the different file structure the mouse data is in, and it uses the 3D (x y time) optical flow [repository](https://github.com/BrainDynamicsUSYD/NeuroPattToolbox) provided by Pulin Gong's group (instead of 4D, x y z time). A lot of the code is commented out because we don't utilize the majority of the very cool feature-extraction code they've compiled for this application. There's also some built-in visualization code if you want to double or triple-check that stuff is processing as expected. There's only one helper script for this one because there's only one session (with 5-6 acquisitions) to submit for each mouse for the main analyses. Just sbatch [sbatch_mouseOpFl.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/sbatch_mouseOpFl.sh) once for each mouse, entering the mouse ID as the only agument. Finally, all mouse optical flow steps are listed for LSD, but equivalent scripts for Diazepam and Dexmedetomidine exist for all of these. They are included in the mouse helper script so they can all be in one place. They run identically, but there are fewer acquisition timepoints for the drugs in some instances. This shouldn't disrupt processing, as scripts are designed to run at the level of individual acquisitions. Individual acquisitions might fail if the data doesn't exist (as one would hope), but the whole processing stream shouldn't be interrupted.


  Human and mouse optical flow both use the same [underlying algorithm](https://en.wikipedia.org/wiki/Horn–Schunck_method). 

### 2B. Magnitudes
  Now that we have optical flow estimations, we can pull out the magnitude of activity displacement. We'll just use some good ol' 2,500 year old math for this one. Specifically we apply Pythagorean theorem by calculating the square root of (x displacement^2 + y displacement^2) as our aggregate magnitude of each vector at each point in space over each point in time. We'll just average this over the entire DMN mask over all timepoints to get a single measurement per scan. 

  Human/spherical optical flow measurements technically start out as 3D (x y and z components for each "between" timepoint), but because we conducted optical flow on the sphere, movement of activity orthogonal to the surface of the sphere is negligible. We take advantage of this redundancy by using cart2sphvec, [a matlab-ordained function](https://www.mathworks.com/help/phased/ref/cart2sphvec.html),to obtain activity movement vectors in a tangent plane (tangential to the spherical surface). We're left with azimuth and elevation, which are equivalent to x and y. This has also been validated previously in our optical flow work in task-fMRI and neurodevelopment.

  For the sake of convenience, note that these extractions (2B and 2C) are integrated into [sbatch_OpFl.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_OpFl.sh) and [sbatch_OpFl_psil.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_OpFl_psil.sh). If you want to run analyses broadly, it might be easier to use these pipeline-level scripts rather than run extractions individually for each instance.

  Extract Magnitudes MDMA: We'll run this for individual subjects and sessions once all optical flow runs for that subject have been completed. This can be done serially for each subject/session individually or at the group-level (running all optical flow delineations and then all vector-field-derived extractions). All it does is get the average magnitude of optical flow vectors within the DMN over the ascquisition. The script for this is [Extract_DMNMag.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Extract_DMNMag.m).
  
  Extract Magnitudes Psil: Again, extremely similar to what we ran for study 1. Just differences in the filepath. The script is called [Extract_DMNMag_psil.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Extract_DMNMag_psil.m).

  Extract Magnitudes Mice: As usual, the script takes some adaptation to fit with the data structure of the mouse data. It's extremely similar under the hood though. sqrt(x^2+y^2) averaged within the DMN over time for each mouse. That script is called [Extract_DMNMag_mice.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Extract_DMNMag_mice.m).
  
  *Aggregating DMN magnitudes for assessment of session-to-session differences*
  
  These script just aggregate the derivatives we've calculated. For MDMA, we'll use [Extract_DMNMag_dif.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Extract_DMNMag_dif.m) to loop over subjects and sessions. If you are running this off of Sherlock (Stanford's slurm cluster), I'd just request an interactive node with sdev, ml matlab, open up matlab with matlab -nodisplay, and run this single script. Should run in <10 seconds. You can then pull the output with either scp or OnDemand to your local machine.

  For psilocybin, we'll use  [Extract_DMNMag_dif_psil.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Extract_DMNMag_dif_psil.m). Again, just a tiny bit more complicated because of the variable number of sessions per participant. Same protocol as above for study 1, but note that many files won't be found! This is because the script is designed to index all *possible* sessions/scans, but subjects all have a different number of sessions and scans per each session. It's a little simpler to script for all *possible* sessions/scans and remove blank rows after than to prescribe the exact number of sessions per subject and scans for each session a priori.
  
  For LSD, we'll use this [Extract_mag_dif_mice.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Extract_mag_dif_mice.m). For LSD there's one mouse missing a 6th session, but it's not as heterogeneously populated as the psilocybin data.
  
### 2C. Bottom-Up Relative Angles
  Extract Relative angles MDMA: This will follow a lot of the same steps as magnitude extractions, and is also called by the same helper scripts. We will again use cart2sphvec to obtain activity movement vectors in a tengent plane, but will instead measure the angle (in degrees) of these movement vectors relative to the gradient of the DMN. The surfaces are loaded in the same way, as are the optical flow vectors. As for magnitude, angles are calculated at each point in the DMN at each timepoint, but are averaged across space and time to yield a single observation per acquisition. Specifically, the percentage of all vectors within the DMN over time that are flowing in the bottom-up direction is saved as our primary measurement. As prior, we define bottom-up as < 90 degrees from the gradient of the DMN (BOLD flowing into the DMN) and top-down as > 90 degrees from the gradient of the DMN (BOLD flowing out of the DMN). That script is called [Extract_RelativeAngles.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Extract_RelativeAngles.m). Note that we are using "gradient" in the [classic calculus definition](https://en.wikipedia.org/wiki/Gradient), not in the sense used in the beautiful work by Paquola, Bernhardt, Valk,  Margulies, Misic, Shafiei, Vogel, Smallwood, too many other brilliant scientists to list, etc. For the rest of this guide, I'll refer to "the gradient of the DMN" as "[nabla](https://en.wikipedia.org/wiki/Nabla_symbol) DMN" for clarity.

  Extract Relative angles Psil: As prior, this is quite similar to the script for MDMA. The only difference is the input and output filepaths. The script is called [Extract_RelativeAngles_psil.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Extract_RelativeAngles_psil.m).
  
  Extract Relative angles mice: This script is simpler (because vectors are already in 2 dimensions), but requires different code than the two above because of the difference in input format. There's lots of visual checks in this one in particular (i.e., steps to visualize intermediate calculations to ensure they are working as expected), but they are commented out for ease of use/computational limits (saving out hundreds to thousands of pngs takes up a lot of space and time). Once the data is loaded in, the operations are the same as in the human fmri scripts. The only difference to note is that the DMN cutoff is set to loadings of .6 rather than .3 to account for the mouse DMNs being two merged NMF components. Here's the script:  [Extract_RelativeAngles_mice.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Extract_RelativeAngles_mice.m)

  *Aggregating Bottom-up angle percentages for assessment of session-to-session differences*

  Again, this uses a parallel structure to what we used for aggregating magnitudes. The script for combining acquisition level measurements of bottom-up percentages for MDMA is [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Extract_ang_dif.m), with a parallel script for psilocybin [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Extract_ang_dif_psil.m). Note it's more complicated for psilocybin because of the non-fixed # of sessions and conditions per participant. For mice, the equivalent script is [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Extract_ang_dif_mice.m)

    
### 2C. Contrast Maps

**2C.I** Humans

*Xue, I'm considering re-doing these magnitude ones in particular to make the plotting fully equivalent across all fMRI scripts (on faces rather than vertices). Won't change results, but maybe you can skip these magnitude contrast maps?*

  To calculate the placebo - drug and no drug - drug contrasts, we'll have to run two other scripts first. The first script averages magnitudes for each suject/session (MDMA [here](https://github.com/WilliamsPanLab/PsychedProps/blob/45ccd2d4f43e71219dd154ea1cab43829151ec77/scripts/OpFl_toVerts.m#L62), psilocyin [here](https://github.com/WilliamsPanLab/PsychedProps/blob/a95be2b4c4c29cbbb93ae411fc31911cdba23ab0/scripts/OpFl_toVerts_psil.m#L62). The second script averages magnitudes over space and time for each subject/session (MDMA [here](https://github.com/WilliamsPanLab/PsychedProps/blob/45ccd2d4f43e71219dd154ea1cab43829151ec77/scripts/Calc_AvgMagnitude.m), psilocybin [here](https://github.com/WilliamsPanLab/PsychedProps/blob/a95be2b4c4c29cbbb93ae411fc31911cdba23ab0/scripts/Calc_AvgMagnitude_psil.m). After calculating these values, we can find the average across subjects/sessions with the AvgMag_x_L/R.mat output files generated by using [this script](https://github.com/WilliamsPanLab/PsychedProps/blob/a95be2b4c4c29cbbb93ae411fc31911cdba23ab0/scripts/Aggregate_AvgMags.m) for MDMA and [this script](https://github.com/WilliamsPanLab/PsychedProps/blob/a95be2b4c4c29cbbb93ae411fc31911cdba23ab0/scripts/Aggregate_AvgMag_psil.m) for psilocybin. These last scripts will print out the images found in figure 2 using [Vis_Vertvec](https://github.com/WilliamsPanLab/PsychedProps/blob/a95be2b4c4c29cbbb93ae411fc31911cdba23ab0/scripts/Vis_Vertvec.m).

  The same procecdure can be applied to generate contrasts maps for % bottom-up. Extract relative angles provides the output needed for [Calc_AvgBup](https://github.com/WilliamsPanLab/PsychedProps/blob/a95be2b4c4c29cbbb93ae411fc31911cdba23ab0/scripts/Calc_AvgBup.m) (MDMA) and [Calc_AvgBup_psil](https://github.com/WilliamsPanLab/PsychedProps/blob/a95be2b4c4c29cbbb93ae411fc31911cdba23ab0/scripts/Calc_AvgBup_psil.m). After calculating these values, we can use scripts equivalent to aggregate_avgMags to produce the same visualizations. The script for doing so is [here](https://github.com/WilliamsPanLab/PsychedProps/blob/a95be2b4c4c29cbbb93ae411fc31911cdba23ab0/scripts/Aggregate_AvgBup.m) for MDMA, and [here](https://github.com/WilliamsPanLab/PsychedProps/blob/a95be2b4c4c29cbbb93ae411fc31911cdba23ab0/scripts/Aggregate_AvgBup_psil.m) for psilocybin. 

**2C.II** Mice

  The scripts are a little simpler for mice because the data is more homogenous, 2d, etc. All aggregation and plotting of average magnitudes in LSD vs nodrug contrasts is in a script called [Aggregate_AvgMag_LSD.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/mice/Aggregate_AvgMag_LSD.m).
  
  Similarly, aggregation and plotting of average % BUP is within a single script that can be found within the script [Aggregate_AvgBup_LSD.m](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Aggregate_AvgBup_LSD.m).
  
  Note there are some additional drugs in the mouse study. For Diazepam, aggregating average magnitudes and plotting is within [Aggregate_AvgMag_Diaz.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/mice/Aggregate_AvgMag_Diaz.m), and the equivalent script for bottom-up % is [Aggregate_AvgBup_Diaz.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/mice/Aggregate_AvgBup_Diaz.m). For Dexmedetomidine, the average magnitudes script is [Aggregate_AvgMag_Dex.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/mice/Aggregate_AvgMag_Dex.m), and the average bottom-up % script is [Aggregate_AvgBup_Dex.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/mice/Aggregate_AvgBup_Dex.m).

## 3. DMN Measurement Derivations

  Great, all of section 3 is a bit more straightforward to calculate. We'll start with FC, then temporal autocrrelation.
  
### 3A. Functional Connectivity (FC)

**3A.I** Humans
  [Extract_DMNSeg](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/Extract_DMNSeg.m) just takes the time series, applies the same masking used in the optical flow pipeline (in terms of the DMN and temporal masking), and calculates the functional connectivity between the DMN and the rest of the cortex. [Here](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/Extract_DMNSeg_psil.m) is the psilocybin version (different filepaths). The easiest way to use these is to uncomment them in the helper scripts [sbatch_OpFl.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_OpFl.sh) and [sbatch_OpFl_psil.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_OpFl_psil.sh)

  The scripts to aggregate the derived FC values is [Extract_DMNSeg_dif.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/Extract_DMNSeg_dif.m) for MDMA, and [Extract_DMNSeg_psil.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/Extract_DMNSeg_psil.m) for psilocybin. The easiest way to run these (on Sherlock if you are at Stanford) is just using an interactive node and running it through command-line matlab. Runs very quickly.

**3A.II** Mice
  We use the same approach in mice, but again have to accomodate the different file formats. Use the script [Extract_DMNSeg_mice](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/mice/Extract_DMNSeg_mice.m). I would use the helper script [sbatch_mouseOpFl.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/sbatch_mouseOpFl.sh), with Extract_DMNSeg_mice uncommented, to run through all of these. That way you'll only have to run it once per mouse and you can run them all in parallel as slurm jobs. Once completed, the script to aggregate extracted DMN segregation values is [Extract_DMNSeg_dif_mice.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/mice/Extract_DMNSeg_dif_mice.m), which again is most-easily run through an interactive node and command-line matlab (matlab -nodisplay) on the computing cluster you are using.

### 3B. Autocorrelation

**3B.I** Humans
  We're following the lead of [this](https://www.nature.com/articles/s41593-023-01299-3) paper because they put forward a case that a lot of effects can be simplified to changes in temporal autocorrelation. Thankfully that does not seem to be the case for this propagation stuff, but we should be thinking about this paper in network analyses broadly (and papers that have the potential to simplify neuroscience rather than complicate it are extremely rare). The operationalization of autocorrelation essentially comes down to the correlation between signal and signal shifted "+1" in time. The script to calculate temporal autocorrelation in the DMN for MDMA is [Extract_TAutoCor.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/Extract_TAutoCor.m), and for psilocybin it is [Extract_TAutoCor_psil.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/Extract_TAutoCor_psil.m). To aggregate the MDMA autocor values, use [Extract_TA_dif.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/Extract_TA_dif.m). To aggregate the psilocybin autocor values, use [Extract_TA_dif_psil.m](https://github.com/WilliamsPanLab/PsychedProps/blob/ccac39b7736991f33c2027a2964eb9d968369a55/scripts/Extract_TA_dif_psil.m) script.

**3B.II** Mice
  [Extract_TAutoCor_mice](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Extract_TAutoCor_mice.m) is pretty simple given the corr(timeseries,timeseries+1) approach, but has a few dozen lines of code just to handle the input filenames. Uses the same approach as in the human data and the same DMN masking as in all other scripts. Will return the average autocorrelation for each run. I would also run this by sbatching the helper script [sbatch_mouseOpFl.sh](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/sbatch_mouseOpFl.sh) for each mouse.
  
  To aggregate the data across runs in mice (for autocorrelation), use [this script](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Extract_TAutoCor_dif_mice.m). Fairly straightforward given all other code preceding this script.

## 4. Main Effects

We are going to run DMN magnitude analyses first, then integration, then autocorrelation, and end with bottom-up because the first three print out metrics used for the AUC curves (ran in the bottom-up scripting). So for a small portion of the bottom-up scripts, magnitude integration and autocorrelation r scripts are a prerequisite. Note that all aggregation scripts should have been run by this point to organize .csvs for r to easily read.

### 4A. Magnitudes
- **4A.I** MDMA
  After running the aggregation scripts, you should be able to download the aggregated files as .csvs to your local machine. I'd reccomend running the r code locally because rstudio is the bomb. The script used to calculate magnitude effects (decreases in DMN magnitude) for study 1 (MDMA) can be found [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz_Mag.Rmd). I'm going to leave the within-rmarkdown comments as standalone instructions because this readme is getting long. Note this script will also saveout magnitude for unified DMN analyses, i.e., the AUC curves in figure 4.
  
- **4A.II** Psilocybin
  This is structured to be parallel to the study 1 (MDMA) analyses, but requires a little extra scripting just to organize all data in an equivalent fashion. The .rmd can be found [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz_psil_Mag.Rmd). As prior, I'm going to leave the within-markdown comments as instructions.
  
- **4A.III** LSD
  In what has become a pattern by this point, the mouse data will be easier to work with. This time it is because of the scan-session structure rather than the 2D vs. 3D data embedding. The markdown that calculates equivalent statistics for mouse LSD data is [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_viz_mice_Mag.Rmd).

### 4B. DMN integration
- **4B.I** MDMA
  Post-aggregation DMN integration evaluation in study 1 is within the bottom-up analysis script (4D.I)
- **4B.II** Psilocybin
  Because the psilocybin scripts take a little more organizing code to get everything in stats-friendly format, the integration and autocorrelation scripts are separated out. You can find the script for evaluating DMN integration [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz_psil_DMNSeg.Rmd). Some variable names and comments are left as they were when copied over from Stats_n_Viz_psil.Rmd, as the script is equivalent but with different input .csvs.
- **4B.III** LSD
  Post-aggregation DMN integration evaluation in study 3 (mice/LSD) is within the bottom-up analysis script (4D.III).
  
### 4C. DMN Autocor
- **4C.I** MDMA
  Post-aggregation DMN autocorrelation evaluation in study 1 is within the bottom-up analysis script (4D.I)
- **4C.II** Psilocybin
  Because the psilocybin scripts take a little more organizing code to get everything in stats-friendly format, the integration and autocorrelation scripts as separated out. You can find the script for evaluating DMN integration [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz_psil_DMNTA.Rmd). Some variable names and comments are left as they were when copied over from Stats_n_Viz_psil.Rmd, as the script is equivalent but with different input .csvs.
- **4C.III** LSD
  Post-aggregation DMN autocorrelation evaluation in study 3 (mice/LSD) is within the bottom-up analysis script (4D.III).

  ### 4D. Bottom-Up Analysis
- **4D.I** MDMA
  Remember this script encapsulates DMN integration/segregation and autocorrelation scripting. It can be found in markdown format [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz.md), or as an .rmd [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz.Rmd).
- **4D.II** Psilocybin
  The equivalent script for psilocybin can be found [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz_psil.md). Note this also has the lasting effects lil' chunk of code (just after line 1,000). The .rmd version is [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz_psil.Rmd)
- **4D.III** LSD
  The equivalent knit markdown for mice can be found [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_viz_mice.md). The .rmd version can be found [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_viz_mice.Rmd)


## 5. Bootstraps and AUC curves
You can find the bootstrap and AUC analyses for MDMA (figure 4) further down the [same markdown file used above](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz.md). The .md file doesn't have line numbers, but you can cntrl-F to "library(pROC)" to find where this section stats. If instead you are looking at [the .Rmd](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz.Rmd), this starts at line 1061.

Same goes for bootstrap and AUC/bootstraps for psilocybin. Check out the [same markdown as prior](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz_psil.md). Corresponding .rmd is [here](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz_psil.Rmd), and AUC stuff starts at line 1037.

Aaaaand same goes for mice. Here's the [link](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_viz_mice.md). [.rmd link](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_viz_mice.Rmd). AUC stuff starts at line 238.

## 6. Self-Report: MDMA
You can find the bootstrap and AUC analyses (figure 4) further down the [same markdown file used above](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz.md). The .md file doesn't have line numbers, but you can cntrl-F to "inter-psychedelic-session to find where this section stats. If instead looking at the [.Rmd](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/Stats_n_Viz.Rmd), this starts at line 1374.

## Appendix: code leveraged from other sources
| Software | Citation | Use(s) in this project|
|----------|----------|-----------------------|
| Non negative matrix factorization adapted for brain data | Li, H. et al. (2017): Large-scale sparse functional networks from resting state fMRI. Neuroimage | NMF for human and mouse DMN delineation |
| Freesurfer | Dale, A. et al. (1999): Cortical surface-baed analysis. I. Segmentation and surface reconstruction. Neuroimage | Cortical surface modeling |
| Connectome Workbench | Marcus, D. et al., (2011) Informatics and data mining tools and strategies for the Human Connectome Project. Front. Neuroinformatics | Image transformation and visualization |
| g_ls.m   | Cui Z., et al. (2013): PANDA: a pipeline toolbox for analyzing brain diffusion images. Front Hum Neurosci | Listing files in a directory when in matlab |
| fmriprep | Esteban, O., et al. (2018): fMRIPrep: a robust preprocessing pipeline for functional MRI. Nat Methods | Preprocessing of fmri data (prior to xcp-d) |
| xcpd | Mehta, K., et al. (2024): XCP-D: A Robust Pipeline for the post-processing of fMRI data. Imaging Neurosci | Post-processing of fmri data |
| Spherical Optical Flow | Kirisits, C., et al. (2013): Decomposition of optical flow on the sphere. International J. on Geomath | Optical flow on spherical surfaces (human data, inflated coritces) |
| "flat" optical flow | Townsend, R. & Gong, P. (2018): Detection and analysis of spatiotemporal patterns in brain activity. PLoS Comp. Biol. | Optical flow on cortical window data (mice, Ca2+) |
| ggplot | Wickham, H. (2016): ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York | Most figures made in R |
| nlme | Pinheiro, J. & Bates, D (2000): Mixed-Effects Models in S and S-PLUS. Springer | Statistical testing in R |

## Thank you x100 to these authors for distributing their code!




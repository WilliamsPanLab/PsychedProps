import nibabel.freesurfer.io as fsio
import numpy as np
import os

# freeSurfer subjects directory on sherlock
fs_dir = '/share/software/user/open/freesurfer/8.1.0/subjects/fsaverage4'
label_dir = os.path.join(fs_dir, 'label')
output_dir = '/home/users/apines/fs4surf'
# mkdir
os.makedirs(output_dir, exist_ok=True)

# define lobe mappings based on Desikan-Killiany parcellation
lobe_mapping = {
    'frontal': [
        'superiorfrontal',
        'rostralmiddlefrontal',
        'caudalmiddlefrontal',
        'parsopercularis',
        'parstriangularis',
        'parsorbitalis',
        'lateralorbitofrontal',
        'medialorbitofrontal',
        'precentral',
        'paracentral',
        'frontalpole',
        'rostralanteriorcingulate',  # Cingulate -> Frontal
        'caudalanteriorcingulate'    # Cingulate -> Frontal
    ],
    'parietal': [
        'superiorparietal',
        'inferiorparietal',
        'supramarginal',
        'postcentral',
        'precuneus',
        'posteriorcingulate',        # Cingulate -> Parietal
        'isthmuscingulate'           # Cingulate -> Parietal
    ],
    'temporal': [
        'superiortemporal',
        'middletemporal',
        'inferiortemporal',
        'bankssts',
        'fusiform',
        'transversetemporal',
        'entorhinal',
        'temporalpole',
        'parahippocampal'
    ],
    'occipital': [
        'lateraloccipital',
        'lingual',
        'cuneus',
        'pericalcarine'
    ]
}

# for each hemi
for hemi in ['lh', 'rh']:
    print(f"Processing {hemi}...")
    
    # aparc annotation
    aparc_path = os.path.join(label_dir, f'{hemi}.aparc.annot')
    labels, ctab, names = fsio.read_annot(aparc_path)
    
    # decode names
    names = [n.decode('utf-8') if isinstance(n, bytes) else n for n in names]
    
    # create new labels array for lobes (init. with 0 = unknown)
    lobe_labels = np.zeros_like(labels)
    
    # create color table for lobes
    lobe_ctab = np.array([
        [0, 0, 0, 0, 0],           # unknown (0)
        [220, 20, 20, 0, 1],       # frontal (1) - red
        [20, 220, 20, 0, 2],       # parietal (2) - green  
        [20, 20, 220, 0, 3],       # temporal (3) - blue
        [220, 220, 20, 0, 4]       # occipital (4) - yellow
    ], dtype=np.int32)
    
    lobe_names = [b'unknown', b'frontal', b'parietal', b'temporal', b'occipital']
    
    # map each ROI to its lobe
    for lobe_idx, (lobe_name, roi_list) in enumerate(lobe_mapping.items(), start=1):
        print(f"  Processing {lobe_name}...")
        for roi_name in roi_list:
            # Find matching ROI in aparc names (case-insensitive partial match)
            matching_indices = [i for i, name in enumerate(names) 
                              if roi_name.lower() in name.lower()]
            
            for roi_idx in matching_indices:
                # Assign all vertices with this ROI label to the lobe
                vertex_mask = (labels == roi_idx)
                lobe_labels[vertex_mask] = lobe_idx
                if np.sum(vertex_mask) > 0:
                    print(f"    Assigned {names[roi_idx]} to {lobe_name}: {np.sum(vertex_mask)} vertices")
    
    # save the lobe annotation
    output_path = os.path.join(output_dir, f'{hemi}.lobes.annot')
    fsio.write_annot(output_path, lobe_labels, lobe_ctab, lobe_names)
    print(f"Saved {output_path}")
    
    # print summary
    print(f"\n{hemi} Summary:")
    for lobe_idx, lobe_name in enumerate(lobe_names):
        if isinstance(lobe_name, bytes):
            lobe_name = lobe_name.decode('utf-8')
        count = np.sum(lobe_labels == lobe_idx)
        print(f"  {lobe_name}: {count} vertices ({100*count/len(lobe_labels):.1f}%)")
    print()

print("Lobe annotations created successfully!")

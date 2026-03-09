# HCP-EP-functional-dysconnectivity
This repo contains code need to run Multivariate Distance Matrix Regression (MDMR) on fMRI data, using cognitive scores as covariate such as age and sex. 

## Dataset
This project uses resting-state functional MRI (rs-fMRI) and behavioral data (version 1.1) from the Human Connectome Project for Early Psychosis (HCP-EP) https://humanconnectome.org/study/human-connectome-project-for-early-psychosis. 
### rs-fMRI 
The rs-fMRI data used are pre-processed using the HCP pipelines https://humanconnectome.org/software/hcp-mr-pipelines. We only included subjects whose T1-weighted structual scans pass the quality control check (on a scale of 0-2; 0 = fail, 1 = questionable, 2 = pass) by 2 independent raters. For each subject, volumetric data are mapped onto the fsLR 32k surface and are concatenated across 2 scanning sessions of both AP and PA directions into one time series. Subject ID of subjects with less than 4 runs of data are saved as "/path/to/missing.csv". The remaining rs-fMRI data is then parcellated using the Schaefer (2018) 400-regions 7-networks atlas "Schaefer2018_7Networks_400.32k.L.label.gii", "Schaefer2018_7Networks_400.32k.R.label.gii" from https://github.com/DiedrichsenLab/fs_LR_32/blob/master/Schaefer_2018.

More infomration on acquring the whole dataset can be found on the HCP-EP website https://humanconnectome.org/study/human-connectome-project-for-early-psychosis. 
### Qualitative (demographic and cognitive)
Subjects with complete demographic and cognitive data of interest are kept. For each cognitive test, subject IDs are filtered to match those in the rs-fMRI list _(lines 23-28 of dist.py)_. 
Qualitative data from HCP-EP are saved as .txt. 

## Workflow
To run the analysis, please do the following:
1. Create a python environment, e.g., cpac_env;
2. Run 
```
module load anaconda | conda activate cpac_env
```
3. Download C-PAC at https://github.com/FCP-INDI/C-PAC/tree/main/CPAC to cpac_env and **keep the original folder name**;
4. Download get_data.R and dist.py in a seperate folder than C-PAC **in the same environment** cpac_env;
5. Get started with get_data.R - this file contains data cleaning and MDMR analysis;
6. Use dist.py to calculate functional connectivity matrix and distance matrix;
You'll notice that some lines of code goes sback and forth between Python and R. This occurs mostly for uploading and downloading files outputted by both files. This is notified in the comments.

CSV files should be in the format: each row is a participant, each column is a variable with two column names (one is the varibale name, the other is a description of the column). For example, for particioant 1, column 1 = participant ID, column 2 = age, column 3 = NIH age-corrected score. The same format applies to the missing ID file.

### Example
The directory structure looks like this:
```
├── get_data
│   ├── clean_phenotypic.R
│   ├── demographics.R
│   └── get_data.R
├── MDMR
│   ├── C_PAC
│   ├── D_flanker_EP.npy
│   ├── D_flanker_HC.npy
│   ├── D_flanker.npy
│   ├── dist_flanker.py
│   ├── flanker_EP_ID.csv
│   ├── flanker_HC_ID.csv
│   ├── flanker_miss_ID.csv
│   └──  subj_non_miss_fMRI.csv
```

# Citation
This is an undergraduate thesis work so there's no publication associated, please cite this GitHub page or the published abstract: Fang, K., Danyluik, M., & Lavigne, K. M. (2025). Understanding Cognitive Impairment in Early Psychosis through Functional Brain Dysconnectivity: A Whole-Brain Voxel-wise Analysis Approach . McGill Science Undergraduate Research Journal, 20(2). https://doi.org/10.26443/msurj.v1i2.323.

# HCP-EP-functional-dysconnectivity
This repo contains code need to run Multivariate Distance Matrix Regression (MDMR) on fMRI data, using cognitive scores as covariate such as age and sex. This project simplifies the workflow of conducting MDMR by reusing the calc_subdists function from [CPAC/cwas/cwas.py](https://github.com/FCP-INDI/C-PAC/blob/main/CPAC/cwas/cwas.py) to calculate the distance matrix and the MDMR package [MDMR package](https://github.com/dmcartor/MDMR/) to conduct the regression using analytical p-values.

MDMR for Connectome Wide Association Studies (CWAS) has been applied to rs-fMRI analysis by Shezad et al. (2014), https://doi.org/10.1016/j.neuroimage.2014.02.024 and Misaki et al. (2018), https://doi.org/10.1016/j.nicl.2018.08.025. Feola et al. (2024), https://doi.org/10.1038/s41380-024-02512-w and Ward et al. (2024), https://doi.org/10.1016/j.biopsych.2024.07.012 used this methods to analyse HCP-EP data.

The mathematics behind MDMR is described in Zapala, & Schork (2012), https://doi.org/10.3389/fgene.2012.00190 and McArtor, Lubke, & Bergeman (2017), https://doi.org/10.1007/s11336-016-9527-8.
## Pre-requisites
Python version 3.13.2, conda 4.12.0, RStudio.

## Dataset
This project uses resting-state functional MRI (rs-fMRI) and behavioral data (version 1.1) from the Human Connectome Project for Early Psychosis (HCP-EP) https://humanconnectome.org/study/human-connectome-project-for-early-psychosis. 
### rs-fMRI 
The rs-fMRI data used are pre-processed using the HCP pipelines https://humanconnectome.org/software/hcp-mr-pipelines. We only included subjects whose T1-weighted structual scans pass the quality control check (on a scale of 0-2; 0 = fail, 1 = questionable, 2 = pass) by 2 independent raters. For each subject, volumetric data are mapped onto the fsLR 32k surface and are concatenated across 2 scanning sessions of both AP and PA phase encoding into one time series. Subject ID of subjects with less than 4 runs of data are saved as "/path/to/missing.csv". The remaining rs-fMRI data is then parcellated using the Schaefer (2018) 400-regions 7-networks atlas "Schaefer2018_7Networks_400.32k.L.label.gii", "Schaefer2018_7Networks_400.32k.R.label.gii" from https://github.com/DiedrichsenLab/fs_LR_32/blob/master/Schaefer_2018.

More infomration on acquring the whole dataset can be found on the HCP-EP website https://humanconnectome.org/study/human-connectome-project-for-early-psychosis. 
### Qualitative
Subjects with complete demographic and cognitive data of interest are kept. This is achieved by matching subject IDs with the ones in the rs-fMRI list _(lines 23-28 of dist.py)_. Qualitative data from HCP-EP are saved as .csv. For example, cognitive data has the format below. Note that age is measured in months.
|Subject_ID | Age | NIH_Age_Corrected |  
|-----------|-----|-------------------|  
| 1001      | 279 | 60                |  
| 1002      | 310 | 50                |  
| 1003      | 250 | 65                |  

Lists of subject IDs are also saved as .csv. missing.csv and subj_non_miss_fMRI.csv have the structure
| 0                               | 
|---------------------------------|
| 1001_01_MR_fsLR_schaefer400.txt |
| 1002_01_MR_fsLR_schaefer400.txt |
| 1003_01_MR_fsLR_schaefer400.txt |

Test-specific subjet ID such as flanker_EP_ID.csv have the structure

|  "x"   | 
|--------|
| "1001" |
| "1002" |
| "1003" |


## Workflow
To run the analysis,
1. Create a python environment using Anaconda, e.g., `conda create -n cpac_env python=3.13.2 anaconda`. 
2. Activate the environment, e.g., `module load anaconda | conda activate cpac_env`.
3. Download C-PAC  https://github.com/FCP-INDI/C-PAC/tree/main/CPAC using `pip install cpac`and **keep the original folder name**.
4. Download get_data.R and dist.py in seperate folders than C-PAC **in the same environment** cpac_env. For example, the project directory structure looks like this:
```
├── get_data
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
5. Get started with get_data.R. This goes through data cleaning and MDMR analysis.
6. Use dist.py to calculate functional connectivity matrix and distance matrix for the Flanker test. To run more than one cognitive tests, e.g., the Flanker test and the Oral Reading Recognition test, name the distance file "dist_flanker.py" ad "dist_orrt.py" and substitute corresponding test-specific information in the code.

Note that some lines of code goes sback and forth between Python and R. This occurs mostly for uploading and downloading files outputted by both files. This is notified in the comments.

# Citation
This is an undergraduate thesis work so there's no publication associated, please cite this GitHub page or the published abstract: Fang, K., Danyluik, M., & Lavigne, K. M. (2025). Understanding Cognitive Impairment in Early Psychosis through Functional Brain Dysconnectivity: A Whole-Brain Voxel-wise Analysis Approach . McGill Science Undergraduate Research Journal, 20(2). https://doi.org/10.26443/msurj.v1i2.323.

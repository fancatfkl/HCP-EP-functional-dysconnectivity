# HCP-EP-functional-dysconnectivity
This repo contains code need to run Multivariate Distance Matrix Regression (MDMR) on fMRI data, using cognitive scores as the covariate (plus age and sex). To run the analysis, please do the following:
1. create a python environment;
2. download C-PAC at https://github.com/FCP-INDI/C-PAC/tree/main/CPAC
3. get started with get_data.R - this file contains data cleaning and MDMR analysis
4. use dist.py to calculate functional connectivity matrix and distance matrix
5. you'll notice that some lines of code goe sback and forth between Python and R. This occurs mostly for uploading and doanlaoding files outputted by both files. 

*Citation*: This is an undergraduate thesis work so there's no publication associated, please cite this GitHub page or the published abstract: Fang, K., Danyluik, M., & Lavigne, K. M. (2025). Understanding Cognitive Impairment in Early Psychosis through Functional Brain Dysconnectivity: A Whole-Brain Voxel-wise Analysis Approach . McGill Science Undergraduate Research Journal, 20(2). https://doi.org/10.26443/msurj.v1i2.323.

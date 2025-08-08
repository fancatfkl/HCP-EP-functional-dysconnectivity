
import os
import nibabel as nb
import numpy as np
import glob
import scipy.stats
import sys
import pandas as pd
import re
sys.path.append(os.path.abspath("C_PAC"))
from CPAC.utils import correlation
from CPAC.utils.typing import ConfigKeyType


datadir = "/path/to/data"
# remove subjects with less than 4 runs of fMRI
all_subj = os.listdir(datadir)
missing_ID = pd.read_csv('/path/to/missing/csv', header = None) # get your list of IDs of missing data
'''
repeat below for all cognitive tests
'''

missing_ID_flanker = pd.read_csv('/path/to/missing/flanker/data') #** get the list of IDs of missing cognitive data from R

missing_names_1 = missing_ID[missing_ID.columns[0]].drop_duplicates().to_numpy() 
    # 9 subjects without full imaging data 

missing_names_2 = missing_ID_flanker[missing_ID_flanker.columns[0]].to_numpy()
missing_names_2 = [str(subj) + '_01_MR' for subj in missing_names_2]

missing_names = np.concatenate((missing_names_1, missing_names_2), axis=0)

pattern = r'\d{4}'
matched_IDs = {re.search(pattern, name).group() for name in missing_names}
flanker_subj = [
    subj for subj in all_subj 
    if re.search(pattern, subj).group() not in matched_IDs
    ] # create a list of subjects for flanker test


'''
total MDMR
'''
# create dictionary
flanker_dict = {
        subjdir : f"/path/to/data/{subjdir}" for
        subjdir in flanker_subj
        }

flanker_files = list(flanker_dict.values())

flanker_data = np.array([
    np.genfromtxt(subject_file).astype('float64').transpose()
    for subject_file in flanker_files
    ])
    
print(flanker_data.shape) # nsubj*nregion*ntimepoints

def calc_subdists(flanker_data): # function to create distance matrix
    nsubs, nregions, _ = flanker_data.shape
    D = np.zeros([nregions, nsubs, nsubs]) # create an empty dataframe first
    for i in range(nregions): # tehn plug in the numbers into the dataframe
        profiles = np.zeros((nsubs, nregions))
        for si in range(nsubs):
            profiles[si] = correlation(flanker_data[si, i], flanker_data[si]) # functional connectivity
        profiles = np.clip(np.nan_to_num(profiles), -0.9999, 0.9999) # ensure it's bounded correctly
        profiles = np.arctanh(np.delete(profiles, i, 1)) # transform
        D[i] = correlation(profiles, profiles) # correlations to calculate distance from

    D = np.sqrt(2.0 * (1.0 - D)) # distance function
    return D
 
D_flanker = calc_subdists(flanker_data)
print(D_flanker.shape) # (400, 158, 158)

#save
np.save("D_flanker.npy", D_flanker) #** use this .npy as Y in get_data.R


'''
between group MDMR
'''
# want to filter the subject_files (subject IDs) according to the group membership

flanker_EP_ID = pd.read_csv('/flanker_EP_ID.csv') # 104 EP

# print(flanker_EP_ID.columns) # Index(['x'], dtype='object') # the elements are loaded as strings

# rename "x" to "ID"
flanker_EP_ID.columns = ['ID'] 

flanker_EP_ID = flanker_EP_ID['ID'] # only keep the ID
flanker_EP_ID = [str(ID) for ID in flanker_EP_ID] # [,,,'1006', '1047', '1072', '1076']
                                                    # <class 'str'>
flanker_HC_ID = pd.read_csv('/flanker_HC_ID.csv') # 54 HC
flanker_HC_ID.columns = ['ID']
flanker_HC_ID = flanker_HC_ID['ID']
flanker_HC_ID = [str(ID) for ID in flanker_HC_ID]

# Subset the dictionary
flanker_EP = [
    subj for subj in flanker_subj
    if re.match(pattern, subj) and re.match(pattern, subj).group() in flanker_EP_ID
]

flanker_EP_dict = {
        subjdir : f"/path/to/data/{subjdir}" for
        subjdir in flanker_EP 
        }

flanker_EP_files = list(flanker_EP_dict.values()) 

flanker_EP_data = np.array([
      np.genfromtxt(subject_file).astype('float64').transpose()
      for subject_file in flanker_EP_files
      ]) # (104, 400, 1640), nsubj*nregion*ntimepoints

flanker_HC = [
    subj for subj in flanker_subj
    if re.match(pattern, subj) and re.match(pattern, subj).group() in flanker_HC_ID
]

flanker_HC_dict = {
        subjdir : f"/path/to/data/{subjdir}" for
        subjdir in flanker_HC 
        }
flanker_HC_files = list(flanker_HC_dict.values())

flanker_HC_data = np.array([
      np.genfromtxt(subject_file).astype('float64').transpose()
      for subject_file in flanker_HC_files
      ]) # (54, 400, 1640)

# Calculate distance matrix
D_flanker_EP = calc_subdists(flanker_EP_data)
print(D_flanker_EP.shape)
D_flanker_EP = np.transpose(D_flanker_EP, (1,2,0)) # transpose to subj*subj*region

D_flanker_HC = calc_subdists(flanker_HC_data)
D_flanker_HC = np.transpose(D_flanker_HC, (1,2,0)) # transpose to subj*subj*region

# Save matrix, then improted into R for MDMR()
np.save("D_flanker_EP.npy", D_flanker_EP)
np.save("D_flanker_HC.npy", D_flanker_HC)

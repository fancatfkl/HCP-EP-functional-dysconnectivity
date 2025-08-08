# install.packages('reticulate')
library(reticulate) # read numpy files
library(tidyverse)
library(table1) # make tables
np <- import('numpy')

# cleaned phenotypic data  ----------------------------------------------------------------------

list_files_pheno <- list.files(path = '/path/to/data')

list_files_pheno #[1] "demographics.csv" "symptom_names.csv" "symptom_sums.csv" "symptoms.csv"


# demographics ------------------------------------------------------------

demographics <- read.csv('/path/to/data') %>%
  as.tibble() # load demographics as a tibble instead of a data.frame

demographics <- demographics %>%
  mutate(
    sex = factor(sex), # factorize the categorical data
    race = factor(race),
    phenotype = factor(phenotype)
  )

demographics$interview_age <- demographics$interview_age/12 # age the age from months to years

EP_demo<- filter(demographics, phenotype == 'Patient') # separate the file according to the phenotype
HC_demo<-filter(demographics, phenotype=='Control')

subj_non_miss_fMRI <- read.csv('/path/to/data', header = T) %>% # **this is manually created in dist.py
  as.data.frame()
length(subj_non_miss_fMRI$X0) # check the numnber of missing data IDs

age.stat <- t.test(EP_demo$interview_age, HC_demo$interview_age)

sex.stat <- chisq.test(table(demographics$phenotype, demographics$sex))

race.stat <- chisq.test(table(demographics$phenotype, demographics$race))

# symptoms ----------------------------------------------------------------
#load in the symptom sums file as a tibble
symp_sums <- read.csv('/path/to/data') %>%
  tibble()

EP_symp <- symp_sums %>% # subset the symp tibble to EP-only
filter(src_subject_id %in% EP_demo$src_subject_id)

# want to find the number of NA's for each column
colSums(is.na(EP_symp)) # note different numbers of NA's for each column, even for the same test

# for each test, want to remove rows with NA's

# flanker test ------------------------------------------------------------

flanker01 <- read.table('/path/to/data', na.strings = c("", "NA")) %>%
  as.data.frame()

flanker <-  flanker01[-c(1:2), c(5, 7, 8, 11)] # extract the ID, age, sex and score

# keep subjects with complete data (imaging and flanker)
# instead of finding indexes, we can straight up match the subject IDs
flanker <- flanker %>%
  filter(V5 %in% subj_non_miss_fMRI$X0) %>%
  mutate(
    V5 = as.numeric(V5), # factorize the categorical data
    V7 = as.numeric(V7)/12,
    V8 = factor(V8),
    V11 = as.numeric(V11)
  ) %>%
  `colnames<-`(c('ID', 'age', 'sex', 'score'))

# print(flanker)

flanker_fMRI_miss <- subj_non_miss_fMRI %>%
  filter(!X0 %in% flanker$ID)
print(flanker_fMRI_miss)
# 13 subjects in the subj_non_miss that are not in flanker

write.csv(flanker_fMRI_miss, '/flanker_miss_ID.csv', row.names = F)
# save the IDs as flanker_miss_ID.csv so that our python dictionary is updated

flanker.bn.EP <- flanker %>%
  filter(ID %in% EP_demo$src_subject_id) # filter out the HC

flanker.bn.EP <- flanker.bn.EP[,-1] # remove the ID column


flanker.bn.HC <- flanker %>%
  filter(ID %in% HC_demo$src_subject_id) # filter out the EP

flanker.bn.HC <- flanker.bn.HC[,-1] %>% # the remove the IDs

## MDMR --------------------------------------------------------------------

D_flanker <-np$load('/D_flanker.npy') #** load from outputs of dist.py
# dim(D_flanker)

# create an empty data frame of 400 by 4
n_rows <- 400
n_cols <- 4
df_flanker <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))

# create a model matrix (optional)
# flanker <- flanker[,-1] # first exclude the ID column
# colnames(flanker) <- c('age', 'sex', 'score')
# flanker <- model.matrix(~ score + age + sex, data = flanker)

# run mdmr for the first region
mdmr(X = flanker, D = D_flanker[1,,])
# %>% summary()

row_index <- 1
for (i in 1:dim(D_flanker)[1]) {
  set.seed(12345) # to be reproducible
  region <- D_flanker[i,,]
  result <- mdmr(X = flanker, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_flanker[row_index, 1] <- F_stats[1]
  df_flanker[row_index, 3] <- F_stats[2]
  df_flanker[row_index, 2] <- p_vals[1]
  df_flanker[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

# extract score_pval
# write.csv(df_flanker$score_pval,
#           'stats_flanker.csv',
#           row.names = T)

### group-wise --------------------------------------------------------------

D_flanker_EP <-np$load('/D_flanker_EP.npy')
D_flanker_HC <-np$load('/D_flanker_HC.npy')

# create an empty data frame of 400 by 4
n_rows <- 400
n_cols <- 4
df_flanker_EP <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))
df_flanker_HC <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))

# loop
row_index <- 1
for (i in 1:dim(D_flanker_EP)[3]) {
  set.seed(12345) # to be reproducible
  region <- D_flanker_EP[,,i]
  result <- mdmr(X = flanker.bn.EP, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_flanker_EP[row_index, 1] <- F_stats[1]
  df_flanker_EP[row_index, 3] <- F_stats[2]
  df_flanker_EP[row_index, 2] <- p_vals[1]
  df_flanker_EP[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

row_index <- 1
for (i in 1:dim(D_flanker_HC)[3]) {
  set.seed(12345)
  region <- D_flanker_HC[,,i]
  result <- mdmr(X = flanker.bn.HC, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_flanker_HC[row_index, 1] <- F_stats[1]
  df_flanker_HC[row_index, 3] <- F_stats[2]
  df_flanker_HC[row_index, 2] <- p_vals[1]
  df_flanker_HC[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

# list sorting ------------------------------------------------------------

lswmt01 <- read.table('/path/to/data', na.strings = c("", "NA"))%>%
  as.data.frame()

lswmt <-  lswmt01[-c(1:2), c(5, 7, 8, 42)] # extract the ID, age, sex and score

sum(is.na(lswmt$V42)) # no missing data
# [1] 0

lswmt <- lswmt %>%
  filter(V5 %in% subj_non_miss_fMRI$X0) %>%
  mutate(
    V5 = as.numeric(V5), # factorize the categorical data
    V7 = as.numeric(V7)/12,
    V8 = factor(V8),
    V42 = as.numeric(V42)
  ) %>%
  `colnames<-`(c('ID', 'age', 'sex', 'score'))

lswmt_fMRI_miss <- subj_non_miss_fMRI %>%
  filter(!X0 %in% lswmt $ID)
print(lswmt_fMRI_miss) # 13 subjects missing
write.csv(lswmt_fMRI_miss, '/lswmt_miss_ID.csv',row.names = F)

lswmt.bn.EP <- lswmt %>%
  filter(ID %in% EP_demo$src_subject_id) # filter out the HC

lswmt.bn.EP <- lswmt.bn.EP[,-1] # remove the ID column

lswmt.bn.HC <- lswmt %>%
  filter(ID %in% HC_demo$src_subject_id) # filter out the EP

lswmt.bn.HC <- lswmt.bn.HC[,-1] # the remove the IDs

write.csv(lswmt.bn.EP$ID, '/lswmt_EP_ID.csv',row.names = F)
write.csv(lswmt.bn.HC$ID, '/lswmt_HC_ID.csv',row.names = F)


## MDMR --------------------------------------------------------------------

D_lswmt <- np$load('/D_lswmt.npy')

# dim(D_lswmt) 400 by 158 by 158 due to the transpose in dist_lswmt.py

# create an empty data frame of 400 by 4
n_rows <- 400
n_cols <- 4
df_lswmt <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))

### group-wise --------------------------------------------------------------

D_lswmt_EP <- np$load('/D_lswmt_EP.npy')
D_lswmt_HC <- np$load('/D_lswmt_HC.npy')

# create an empty data frame of 400 by 4
n_rows <- 400
n_cols <- 4
df_lswmt_EP <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))
df_lswmt_HC <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))


# loop
row_index <- 1
for (i in 1:dim(D_lswmt_EP)[3]) {
  set.seed(12345) # to be reproducible
  region <- D_lswmt_EP[,,i]
  result <- mdmr(X = lswmt.bn.EP, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_lswmt_EP[row_index, 1] <- F_stats[1]
  df_lswmt_EP[row_index, 3] <- F_stats[2]
  df_lswmt_EP[row_index, 2] <- p_vals[1]
  df_lswmt_EP[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

row_index <- 1
for (i in 1:dim(D_lswmt_HC)[3]) {
  set.seed(12345)
  region <- D_lswmt_HC[,,i]
  result <- mdmr(X = lswmt.bn.HC, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_lswmt_HC[row_index, 1] <- F_stats[1]
  df_lswmt_HC[row_index, 3] <- F_stats[2]
  df_lswmt_HC[row_index, 2] <- p_vals[1]
  df_lswmt_HC[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

# oral recognition --------------------------------------------------------

orrt01 <- read.table('/path/to/data', na.strings = c("", "NA"))%>%
  as.data.frame()
orrt <-  orrt01[-c(1:2), c(5, 7, 8, 10)] # extract the ID, age, sex and score
sum(is.na(orrt$V10)) # no missing data
# [1] 0
orrt <-
  orrt %>%
  filter(V5 %in% subj_non_miss_fMRI$X0) %>%
  mutate(
    V5 = as.numeric(V5), # factorize the categorical data
    V7 = as.numeric(V7)/12,
    V8 = factor(V8),
    V10 = as.numeric(V10)
  ) %>%
  `colnames<-`(c('ID', 'age', 'sex', 'score'))

orrt_fMRI_miss <- subj_non_miss_fMRI %>%
  filter(!X0 %in% orrt$ID)

print(orrt_fMRI_miss) # 13 subjects missing
write.csv(orrt_fMRI_miss, '/orrt_miss_ID.csv',row.names = F)
# save the IDs as miss_ID.csv so that our python dictionary is updated

orrt.bn.EP <- orrt %>%
  filter(ID %in% EP_demo$src_subject_id) # use V10
write.csv(orrt.bn.EP$ID, '/orrt_EP_ID.csv',row.names = F)
orrt.bn.EP <- orrt.bn.EP[,-1]

orrt.bn.HC <- orrt %>%
  filter(ID %in% HC_demo$src_subject_id) # use V10
write.csv(orrt.bn.HC$ID, '/orrt_HC_ID.csv',row.names = F)
orrt.bn.HC <- orrt.bn.HC[,-1]

## MDMR --------------------------------------------------------------------

D_orrt <-np$load('/D_orrt.npy')

# dim(D_lswmt) 400 by 158 by 158 due to the transpose in dist_lswmt.py

# create an empty data frame of 400 by 4
n_rows <- 400
n_cols <- 4
df_orrt <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))

### group-wise --------------------------------------------------------------

D_orrt_EP <- np$load('/D_orrt_EP.npy')
D_orrt_HC <- np$load('/D_orrt_HC.npy')

# create an empty data frame of 400 by 4
n_rows <- 400
n_cols <- 4
df_orrt_EP <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))
df_orrt_HC <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))


# loop
row_index <- 1
for (i in 1:dim(D_orrt_EP)[3]) {
  set.seed(12345) # to be reproducible
  region <- D_orrt_EP[,,i]
  result <- mdmr(X = orrt.bn.EP, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_orrt_EP[row_index, 1] <- F_stats[1]
  df_orrt_EP[row_index, 3] <- F_stats[2]
  df_orrt_EP[row_index, 2] <- p_vals[1]
  df_orrt_EP[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

row_index <- 1
for (i in 1:dim(D_orrt_HC)[3]) {
  set.seed(12345)
  region <- D_orrt_HC[,,i]
  result <- mdmr(X = orrt.bn.HC, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_orrt_HC[row_index, 1] <- F_stats[1]
  df_orrt_HC[row_index, 3] <- F_stats[2]
  df_orrt_HC[row_index, 2] <- p_vals[1]
  df_orrt_HC[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

# picture vocab -----------------------------------------------------------

tpvt01 <- read.table('/path/to/data', na.strings = c("", "NA"))%>%
  as.data.frame()

# want the uncorrected standrad score, V12
tpvt <-  tpvt01[-c(1:2), c(5, 7, 8, 12)] # extract the ID, age, sex and score

sum(is.na(tpvt$V12)) # no missing data
# [1] 0

tpvt <-
  tpvt %>%
  filter(V5 %in% subj_non_miss_fMRI$X0) %>%
  mutate(
    V5 = as.numeric(V5), # factorize the categorical data
    V7 = as.numeric(V7)/12,
    V8 = factor(V8),
    V12 = as.numeric(V12)
  ) %>%
  `colnames<-`(c('ID', 'age', 'sex', 'score'))

print(tpvt)

tpvt_fMRI_miss <- subj_non_miss_fMRI %>%
  filter(!X0 %in% tpvt$ID)
print(tpvt_fMRI_miss) # 13 subjects missing
write.csv(tpvt_fMRI_miss, '/tpvt_miss_ID.csv',row.names = F)
# save the IDs as miss_ID.csv so that our python dictionary is updated

tpvt.bn.EP <- tpvt %>%
  filter(ID %in% EP_demo$src_subject_id)
write.csv(tpvt.bn.EP$ID, '/tpvt_EP_ID.csv',row.names = F)
tpvt.bn.EP <- tpvt.bn.EP[,-1]

tpvt.bn.HC <- tpvt %>%
  filter(ID %in% HC_demo$src_subject_id)
write.csv(tpvt.bn.HC$ID, '/tpvt_HC_ID.csv',row.names = F)
tpvt.bn.HC <- tpvt.bn.HC[,-1]


## MDMR --------------------------------------------------------------------

D_tpvt <-np$load('/D_tpvt.npy')

# create an empty data frame of 400 by 4
df_tpvt <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))

### group-wise --------------------------------------------------------------

D_tpvt_EP <-np$load('/D_tpvt_EP.npy')
D_tpvt_HC <-np$load('/D_tpvt_HC.npy')

# create an empty data frame of 400 by 4

df_tpvt_EP <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))

df_tpvt_HC <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))


# loop
row_index <- 1
for (i in 1:dim(D_tpvt_EP)[3]) {
  set.seed(12345) # to be reproducible
  region <- D_tpvt_EP[,,i]
  result <- mdmr(X = tpvt.bn.EP, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_tpvt_EP[row_index, 1] <- F_stats[1]
  df_tpvt_EP[row_index, 3] <- F_stats[2]
  df_tpvt_EP[row_index, 2] <- p_vals[1]
  df_tpvt_EP[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

row_index <- 1
for (i in 1:dim(D_tpvt_HC)[3]) {
  set.seed(12345)
  region <- D_tpvt_HC[,,i]
  result <- mdmr(X = tpvt.bn.HC, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_tpvt_HC[row_index, 1] <- F_stats[1]
  df_tpvt_HC[row_index, 3] <- F_stats[2]
  df_tpvt_HC[row_index, 2] <- p_vals[1]
  df_tpvt_HC[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

# pattern comparison ------------------------------------------------------

pcps01 <- read.table('/path/to/data', na.strings = c("", "NA"))
# View(pcps01)
# want the NIH uncorrected scaled score, V10
pcps <- pcps01[-c(1:2), c(5, 7, 8, 10)] # extract the ID, age, sex and score
sum(is.na(pcps$V10)) # no missing data
# [1] 0
# which(pcps$V10==999)
# integer(0)

pcps <-
  pcps %>%
  filter(V5 %in% subj_non_miss_fMRI$X0) %>% # 157 subjects
  mutate(
    V5 = as.numeric(V5), # factorize the categorical data
    V7 = as.numeric(V7)/12,
    V8 = factor(V8),
    V10 = as.numeric(V10)
  ) %>%
  `colnames<-`(c('ID', 'age', 'sex', 'score'))

# print(pcps)

pcps_fMRI_miss <- subj_non_miss_fMRI %>%
  filter(!X0 %in% pcps$ID)
print(pcps_fMRI_miss) # 14 subjects missing
write.csv(pcps_fMRI_miss, '/pcps_miss_ID.csv',row.names = F)
# save the IDs as miss_ID.csv so that our python dictionary is updated

pcps.bn.EP <- pcps %>%
  filter(ID %in% EP_demo$src_subject_id)
write.csv(pcps.bn.EP$ID, '/pcps_EP_ID.csv',row.names = F)
pcps.bn.EP <- pcps.bn.EP[,-1]

pcps.bn.HC <- pcps %>%
  filter(ID %in% HC_demo$src_subject_id)
write.csv(pcps.bn.HC$ID, '/pcps_HC_ID.csv',row.names = F)
pcps.bn.HC <- pcps.bn.HC[,-1]

## MDMR --------------------------------------------------------------------

D_pcps <-np$load('/D_pcps.npy')

# create an empty data frame of 400 by 4
df_pcps <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))

### group-wise --------------------------------------------------------------

D_pcps_EP <-np$load('/D_pcps_EP.npy')
D_pcps_HC <-np$load('/D_pcps_HC.npy')

# create an empty data frame of 400 by 4

df_pcps_EP <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))

df_pcps_HC <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))


# loop
row_index <- 1
for (i in 1:dim(D_pcps_EP)[3]) {
  set.seed(12345) # to be reproducible
  region <- D_pcps_EP[,,i]
  result <- mdmr(X = pcps.bn.EP, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_pcps_EP[row_index, 1] <- F_stats[1]
  df_pcps_EP[row_index, 3] <- F_stats[2]
  df_pcps_EP[row_index, 2] <- p_vals[1]
  df_pcps_EP[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

row_index <- 1
for (i in 1:dim(D_pcps_HC)[3]) {
  set.seed(12345)
  region <- D_pcps_HC[,,i]
  result <- mdmr(X = pcps.bn.HC, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_pcps_HC[row_index, 1] <- F_stats[1]
  df_pcps_HC[row_index, 3] <- F_stats[2]
  df_pcps_HC[row_index, 2] <- p_vals[1]
  df_pcps_HC[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

# picture sequence --------------------------------------------------------

psm01 <- read.table('/path/to/data', na.strings = c("", "NA"))
# View(psm01)

# want the NIH uncorrected scaled score, V11
psm <- psm01[-c(1:2), c(5, 7, 8, 11)] # extract the ID, age, sex and score
sum(is.na(psm$V11)) # no missing data
# [1] 0
which(psm$V11==999)
# [1] 64

psm <-
  psm[-64,] %>%
  filter(V5 %in% subj_non_miss_fMRI$X0) %>% # 157 subjects
  mutate(
    V5 = as.numeric(V5), # factorize the categorical data
    V7 = as.numeric(V7)/12,
    V8 = factor(V8),
    V11 = as.numeric(V11)
  ) %>%
  `colnames<-`(c('ID', 'age', 'sex', 'score'))

# print(psm)

psm_fMRI_miss <- subj_non_miss_fMRI %>%
  filter(!X0 %in% psm$ID)
print(psm_fMRI_miss) # 13 subjects missing
write.csv(psm_fMRI_miss, '/psm_miss_ID.csv',row.names = F)
# save the IDs as miss_ID.csv so that our python dictionary is updated

psm.bn.EP <- psm %>% # one less EP subject
  filter(ID %in% EP_demo$src_subject_id)
write.csv(psm.bn.EP$V5, '/psm_EP_ID.csv',row.names = F)
psm.bn.EP <- psm.bn.EP[,-1]

psm.bn.HC <- psm %>%
  filter(ID %in% HC_demo$src_subject_id)
write.csv(psm.bn.HC$V5, '/psm_HC_ID.csv',row.names = F)
psm.bn.HC <- psm.bn.HC[,-1]

## MDMR --------------------------------------------------------------------

D_psm <-np$load('/D_psm.npy')

# create an empty data frame of 400 by 4
df_psm <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))

### group-wise --------------------------------------------------------------

D_psm_EP <- np$load('/D_psm_EP.npy')
D_psm_HC <- np$load('/D_psm_HC.npy')

# create an empty data frame of 400 by 4

df_psm_EP <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))

df_psm_HC <- matrix(data = NA, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  `colnames<-`(c('omibus_F', 'F_pval', 'score_F', 'score_pval'))


# loop
row_index <- 1
for (i in 1:dim(D_psm_EP)[3]) {
  set.seed(12345) # to be reproducible
  region <- D_psm_EP[,,i]
  result <- mdmr(X = psm.bn.EP, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_psm_EP[row_index, 1] <- F_stats[1]
  df_psm_EP[row_index, 3] <- F_stats[2]
  df_psm_EP[row_index, 2] <- p_vals[1]
  df_psm_EP[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

row_index <- 1
for (i in 1:dim(D_psm_HC)[3]) {
  set.seed(12345)
  region <- D_psm_HC[,,i]
  result <- mdmr(X = psm.bn.HC, D = region)
  stats <- result[[1]] # returns a data frame instead of a list (called by result[1])
  p <- result[[3]]
  F_stats <- stats[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  p_vals <- p[1:2,] %>% # get the rows of interest (omnibus and score)
    t()
  df_psm_HC[row_index, 1] <- F_stats[1]
  df_psm_HC[row_index, 3] <- F_stats[2]
  df_psm_HC[row_index, 2] <- p_vals[1]
  df_psm_HC[row_index, 4] <- p_vals[2]
  row_index <- row_index + 1
}

# optional - design matrix X -----------------------------------------------

# # it must have same number of rows as Y (so n = 158, number of subjects w/ complete data)
# flanker <- flanker[,-c(1,3)]

# MDMR automatically formats the design matrix if it's not a model.matrix already, but it'd be best if we
# explicitly state the design matrix here
# given my X here is a data frame (str(flanker) -- data frame), but all columns are character strings!
# first exclude the ID column
flanker <- flanker[,-1]
colnames(flanker) <- c('age', 'sex', 'score')
flanker$age <- as.numeric(flanker$age)
flanker$age <- flanker$age/12
flanker$score <- as.numeric(flanker$score)
# create a model matrix
flanker <- model.matrix(~ score + age + sex, data = flanker)

# optional - MDMR ----------------------------------------------------------------
# import G (been converted to csv in Python)
# Gower = read.csv('/data/lavlab/students/fancat/MDMR/Gower_test.csv', header = F)
# View(Gower)
# Gower <- as.matrix(Gower)



library(tidyverse)
library(rmarkdown)
library(rlang)
library(parameters)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
library(glue)
library(ggforestplot)
library(ggbeeswarm)
library(patchwork)
library(lme4)
library(ggstance)
library(knitr)
library(grid)

setwd('/projects/home/sramesh/github/myocarditis')
# source('functions/de.R')
source('functions/masc.R')
source('functions/plot_masc.R')
source('functions/blood_abundance.R')
source('functions/blood_troponin.R')

### Read in obs matrix and cluster annotations
obs <- read_csv(glue('/projects/home/sramesh/myo_final/blood/final/myo_blood_global_obs.csv'))
cluster_annots <- readxl::read_excel('/projects/home/sramesh/myo_final/blood/other_stuff/cluster_annotations.xlsx')

## Figure 2C
# run masc for case vs control by cluster (or read in file if it exists)
if (!file.exists('output/masc_for_deg_case_control_by_cluster.csv')) {
  cluster_masc_res <- run_masc(obs,
                               c('deg_case_control', 'control', 'case'),
                               'lineage_cluster',
                               fixed_effects = c('sex', 'ici_type'))
  write_csv(cluster_masc_res, 'output/masc_for_deg_case_control_by_cluster.csv')
} else {
  cluster_masc_res <- read_csv('output/masc_for_deg_case_control_by_cluster.csv')
}

# rename clusters
cluster_name_dict <- setNames(cluster_annots$cluster_name_w_num,
                              cluster_annots$lineage_cluster)
cluster_masc_res$cluster_name_w_num <- cluster_masc_res$cluster %>%
  str_remove('cluster') %>%
  recode(!!!cluster_name_dict)

# filter for only cd8 and nk clusters and plot
cd8_clusters <- obs %>%
  filter((lineage == 'CD8 and NK') & (lineage != 'Doublets and RBCs')) %>%
  pull(cluster_name_w_num) %>%
  unique()
masc_helper(obs,
            c('deg_case_control', 'control', 'case'),
            'cluster_name_w_num',
            fixed_effects = c('sex', 'ici_type'),
            group_subset=cd8_clusters,
            masc_res=cluster_masc_res,
            colors=c('tomato4', 'slategray'))

## Figure 2F
cd4_clusters <- obs %>%
  filter((lineage == 'CD4') & (lineage != 'Doublets and RBCs')) %>%
  pull(cluster_name_w_num) %>%
  unique()
masc_helper(obs,
            c('deg_case_control', 'control', 'case'),
            'cluster_name_w_num',
            fixed_effects = c('sex', 'ici_type'),
            group_subset=cd4_clusters,
            masc_res=cluster_masc_res,
            colors=c('tomato4', 'slategray'))

## Figure 2I
mnp_clusters <- c("11. b-DC1: CLEC9A, IDO1", "6. b-DC2: CD1C, CLEC10A", "8. b-DC3: RNASE2, F13A1",
                  "10. b-pDC: LILRA4, JCHAIN", "7. b-MNP: IFI44L, IFI6", "5. b-MNP: C1QB, CTSL",
                  "9. b-MNP: CCL3, IL1B", "1. b-MNP: CD14, CSF3R", "2. b-MNP: S100A12, CTSD",
                  "4. b-MNP: DUSP6, MARCKS", "3. b-MNP: FCGR3A, CDKN1C", "12. b-Mast: PRSS57, CD34")
mnp_clusters <- mnp_clusters[length(mnp_clusters):1]
mnp_cluster_masc_res <- cluster_masc_res %>%
  filter(cluster_name_w_num %in% mnp_clusters) %>%
  # remove results for mast cells
  mutate(across(c(model.pvalue, OR, OR.95pct.ci.lower, OR.95pct.ci.upper),
                ~ ifelse(cluster_name_w_num == '12. b-Mast: PRSS57, CD34', NA, .)),
  # make cluster names as factor so they are ordered correctly
         cluster_name_w_num = factor(cluster_name_w_num, levels=mnp_clusters)) %>%
  arrange(cluster_name_w_num)
masc_helper(obs,
            c('deg_case_control', 'control', 'case'),
            'cluster_name_w_num',
            fixed_effects = c('sex', 'ici_type'),
            group_subset=mnp_clusters,
            masc_res=mnp_cluster_masc_res,
            colors=c('tomato4', 'slategray'),
            row_order=T)

## Figure 2J
# read in troponin data and format obs data for troponin analysis
troponin_data <- read_csv('/projects/home/sramesh/myo_final/blood/other_stuff/blood_troponin.csv')
filt_obs <- obs %>%
  filter(steroid_treatment == 'pre_steroid') %>% # only one sample per donor bc only looking at pre-steroid samples
  filter(!str_detect(lineage, "Doublet")) %>%
  left_join(troponin_data, by='sample_id') %>%
  filter(abs(days_sample_to_troponin) <= 3) %>%
  rename(nearest_troponin = troponin_nearest_sample,
         cluster_names = lineage_cluster,
         lineage_names = lineage,
         donor = sample_id) %>%
  select(donor, nearest_troponin, cluster_names, lineage_names)

# select significantly abundant lineages
lineage_masc_res <- read_csv('output/masc_for_deg_case_control_by_lineage.csv')
sig_lineages <- lineage_masc_res %>%
  filter((model.pvalue < 0.1)) %>%
  pull(lineage)

# calculate troponin trends
troponin_percents <- troponin_get_percent_per_level(filt_obs, level='lineage')
troponin_model <- troponin_fit_model(troponin_percents, level='lineage')
troponin_model <- troponin_model %>%
  mutate(trop_coef = unlist(trop_coef),
         trop_se = unlist(trop_se),
         trop_pval = unlist(trop_pval)) %>%
  select(-c(data, model))
troponin_percents <- troponin_percents %>% filter(lineage_names %in% sig_lineages)
troponin_model <- troponin_model %>% filter(lineage_names %in% sig_lineages)
troponin_model$padj <- p.adjust(troponin_model$trop_pval, method='fdr')

write_csv(troponin_model, 'output/troponin_by_sig_lineage.csv')
troponin_plot_model(troponin_model, troponin_percents, 'Significant Lineages',
               level='lineage', point_size = 2.2, type='detailed')

## Figure 2K

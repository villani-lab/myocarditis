Supplemental Figure 1
================

## Set up

Load R libraries

``` r
library(reticulate)
library(tidyverse)
library(rmarkdown)
library(rlang)
library(parameters)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
library(glue)
library(fgsea)
library(ggforestplot)
library(ggbeeswarm)
library(patchwork)
library(lme4)
library(ggstance)
library(knitr)
library(grid)
library(readxl)

setwd('/projects/home/ikernin/github_code/myocarditis/functions')
source('de.R')
source('gsea.R')
source('masc.R')
source('plot_masc.R')
source('blood_abundance.R')
source('blood_troponin.R')

use_python("/projects/home/nealpsmith/.conda/envs/updated_pegasus/bin/python")
wd = '/projects/home/sramesh/github/myocarditis'
```

Load Python packages

``` python
import pandas as pd
import pegasus as pg
import os
import sys
import warnings

sys.path.append("/projects/home/ikernin/github_code/myocarditis/functions")
import python_functions

sys.path.append("/projects/home/sramesh/github/myocarditis/functions")
import python_functions as pyfun

wd = '/projects/home/sramesh/github/myocarditis'
warnings.filterwarnings('ignore')
```

Read in single-cell data and cluster names

``` python
# read tissue data
tissue_global = pg.read_input('/projects/home/ikernin/projects/myocarditis/updated_datasets/tissue_all_with_dcs.zarr')

# read blood data
```

``` python
adata = pg.read_input('/projects/home/sramesh/myo_final/blood/final/myo_blood_global.zarr')
```

``` python
cluster_annots = pd.read_excel('/projects/home/sramesh/myo_final/blood/other_stuff/cluster_annotations.xlsx')
name_dict = dict(zip(cluster_annots['lineage_cluster_mod'], cluster_annots['cluster_name_w_num']))
```

``` r
# read tissue data
tissue_global_obs <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/metadata/tissue_full_obs.csv')

# read blood data
obs <- read_csv('/projects/home/sramesh/myo_final/blood/final/myo_blood_global_obs.csv')
cluster_annots <- readxl::read_excel('/projects/home/sramesh/myo_final/blood/other_stuff/cluster_annotations.xlsx')
```

## Supplemental Figure 1A

``` python
python_functions.make_gene_dotplot(tissue_global.to_anndata(),
             cluster_order=['Endothelial cells', 'cDC', 'pDC',
                            'Myeloid cells', 'B and plasma cells',
                            'T and NK cells', 'Cardiomyocytes', 'Mural cells',
                            'Fibroblasts',
                            'Neuronal cells'
                            ],
             gene_order=['VWF', 'AQP1', 'CA4', 'HEY1', 'HLA-DQA1', 'CD1C',
                         'CLEC9A', 'LILRA4', 'IL3RA', 'CD14', 'CD68', 'C1QA', 'FCGR3A',
                         'CD79A', 'MZB1', 'MS4A1', 'CD3D', 'CD8A', 'KLRB1', 'TNNI3',
                         'MB', 'RGS5', 'MYH11', 'KCNJ8', 'ACTA2', 'DCN',
                         'PDGFRA', 'PLP1', 'NRXN1'],
             title='All Cells'
             )
```

<img src="supplemental_figure_1_files/figure-gfm/supp_1a-1.png" width="1152" />

## Supplemental Figure 1B

``` r
p_1b <- ggplot(tissue_global_obs, aes(y = umap_name, fill = institution)) +
  geom_bar(position = 'fill') +
  scale_fill_manual(values = c('#0273b5', '#e18728')) +
  scale_x_continuous(labels = scales::percent,
                     expand = c(0,0),
                     limits = c(0,1)) +
  theme_classic()
print(p_1b)
```

![](supplemental_figure_1_files/figure-gfm/supp_1b-3.png)<!-- -->

## Supplemental Figure 1C

``` r
melanoma <- c("SIC_171",
              "SIC_182",
              "SIC_237",
              "SIC_319")

renal <- c("SIC_197",
           "SIC_199",
           "SIC_217",
           "SIC_232",
           "SIC_317")

other <- c("SIC_48",
           "SIC_153",
           "SIC_164",
           "SIC_175",
           "SIC_177",
           "SIC_258",
           "SIC_264",
           "SIC_333")
non_sanger <- unique(tissue_global_obs$donor)[!str_detect(unique(tissue_global_obs$donor), 'Sanger')]

p_1c <- tissue_global_obs %>%
  filter(donor %in% non_sanger) %>%
  filter(condition == 'myocarditis') %>%
  filter(on_steroids == 'False') %>%
  mutate(tumor_type = case_when(
    donor %in% melanoma ~ 'Melanoma',
    donor %in% renal ~ 'Renal Cell',
    donor %in% other ~ 'Other',
    TRUE ~ "N/A"
  )) %>%
  ggplot(aes(y = donor, fill = umap_name)) +
  geom_bar(position = 'fill') +
  scale_fill_manual(values = c(
    "#b3e900", # b and plasma
    "#ff7f00", # cardiomyocytes
    "#7f80cd", # cDC
    "#f781bf", # doublets
    "#ff0029", # endothelial
    "#00d2d5", # fibroblast
    "#377eb8", # mural
    "#984ea3", # mnp
    "#af8d00", # neuronal
    "#c42e60", # pdf
    '#66a61e' # t and nk

  )) +
  scale_x_continuous(labels = scales::percent,
                     expand = c(0,0),
                     limits = c(0,1)) +
  scale_y_discrete(position = "right") +
  facet_grid(tumor_type ~ ., scales = 'free_y', space = 'free', switch = 'y') +
  theme_classic() +
  labs(
    x = '% of donor cells',
    y = element_blank(),
    fill = 'Lineage'
  )
print(p_1c)
```

![](supplemental_figure_1_files/figure-gfm/supp_1c-1.png)<!-- -->

## Supplemental Figure 1D

``` r
abund_data <- read_excel("/projects/home/nealpsmith/projects/myocarditis/data/updated_supplemental_tables/abundance_analysis_aggregated.xlsx",
                         sheet = "masc_res")

tissue_data <- abund_data %>%
  dplyr::filter(tissue == "heart")

withmets <- tissue_data %>%
  dplyr::filter(contrast == "irMyocarditis vs. control") %>%
  dplyr::select(cluster, OR) %>%
  `colnames<-`(c("cluster", "OR_mets"))

nomets <- tissue_data %>%
  dplyr::filter(contrast == "irMyocarditis vs. control without heart metastasis") %>%
  dplyr::select(cluster, OR) %>%
  `colnames<-`(c("cluster", "OR_nomets"))

plot_df <- left_join(withmets, nomets, by = "cluster") %>%
  dplyr::filter(cluster != "h-pDC: LILRA4 IRF8")

ggplot(plot_df, aes(x = OR_mets, y = OR_nomets)) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline() +
  geom_vline(xintercept = 1) +
  geom_hline(yintercept = 1) +
  ylab("OR : without cardiac metastases") +
  xlab("OR : all samples") +
  theme_classic(base_size = 20)
```

![](supplemental_figure_1_files/figure-gfm/supp_1d-1.png)<!-- -->

## Supplemental Figure 1E

``` r
all_data <- read.csv("/projects/home/nealpsmith/projects/myocarditis/data/supplemental_tables/all_degs_concatenated.csv")
tcell_no_mm <- read.csv("/projects/home/nealpsmith/projects/myocarditis/data/t_no_micromets_de_by_condition_all_results.csv")
tcell_no_mm$X <- NULL
tcell_no_mm$cluster <- NULL
colnames(tcell_no_mm) <- c(paste("no_mm", colnames(tcell_no_mm)[1:6], sep = "_"), "gene_symbol", "cluster")


### MNP ###
mnp_no_mm <- read.csv("/projects/home/nealpsmith/projects/myocarditis/data/myeloid_no_micromets_de_by_condition_all_results.csv")
mnp_no_mm$X <- NULL
mnp_no_mm$cluster <- NULL
colnames(mnp_no_mm) <- c(paste("no_mm", colnames(mnp_no_mm)[1:6], sep = "_"), "gene_symbol", "cluster")

### Non-immune ###
nonimmune_no_mm <- read.csv("/projects/home/nealpsmith/projects/myocarditis/data/nonimmune_no_micromets_de_by_condition_all_results.csv")
nonimmune_no_mm$X <- NULL
nonimmune_no_mm$cluster <- NULL
colnames(nonimmune_no_mm) <- c(paste("no_mm", colnames(nonimmune_no_mm)[1:6], sep = "_"), "gene_symbol", "cluster")

### How about one for all DEGs ###
no_mm_data_all <- rbind(tcell_no_mm, mnp_no_mm, nonimmune_no_mm)

plot_degs_all <- all_data %>%
  dplyr::filter(tissue == "heart", padj < 0.1) %>%
  dplyr::left_join(no_mm_data_all, by = c("gene_symbol", "cluster"))

ggplot(plot_degs_all, aes(x = log2FoldChange, y = no_mm_log2FoldChange)) +
  # geom_point(data = plot_degs_all %>% dplyr::filter(label == ""), color = "grey") +
  geom_point(data = plot_degs_all, fill = "black", pch = 21, size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline() +
  ylab("log2FC : without cardiac metastases") +
  xlab("log2FC : with cardiac metastases") +
  # geom_text_repel(aes(label = label), min.segment.length = 0.1, max.overlaps = 100) +
  theme_classic(base_size = 20)
```

    ## Warning: Removed 154 rows containing missing values (`geom_point()`).

![](supplemental_figure_1_files/figure-gfm/supp_1e-1.png)<!-- -->

## Supplemental Figure 1F

``` r
# read in serum values per sample
serum_df <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/metadata/tidy_serum.csv')

# get between condition statistics for log1p values
serum_model <- serum_df %>%
  select(!c(patient, timepoint, Sample)) %>%
  pivot_longer(cols=!c(condition),
               names_to='protein',
               values_to='value') %>%
  mutate(log1p_value = log1p(value)) %>%
  group_by(condition, protein) %>%
  summarize(log1p_values = list(log1p_value)) %>%
  pivot_wider(names_from = condition,
              values_from = log1p_values) %>%
  group_by(protein) %>%
  summarize(mean_case = mean(unlist(myocarditis)),
            std_case = sd(unlist(myocarditis)),
            n_case = length(unlist(myocarditis)),
            mean_control = mean(unlist(control)),
            std_control = sd(unlist(control)),
            n_control = length(unlist(control)),
            p_value = t.test(unlist(control), unlist(myocarditis))$p.value,
            t_stat = t.test(unlist(control), unlist(myocarditis))$statistic,
            log2FC = log2(mean_case/mean_control))
kable(serum_model)

# proteins to annotate in figures
serum_annot <-  c('MIG/CXCL9', 'IP-10', 'TNFα','IFNγ', 'IL-2',
                  'IL-12p40', 'IL-15', 'IL-18', 'IL-27', '6CKine')

micromet_donors <- c("SIC_136", "SIC_171", "SIC_232", "SIC_182")

# get data for ext fig 1e
ext_1e_df <- serum_df %>%
  select(!c(Sample, timepoint)) %>%
  pivot_longer(cols = !c(patient, condition),
               names_to='protein') %>%
  filter(protein %in% serum_annot) %>%
  left_join(serum_model %>% select(protein, p_value)) %>%
  mutate(
    protein = case_when(
      protein == '6CKine' ~ 'CCL21',
      protein == 'IP-10' ~ 'CXCL10',
      protein == 'MIG/CXCL9' ~ 'CXCL9',
      TRUE ~ protein),
    condition2 =
      case_when(
        patient %in% micromet_donors ~ 'micromet',
        TRUE ~ condition)
    )

ext_1e_labels <- ext_1e_df %>%
  select(protein, p_value) %>%
  mutate(facet_name = str_c(protein, str_c("p=", signif(p_value, 2)), sep = '\n')) %>%
  select(protein, facet_name) %>%
  distinct() %>%
  deframe()


p_1e <- ggplot(ext_1e_df, aes(x = condition, y = value + 1, fill = condition)) +
  geom_boxplot(alpha = 0.6,
               outlier.shape = NA) +
  geom_jitter(aes(fill = condition2,
                  size = condition2 == 'micromet'),
              shape=21) +
  scale_y_log10()+
  scale_color_manual(values = c('slategray', 'tomato4')) +
  scale_fill_manual(values = c('slategray', 'yellow', 'tomato4')) +
  scale_size_manual(values = c(1,3)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  guides(size = "none") +
  facet_wrap(~protein,
             labeller = labeller(protein = ext_1e_labels),
             scales = 'free_y',
             nrow=2) +
  labs(fill = 'Condition',
       color = 'Condition',
       x = element_blank(),
       y = "pg/mL + 1")
print(p_1e)
```

![](supplemental_figure_1_files/figure-gfm/supp_1f-1.png)<!-- -->

| protein      | mean\_case | std\_case | n\_case | mean\_control | std\_control | n\_control |  p\_value |     t\_stat |      log2FC |
| :----------- | ---------: | --------: | ------: | ------------: | -----------: | ---------: | --------: | ----------: | ----------: |
| 6CKine       |  6.4363743 | 0.7101943 |      16 |     4.5567455 |    2.5437937 |         10 | 0.0459437 | \-2.2817143 |   0.4982444 |
| BCA-1        |  4.8896005 | 1.0394112 |      16 |     3.8822775 |    0.5560748 |         10 | 0.0037935 | \-3.2104880 |   0.3328134 |
| CTACK        |  7.3854061 | 0.3856017 |      16 |     7.0697631 |    0.6390623 |         10 | 0.1818264 | \-1.4097232 |   0.0630154 |
| EGF          |  5.2334156 | 0.7222864 |      16 |     4.7273420 |    1.7544192 |         10 | 0.4043610 | \-0.8673928 |   0.1467236 |
| ENA-78       |  6.4743338 | 0.6180568 |      16 |     6.7690392 |    0.8540781 |         10 | 0.3586866 |   0.9471231 | \-0.0642193 |
| Eotaxin      |  3.9514577 | 0.6351661 |      16 |     4.0296076 |    0.9319650 |         10 | 0.8187338 |   0.2334441 | \-0.0282544 |
| Eotaxin-2    |  6.3029821 | 0.8643721 |      16 |     6.0331097 |    0.2604114 |         10 | 0.2576432 | \-1.1670045 |   0.0631328 |
| Eotaxin-3    |  1.2128693 | 1.4846713 |      16 |     0.6737767 |    1.3778945 |         10 | 0.3573158 | \-0.9418352 |   0.8480817 |
| FGF-2        |  5.0137612 | 1.2062506 |      16 |     3.8426111 |    1.4997129 |         10 | 0.0533933 | \-2.0838641 |   0.3838063 |
| FLT-3L       |  3.8099334 | 0.7345215 |      16 |     3.3411918 |    1.1989049 |         10 | 0.2855658 | \-1.1127215 |   0.1894030 |
| Fractalkine  |  5.3353875 | 0.7747217 |      16 |     4.3582980 |    1.5767914 |         10 | 0.0932337 | \-1.8266092 |   0.2918282 |
| G-CSF        |  2.4111508 | 2.2985718 |      16 |     1.9313054 |    1.9185079 |         10 | 0.5716786 | \-0.5742302 |   0.3201456 |
| GM-CSF       |  4.1298428 | 1.4318275 |      16 |     3.5633644 |    1.3597476 |         10 | 0.3233776 | \-1.0124982 |   0.2128469 |
| GROα         |  3.4976916 | 0.7170439 |      16 |     3.0594765 |    1.2331317 |         10 | 0.3260535 | \-1.0210509 |   0.1931183 |
| I-309        |  1.7733368 | 0.5745722 |      16 |     1.2702573 |    0.4969847 |         10 | 0.0276843 | \-2.3628207 |   0.4813458 |
| IFN-α2       |  4.3034714 | 1.9341091 |      16 |     3.0050394 |    2.2079845 |         10 | 0.1443742 | \-1.5288179 |   0.5181170 |
| IFNγ         |  1.9435044 | 1.4059266 |      16 |     0.9612081 |    0.6652875 |         10 | 0.0250747 | \-2.3979847 |   1.0157397 |
| IL-10        |  2.7738412 | 1.7354165 |      16 |     1.2456889 |    1.3806524 |         10 | 0.0209809 | \-2.4827516 |   1.1549414 |
| IL-12p40     |  6.1262569 | 1.3028557 |      16 |     4.4349363 |    1.2923763 |         10 | 0.0042644 | \-3.2363434 |   0.4660925 |
| IL-12p70     |  2.3888904 | 2.1255133 |      16 |     1.2592322 |    0.9518119 |         10 | 0.0776223 | \-1.8497685 |   0.9237963 |
| IL-13        |  4.4144231 | 1.0829592 |      16 |     3.6266109 |    1.5080893 |         10 | 0.1716196 | \-1.4365869 |   0.2836030 |
| IL-15        |  3.7151701 | 1.0837076 |      16 |     2.2245260 |    1.1138127 |         10 | 0.0033590 | \-3.3545572 |   0.7399303 |
| IL-16        |  1.2938391 | 2.4059516 |      16 |     1.3318474 |    2.8292304 |         10 | 0.9722888 |   0.0352558 | \-0.0417706 |
| IL-17A       |  2.8404342 | 1.7656720 |      16 |     1.5745222 |    1.5799273 |         10 | 0.0714845 | \-1.8988200 |   0.8511974 |
| IL-17E/IL-25 |  5.7418729 | 2.6941808 |      16 |     3.9094564 |    2.8434321 |         10 | 0.1198353 | \-1.6310400 |   0.5545534 |
| IL-17F       |  2.9082138 | 0.9201737 |      16 |     2.3049660 |    1.6178459 |         10 | 0.3022050 | \-1.0754080 |   0.3353879 |
| IL-18        |  4.5529270 | 0.8469335 |      16 |     3.1609961 |    1.7582195 |         10 | 0.0379982 | \-2.3395791 |   0.5264150 |
| IL-1RA       |  3.9844177 | 1.3788921 |      16 |     2.9077550 |    1.0645687 |         10 | 0.0355673 | \-2.2345065 |   0.4544632 |
| IL-1α        |  3.7092153 | 1.7665737 |      16 |     2.7010582 |    1.5049316 |         10 | 0.1350191 | \-1.5527928 |   0.4575893 |
| IL-1β        |  3.7776029 | 1.6315738 |      16 |     2.5616238 |    1.4620230 |         10 | 0.0619753 | \-1.9722459 |   0.5604124 |
| IL-2         |  1.8627309 | 1.6871037 |      16 |     0.6756830 |    0.7037706 |         10 | 0.0209716 | \-2.4891432 |   1.4630008 |
| IL-20        |  2.5522352 | 3.1023734 |      16 |     1.2752235 |    2.7138124 |         10 | 0.2819649 | \-1.1039848 |   1.0010112 |
| IL-21        |  1.4868531 | 1.6123368 |      16 |     1.0480728 |    1.8236082 |         10 | 0.5409368 | \-0.6236353 |   0.5045232 |
| IL-22        |  3.0965614 | 2.0476407 |      16 |     2.7713745 |    2.0495117 |         10 | 0.6981034 | \-0.3937390 |   0.1600654 |
| IL-23        |  2.2091780 | 3.4142322 |      16 |     2.4069013 |    4.1287228 |         10 | 0.9006628 |   0.1267562 | \-0.1236673 |
| IL-27        |  8.1959283 | 0.4365125 |      16 |     7.1586560 |    0.8112040 |         10 | 0.0028027 | \-3.7208548 |   0.1952186 |
| IL-28A       |  1.9237581 | 2.8072062 |      16 |     1.3484567 |    2.8816697 |         10 | 0.6227331 | \-0.5001813 |   0.5126182 |
| IL-3         |  0.3936793 | 0.2708367 |      16 |     0.3364586 |    0.2951412 |         10 | 0.6257284 | \-0.4962530 |   0.2265919 |
| IL-33        |  1.5510885 | 1.9664495 |      16 |     2.0807321 |    2.6816080 |         10 | 0.5968677 |   0.5403444 | \-0.4238102 |
| IL-4         |  1.2619752 | 1.1323697 |      16 |     0.9644388 |    0.8338522 |         10 | 0.4496054 | \-0.7690763 |   0.3879219 |
| IL-5         |  2.6155659 | 0.9475224 |      16 |     1.9732390 |    1.1590662 |         10 | 0.1600244 | \-1.4718343 |   0.4065574 |
| IL-6         |  2.4841638 | 0.9811992 |      16 |     2.0299885 |    1.7216854 |         10 | 0.4607923 | \-0.7605674 |   0.2912888 |
| IL-7         |  2.4569886 | 1.0001989 |      16 |     1.8172582 |    0.7737255 |         10 | 0.0806058 | \-1.8286225 |   0.4351277 |
| IL-8         |  3.1816575 | 1.3405947 |      16 |     2.8762683 |    1.3367150 |         10 | 0.5778419 | \-0.5661131 |   0.1455803 |
| IL-9         |  2.3737289 | 1.3995452 |      16 |     0.8601500 |    1.3993118 |         10 | 0.0146041 | \-2.6830949 |   1.4644951 |
| IP-10        |  7.3232318 | 1.7088805 |      16 |     5.2281214 |    1.9587047 |         10 | 0.0125870 | \-2.7844134 |   0.4861878 |
| LIF          |  1.7842529 | 1.5358532 |      16 |     1.1075488 |    1.7872170 |         10 | 0.3358263 | \-0.9904071 |   0.6879499 |
| M-CSF        |  5.4138877 | 0.8799103 |      16 |     4.4735527 |    1.8044625 |         10 | 0.1506946 | \-1.5376150 |   0.2752439 |
| MCP-1        |  5.6576888 | 0.8842821 |      16 |     5.4289545 |    1.5628755 |         10 | 0.6797663 | \-0.4224752 |   0.0595385 |
| MCP-2        |  3.9274040 | 0.3776841 |      16 |     3.9516614 |    0.1769999 |         10 | 0.8270705 |   0.2209949 | \-0.0088833 |
| MCP-3        |  3.2201558 | 1.0213659 |      16 |     2.5523514 |    1.0656746 |         10 | 0.1310709 | \-1.5794558 |   0.3353036 |
| MCP-4        |  4.0842874 | 1.7666237 |      16 |     4.2803509 |    1.5920842 |         10 | 0.7726134 |   0.2927517 | \-0.0676447 |
| MDC          |  6.5570121 | 0.3519171 |      16 |     5.9321419 |    2.1367101 |         10 | 0.3822666 | \-0.9170513 |   0.1444855 |
| MIG/CXCL9    |  9.1631650 | 0.6288525 |      16 |     8.2298567 |    1.1837104 |         10 | 0.0399116 | \-2.2988098 |   0.1549787 |
| MIP-1α       |  4.3172488 | 1.1287316 |      16 |     3.2865776 |    1.2553439 |         10 | 0.0487888 | \-2.1161581 |   0.3935262 |
| MIP-1β       |  4.3389583 | 0.5230917 |      16 |     3.5320302 |    1.3570783 |         10 | 0.1003234 | \-1.7986546 |   0.2968510 |
| MIP-1δ       |  8.6986187 | 0.6016228 |      16 |     8.0169640 |    1.1240404 |         10 | 0.1022145 | \-1.7661075 |   0.1177303 |
| PDGF-AA      |  7.2625498 | 0.4149287 |      16 |     6.8714718 |    2.4434248 |         10 | 0.6275551 | \-0.5016323 |   0.0798570 |
| PDGF-AB/BB   |  9.8534586 | 0.2195944 |      16 |     9.5292300 |    1.2586629 |         10 | 0.4397492 | \-0.8069556 |   0.0482706 |
| RANTES       |  7.4055882 | 0.4235097 |      16 |     6.8356067 |    2.4265108 |         10 | 0.4798918 | \-0.7358396 |   0.1155449 |
| SCF          |  3.2519765 | 2.0302255 |      16 |     1.6751457 |    2.4194080 |         10 | 0.1043613 | \-1.7174330 |   0.9570302 |
| SDF-1α+β     |  8.8320183 | 0.5500243 |      16 |     8.6591937 |    0.5690661 |         10 | 0.4549002 | \-0.7631023 |   0.0285105 |
| TARC         |  4.5425314 | 0.5415717 |      16 |     4.6026783 |    0.9859257 |         10 | 0.8623928 |   0.1769515 | \-0.0189771 |
| TGFα         |  2.4649651 | 0.8532875 |      16 |     2.2624766 |    1.3603165 |         10 | 0.6798823 | \-0.4217117 |   0.1236643 |
| TNFα         |  5.0640180 | 0.7379993 |      16 |     4.2743596 |    0.9044013 |         10 | 0.0335644 | \-2.3201738 |   0.2445742 |
| TNFβ         |  2.8253295 | 0.9145280 |      16 |     2.4683577 |    0.9559061 |         10 | 0.3583456 | \-0.9418542 |   0.1948676 |
| TPO          |  5.6742710 | 1.0042720 |      16 |     5.4147425 |    2.8027303 |         10 | 0.7836404 | \-0.2817361 |   0.0675423 |
| TRAIL        |  3.7510727 | 0.4668613 |      16 |     3.7174592 |    0.5276374 |         10 | 0.8707802 | \-0.1650769 |   0.0129863 |
| TSLP         |  0.7325628 | 1.1914176 |      16 |     1.4436622 |    2.9711651 |         10 | 0.4859080 |   0.7214547 | \-0.9787090 |
| VEGF-A       |  5.2810634 | 0.6310147 |      16 |     5.1891024 |    1.8756073 |         10 | 0.8837901 | \-0.1498372 |   0.0253435 |
| sCD40L       |  8.4292255 | 0.4433345 |      16 |     8.1623234 |    1.8579870 |         10 | 0.6651690 | \-0.4463922 |   0.0464202 |

## Supplemental Figure 1G

``` r
serum_data <- read.csv("/projects/home/nealpsmith/projects/myocarditis/data/updated_supplemental_tables/all_serum_stats.csv")

withmets <- serum_data %>%
  dplyr::filter(contrast == "irMyocarditis vs. control with mets") %>%
  dplyr::select(protein, t_statistic) %>%
  `colnames<-`(c("protein", "t_stat_mets"))

nomets <- serum_data %>%
  dplyr::filter(contrast == "irMyocarditis vs. control without mets") %>%
  dplyr::select(protein, t_statistic) %>%
  `colnames<-`(c("protein", "t_stat_nomets"))

plot_df <- left_join(withmets, nomets, by = "protein")

ggplot(plot_df, aes(x = -t_stat_mets, y = -t_stat_nomets)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline() +
  ylab("T-statistic : without cardiac metastases") +
  xlab("T-statistic : all samples") +
  theme_classic(base_size = 20)
```

![](supplemental_figure_1_files/figure-gfm/supp_1g-1.png)<!-- -->

## Supplemental Figure 1h

``` python
cluster_order = ['pDCs', 'B and plasma', 'cDCs', 'CD4', 'MNP', 'CD8 and NK']
gene_order = ['LILRA4', 'IL3RA', 'CD79A', 'MS4A1', 'JCHAIN', 'CD1C', 'CLEC10A', 'CLEC9A', 'CD4', 'CD3D', 'IL7R', 'LYZ',
              'CD14', 'AIF1', 'CD8A', 'NKG7', 'KLRF1']
pyfun.make_gene_dotplot(adata.to_anndata(), cluster_order, gene_order, 'All cells', figsize=(5, 4), names='lineage')
```

<img src="supplemental_figure_1_files/figure-gfm/supp_1h-1.png" width="1152" />

## Supplemental Figure 1I

``` r
# run masc for pre- vs post-steroid by lineage (or read in file if it exists)
if (!file.exists(glue('{wd}/output/masc_for_steroid_treatment_by_lineage.csv'))) {
  lineage_masc_res <- run_masc(obs,
                               c('steroid_treatment', 'post_steroid', 'pre_steroid'),
                               'lineage')
  write_csv(lineage_masc_res, glue('{wd}/output/masc_for_steroid_treatment_by_lineage.csv'))
} else {
  lineage_masc_res <- read_csv(glue('{wd}/output/masc_for_steroid_treatment_by_lineage.csv'))
}

lineage_masc_res$lineage <- str_replace_all(lineage_masc_res$lineage, '_', ' ')

masc_helper(obs,
            c('steroid_treatment', 'post_steroid', 'pre_steroid'),
            'lineage',
            masc_res = lineage_masc_res,
            colors = c('darkolivegreen4', 'darkolivegreen2'))
```

![](supplemental_figure_1_files/figure-gfm/supp_1I-3.png)<!-- -->

## Supplemental Figure 1J

``` r
#### first run de
if (!file.exists(glue('{wd}/output/lineage_de_by_steroid_treatment_all_results.csv'))) {
  counts <- read_counts(glue("{wd}/output/pb_counts_by_sample_id_and_lineage.csv"))
  meta <- read_meta(glue("{wd}/output/pb_meta_by_sample_id_and_lineage.csv"))
  meta <- meta %>%
    filter(steroid_treatment != "NA") %>%
    filter(!str_detect(lineage, "Doublet"))
  steroid_contrast_vec <- c('steroid_treatment', 'pre_steroid', 'post_steroid')
  steroid_lineage_results <- run_de_by_comp_var(counts, meta, glue('{wd}/output/lineage'), steroid_contrast_vec,
                                                deseq_formula = formula("~ steroid_treatment"))
} else {
  steroid_lineage_results <- read_csv(glue('{wd}/output/lineage_de_by_steroid_treatment_all_results.csv'))
}

#### now plot barplot
colors <- c('darkolivegreen2', 'darkolivegreen4')
cluster_order <- c('pDCs', 'B and plasma', 'cDCs', 'CD4', 'MNP', 'CD8 and NK')
res <- steroid_lineage_results %>%
  mutate(cluster = str_replace_all(cluster, '_', ' '))
res <- res[res$padj < 0.1, ]
ups <- downs <- numeric(length(cluster_order))

for (i in seq_along(cluster_order)) {
  clust <- cluster_order[i]
  up_genes <- res[res$cluster == clust & res$log2FoldChange > 0, ]
  down_genes <- res[res$cluster == clust & res$log2FoldChange < 0, ]
  ups[i] <- nrow(up_genes)
  downs[i] <- -nrow(down_genes)
}

df <- data.frame(cluster = cluster_order, up = ups, down = downs)
df <- df %>% mutate(cluster = factor(cluster, levels = rev(cluster_order)))

ggplot(df, aes(x = up, y = cluster)) +
  geom_bar(stat = "identity", fill = colors[1]) +
  geom_bar(stat = "identity", aes(x = down, y = cluster), fill = colors[2]) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        axis.text = element_text(color = 'black')) +
  labs(x = element_blank(),
       y = element_blank())
```

![](supplemental_figure_1_files/figure-gfm/supp_1j-1.png)<!-- -->

## Supplemental Figure 1k

``` python
bcells = pg.read_input('/projects/home/sramesh/myo_final/blood/final/myo_blood_bcells.h5ad')
bcells.obs['cluster_name_w_num'] = adata.obs['cluster_name_w_num']
bcells.obs[['umap_number', 'umap_name']] = bcells.obs['cluster_name_w_num'].str.split('\. ', expand=True)
bcells.obs['umap_number'] = bcells.obs['umap_number'].astype(int).astype('category')
bcells.obs['umap_name'] = bcells.obs['umap_name'].astype('category')
pyfun.plot_umap(bcells, 'Blood: B and Plasma Cells', pyfun.blood_b_pal, marker_multiplier=13)
```

    ## 2024-02-14 16:38:27,500 - pegasusio.readwrite - INFO - h5ad file '/projects/home/sramesh/myo_final/blood/final/myo_blood_bcells.h5ad' is loaded.
    ## 2024-02-14 16:38:27,500 - pegasusio.readwrite - INFO - Function 'read_input' finished in 2.08s.

<img src="supplemental_figure_1_files/figure-gfm/supp_1k-1.png" width="960" />

## Supplemental Figure 1L

``` python
cluster_order = [f'bcells_{i}' for i in [4, 5, 3, 2, 6, 1]]
cluster_order = [name_dict[i] for i in cluster_order]
gene_order = ['FOXP1', 'LINC01857', 'BACH2', 'CIB1', 'FGR', 'AIM2', 'COCH', 'XBP1', 'JCHAIN', 'TCL1A', 'CXCR4']
pyfun.make_gene_dotplot(bcells.to_anndata(), cluster_order, gene_order, 'B and plasma', figsize=(5, 4),
                        names='cluster_name_w_num')
```

<img src="supplemental_figure_1_files/figure-gfm/supp_1l-3.png" width="1152" />

## Supplemental Figure 1M

``` r
# run masc for case vs control by cluster (or read in file if it exists)
if (!file.exists(glue('{wd}/output/masc_for_deg_case_control_by_cluster.csv'))) {
  cluster_masc_res <- run_masc(obs,
                               c('deg_case_control', 'control', 'case'),
                               'lineage_cluster',
                               fixed_effects = c('sex', 'ici_type'))
  write_csv(cluster_masc_res, glue('{wd}/output/masc_for_deg_case_control_by_cluster.csv'))
} else {
  cluster_masc_res <- read_csv(glue('{wd}/output/masc_for_deg_case_control_by_cluster.csv'))
}

# rename clusters
cluster_name_dict <- setNames(cluster_annots$cluster_name_w_num,
                              cluster_annots$lineage_cluster)
cluster_masc_res$cluster_name_w_num <- cluster_masc_res$lineage_cluster %>% recode(!!!cluster_name_dict)

# filter for only b and plasma clusters and plot
bcells_clusters <- obs %>%
  filter((lineage == 'B and plasma') & (lineage != 'Doublets and RBCs')) %>%
  pull(cluster_name_w_num) %>%
  unique()
masc_helper(obs,
            c('deg_case_control', 'control', 'case'),
            'cluster_name_w_num',
            fixed_effects = c('sex', 'ici_type'),
            group_subset = bcells_clusters,
            masc_res = cluster_masc_res,
            colors = c('slategray', 'tomato4'))
```

![](supplemental_figure_1_files/figure-gfm/supp_1m-5.png)<!-- -->

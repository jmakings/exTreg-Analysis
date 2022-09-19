library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")
library(collections)
library(tibble)
library(janitor)
library(clipr)

##### Subsetting patterns for exTregs with log2FC > 1 compared to other cell types #####

# opens csv file, subsets data for what you want, and writes files to the clipboard

read.csv("AntoineGSEA/exTreg_vs_Tn_Stats.csv")  %>% filter(expansion_score<0.05)  %>% filter(P_value < 0.05) %>% 
  filter(exTregs_by_Tn_log2FC > 1 ) %>% dplyr::select(pattern) %>% clipr::write_clip()

read.csv("AntoineGSEA/exTreg_vs_Treg_Stats.csv")  %>% filter(expansion_score<0.05)  %>% filter(P_value < 0.05) %>% 
  filter(exTregs_by_Tregs_log2FC > 1) %>% dplyr::select(pattern) %>% clipr::write_clip()

read.csv("AntoineGSEA/exTreg_vs_Th1_Stats.csv") %>% filter(expansion_score<0.05)  %>% filter(P_value < 0.05) %>% 
  filter(exTregs_by_Th1_log2FC > 1 ) %>% dplyr::select(pattern) %>% clipr::write_clip()


# This one also filters for log2FC values less than 1 (not used in analysis)

read.csv("AntoineGSEA/exTreg_vs_Tn_Stats.csv")  %>% filter(expansion_score<0.05)  %>% filter(P_value < 0.05) %>% 
  filter(exTregs_by_Tn_log2FC > 1 | exTregs_by_Tn_log2FC < -1) %>% dplyr::select(pattern) %>% clipr::write_clip()

read.csv("AntoineGSEA/exTreg_vs_Treg_Stats.csv")  %>% filter(expansion_score<0.05)  %>% filter(P_value < 0.05) %>% 
  filter(exTregs_by_Tregs_log2FC > 1 | exTregs_by_Tregs_log2FC < -1) %>% dplyr::select(pattern) %>% clipr::write_clip()

read.csv("AntoineGSEA/exTreg_vs_Th1_Stats.csv") %>% filter(expansion_score<0.05)  %>% filter(P_value < 0.05) %>% 
  filter(exTregs_by_Th1_log2FC > 1 | exTregs_by_Th1_log2FC < -1) %>% dplyr::select(pattern) %>% clipr::write_clip()

# Then, these clips are used at https://www.biotools.fr/misc/venny to create Venn Diagrams

##### Subsetting for all patterns that appear in exTregs #####

read.csv('AntoineGSEA/Raw_TCR_Contribution.csv') %>% filter(is.na(exTregs) == FALSE) %>% filter(is.na(Tregs) == FALSE) %>%
  select(pattern) %>% unique() -> tregsOnly

tregsOnly$pattern %>% clipr::write_clip()

read.csv('AntoineGSEA/Raw_TCR_Contribution.csv') %>% filter(is.na(exTregs) == FALSE) %>% filter(is.na(Th1) == FALSE) %>%
  select(pattern) %>% unique() -> th1Only
th1Only$pattern  %>% clipr::write_clip()

read.csv('AntoineGSEA/Raw_TCR_Contribution.csv') %>% filter(is.na(exTregs) == FALSE) %>% filter(is.na(Tn) == FALSE) %>%
  select(pattern) %>% unique() -> tnOnly
tnOnly$pattern  %>% clipr::write_clip()

read.csv('AntoineGSEA/Raw_TCR_Contribution.csv') %>% filter(is.na(exTregs) == FALSE) %>% filter(is.na(Tn) == TRUE) %>%
  filter(is.na(Th1) == TRUE) %>% filter(is.na(Tregs) == TRUE) %>%select(pattern) %>% unique() -> exTregExclusive

read.csv('AntoineGSEA/Raw_TCR_Contribution.csv') %>% filter(is.na(exTregs) == FALSE) %>%
  select(pattern) %>% unique() -> exTregOnly
exTregOnly$pattern  %>% clipr::write_clip()








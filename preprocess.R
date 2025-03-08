# --- Load Libraries -----------------
library(tidyverse)
library(readxl)
invisible(lapply(
  list.files(path = "~/jhu/projects/tidyexposomics/R/",
             pattern="*.R",full.names = TRUE),
  source))

# --- Exposure Definitions -------
des <- read.csv("~/jhu/projects/exposomePilot/results/description.csv")

# --- Exposures ----
# cytokine levels
aw_cytokines <- read_excel("~/jhu/projects/exposomePilot/data/AIRWEIGHS_Cytokine Results Sheet_REDCAP.xlsx",sheet = "Sheet2") %>% 
  dplyr::select(colnames(.)[grepl("id|visit|conc",colnames(.))]) %>% 
  filter(visit_cv=="CV2") %>% 
  mutate(aw_id=as.character(aw_id)) %>% 
  dplyr::select(-visit_cv)

# creatine levels
aw_creatine <- read_excel("~/jhu/projects/exposomePilot/data/Exposome_Creatinine Request_AIRWEIGHS MEBRL ESA Results 2021_Exposome Subset.xlsx") %>% 
  mutate(id=gsub("AW-|AW-0","",ID)) %>% 
  dplyr::select(id,`Creatinine (mg/dL)`)

# indoor pm 2.5
aw_pm25 <- read.csv("~/jhu/projects/exposomePilot/data/AW_PM25_2019-11-04.csv") %>% 
  pivot_wider(names_from = "redcap_event_name",
              values_from = "indoor_pm25_final") %>%
  setNames(c("id",
             "pm25_home_visit_1_arm_1",
             "pm25_home_visit_4_arm_1")) %>% 
  mutate(id=as.character(id))

# second hand smoke
aw_shs <- read.csv("~/jhu/projects/exposomePilot/data/AW_SHS_2019-11-04.csv") %>%
  pivot_wider(names_from = "redcap_event_name",
              values_from = "indoor_airnic") %>%
  setNames(c("id",
             "shs_home_visit_1_arm_1",
             "shs_home_visit_4_arm_1")) %>% 
  mutate(id=as.character(id))

# dust study
aw_dust <- read_excel("~/jhu/projects/exposomePilot/data/Shared AirWeighs Dust Study.xlsx") %>% 
  .[-1,] %>% 
  dplyr::select(colnames(.)[grepl("Record_ID|1-1",colnames(.))]) %>% 
  mutate(across(ends_with("-1"),as.numeric))

# chemicals
aw_chem <- read_excel("~/jhu/projects/exposomePilot/data/AWs_Chem_copy_w_counts.xlsx") %>% 
  .[1:195,] %>%  
  dplyr::select(colnames(.)[grepl(
    "\\bid\\b|\\bvisit\\b|_sg\\b|_cr\\b",colnames(.)
  )]) %>% 
  filter(visit==2) %>% 
  dplyr::select(-visit) %>% 
  mutate(across(ends_with("cr"),as.numeric)) %>% 
  mutate(across(ends_with("sg"),as.numeric)) 

# urine metals
aw_urine_metals <- read_excel(
  "~/jhu/projects/exposomePilot/data/Exposure final results.xlsx",
  sheet = "Urine") |> 
  dplyr::select(-`...2`) |> 
  # remove top two rows
  (\(df){df=df[-c(1,2),];df})() |> 
  mutate(id=gsub("AW-","",`...1`)) |> 
  filter(!grepl("DISC",id)) |> 
  # replace <LOD with NA
  mutate(across(everything(),~ifelse(.=="<LOD",NA,.))) |> 
  dplyr::select(-`...1`) |> 
  (\(df){names(df)=paste0("urine_",names(df));df})() |> 
  mutate_all(as.numeric) |> 
  mutate(urine_id=as.character(urine_id))

# serum metals
aw_serum_metals <- read_excel(
  "~/jhu/projects/exposomePilot/data/Exposure final results.xlsx",
  sheet = "Serum") |> 
  dplyr::select(-`...2`) |> 
  # remove top two rows
  (\(df){df=df[-c(1,2),];df})() |> 
  mutate(id=gsub("AW-","",`...1`)) |> 
  filter(!grepl("DISC",id)) |> 
  # replace <LOD with NA
  mutate(across(everything(),~ifelse(.=="<LOD",NA,.))) |> 
  dplyr::select(-c(`...1`,`Dilution factor in total`,`...20`,`NOTES (urine samples use)`)) |> 
  (\(df){names(df)=paste0("serum_",names(df));df})() |> 
  mutate_all(as.numeric) |> 
  mutate(serum_id=as.character(serum_id))
  
aw_dataset_cv2 <- read_excel(
  "~/jhu/projects/exposomePilot/data/Exposome_aw_dataset.xlsx",
  sheet = "CV2") 

# filter columns where more than X percent of the data are NA
# Function to filter columns with more than 20% NA values
filter_na_columns <- function(df, threshold = 0.2) {
  na_percent <- colMeans(is.na(df))
  df_filtered <- df[, na_percent <= threshold]
  return(df_filtered)
}

# Filter columns with more than 20% NA values
aw_dataset_cv2_filt <- aw_dataset_cv2 %>% 
  mutate( id=as.character(id)) %>% 
  left_join(.,
            aw_cytokines,
            by=c("id"="aw_id")) %>% 
  left_join(.,
            aw_creatine,
            by="id") %>% 
  left_join(.,
            aw_pm25,
            by="id") %>% 
  left_join(.,
            aw_shs,
            by="id") %>% 
  left_join(.,
            aw_dust,
            by=c("id"="Record_ID")) %>% 
  left_join(.,
            aw_chem,
            by="id") %>% 
  left_join(.,
            aw_urine_metals,
            by=c("id"="urine_id")) %>%
  left_join(.,
            aw_serum_metals,
            by=c("id"="serum_id")) %>%
  filter_na_columns(., threshold = 0.2) %>% 
  # https://www.aafp.org/pubs/afp/issues/2014/0301/p359.html
  mutate(fev1fvc_category=case_when(
    (pftfev1fvc_actual >= 0.75) & (pftfev1fvc_actual < 0.85) ~ "Moderate",
    (pftfev1fvc_actual < 0.75)  ~ "Severe",
    (pftfev1fvc_actual >= 0.85)  ~ "Normal"
  )) %>% 
  # https://pubmed.ncbi.nlm.nih.gov/17983880/
  mutate(ataq_category=case_when(
    ataqscore1<1~ "Normal",
    ataqscore1>=1 & ataqscore1 < 3  ~ "Uncontrolled_Asthma",
    ataqscore1>=3  ~ "Very_Uncontrolled_Asthma"
  )) %>% 
  # https://ww2.arb.ca.gov/resources/inhalable-particulate-matter-and-health
  mutate(pm25_category=case_when(
    pm25f<15 & !is.na(pm25f)~"Low",
    pm25f>15 & !is.na(pm25f)~"High",
    .default = NA
  )) %>% 
  mutate(pm25_stages=case_when(
    pm25f<10 & !is.na(pm25f)~"Low",
    pm25f>10 & pm25f<25 & !is.na(pm25f) ~"Medium",
    pm25f>25 & !is.na(pm25f)~"High",
    .default = NA
  )) %>% 
  mutate(pm10_category=case_when(
    pm10f<150  & !is.na(pm10f)~"Low",
    pm10f>150  & !is.na(pm10f)~"High",
    .default = NA
  )) %>% 
  as.data.frame() %>% 
  mutate(aw_id=paste0("s",id)) %>%
  `rownames<-`(.$aw_id) |> 
  mutate(female=case_when(
    aw_id == "s2809" ~ "male",
    .default=female
  )) |> 
  mutate(race=case_when(
    black=="black" ~ "black",
    black=="non-black" ~ "non_black",
  )) |> 
  mutate(fev1fvc_category=case_when(
    fev1fvc_category == "Normal" ~ "Unobstructed",
    fev1fvc_category == "Moderate" ~ "Moderate",
    fev1fvc_category == "Severe" ~ "Severe"
  )) |>
  mutate(fev1fvc_category=factor(fev1fvc_category,
                                 levels=c("Unobstructed","Moderate","Severe"))) |> 
  mutate(age_scaled=as.numeric(scale(age))) |> 
  mutate(sex=female) |> 
  mutate(sex=factor(sex,levels=c("female","male"))) |>
  mutate(race=factor(race,levels=c("non_black","black"))) |> 
  mutate(height=gli_height/100) |>
  mutate(fev_height=pftprefev1best/(height^2))


# --- RNA ----------
# get gene counts and meta data
gene_counts <- read.csv("~/jhu/projects/exposomePilot/results/gene_counts.csv") %>% 
  column_to_rownames("X") %>% 
  `colnames<-`(gsub(".*_","",colnames(.))) %>% 
  as.data.frame()

cd16_gene_meta <- read_excel(
  "~/jhu/projects/exposomePilot/data/AIRWEIGHS_RNA Extractions_20211109.xlsx",
  sheet = "CD16+ Non-classical ") %>% 
  dplyr::select(`Lab ID`,`Study ID`) %>% 
  mutate(cell_type="CD16+ Non-classical") %>% 
  mutate(cell_simple="CD16") %>% 
  mutate(`Lab ID`=paste0("S",`Lab ID`)) %>% 
  mutate(aw_id=paste0("s",`Study ID`))

cd4_gene_meta <- read_excel(
  "~/jhu/projects/exposomePilot/data/AIRWEIGHS_RNA Extractions_20211109.xlsx",
  sheet = "CD4+ T cells",skip = 1) %>% 
  dplyr::select(`Lab ID`,`Study ID`)%>% 
  mutate(cell_type="CD4+ T cells") %>% 
  mutate(cell_simple="CD4")%>% 
  mutate(`Lab ID`=paste0("S",`Lab ID`)) %>% 
  mutate(aw_id=paste0("s",`Study ID`))

# extract gene counts for cd16 and cd4
cd16_rna <- gene_counts %>% 
  dplyr::select(cd16_gene_meta$`Lab ID` %>% .[. %in% colnames(gene_counts)]) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(`Lab ID`=rownames(.)) %>%
  inner_join(.,
             cd16_gene_meta %>% 
               dplyr::select(`Lab ID`,`aw_id`),
             by="Lab ID") %>%
  column_to_rownames("aw_id") %>% 
  dplyr::select(-`Lab ID`) %>% 
  t() %>% 
  as.data.frame()

cd4_rna <- gene_counts %>% 
  dplyr::select(cd4_gene_meta$`Lab ID` %>% .[. %in% colnames(gene_counts)]) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(`Lab ID`=rownames(.)) %>%
  inner_join(.,
             cd4_gene_meta %>% 
               dplyr::select(`Lab ID`,`aw_id`),
             by="Lab ID") %>%
  column_to_rownames("aw_id") %>% 
  dplyr::select(-`Lab ID`) %>% 
  t() %>% 
  as.data.frame()

cd4_rna_fdata <- data.frame(gene=rownames(cd4_rna))
rownames(cd4_rna_fdata) <- cd4_rna_fdata$gene

cd16_rna_fdata <- data.frame(gene=rownames(cd16_rna))
rownames(cd16_rna_fdata) <- cd16_rna_fdata$gene

# --- Isoform ----------
isoform_counts <- read.csv("~/jhu/projects/exposomePilot/results/tx_counts.csv") %>% 
  column_to_rownames("X") %>% 
  `colnames<-`(gsub(".*_","",colnames(.))) %>% 
  as.data.frame()

# extract isoform counts for cd16 and cd4
cd16_isoform <- isoform_counts %>% 
  dplyr::select(cd16_gene_meta$`Lab ID` %>% .[. %in% colnames(isoform_counts)]) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(`Lab ID`=rownames(.)) %>%
  inner_join(.,
             cd16_gene_meta %>% 
               dplyr::select(`Lab ID`,`aw_id`),
             by="Lab ID") %>%
  column_to_rownames("aw_id") %>% 
  dplyr::select(-`Lab ID`) %>% 
  t() %>% 
  as.data.frame()

cd4_isoform <- isoform_counts %>% 
  dplyr::select(cd4_gene_meta$`Lab ID` %>% .[. %in% colnames(isoform_counts)]) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(`Lab ID`=rownames(.)) %>%
  inner_join(.,
             cd4_gene_meta %>% 
               dplyr::select(`Lab ID`,`aw_id`),
             by="Lab ID") %>%
  column_to_rownames("aw_id") %>% 
  dplyr::select(-`Lab ID`) %>% 
  t() %>% 
  as.data.frame()

# --- miRNA ----------

mirna_counts <- read.csv("~/jhu/projects/exposomePilot/results/mirna.csv") %>% 
  column_to_rownames("X") %>% 
  as.data.frame() %>% 
  `colnames<-`(gsub(".*_","",colnames(.)))

# extract miRNA counts for cd16 and cd4
cd16_mirna <- mirna_counts %>% 
  dplyr::select(cd16_gene_meta$`Lab ID` %>% .[. %in% colnames(mirna_counts)]) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(`Lab ID`=rownames(.)) %>%
  inner_join(.,
             cd16_gene_meta %>% 
               dplyr::select(`Lab ID`,`aw_id`),
             by="Lab ID") %>%
  column_to_rownames("aw_id") %>% 
  dplyr::select(-`Lab ID`) %>% 
  t() %>% 
  as.data.frame()

cd4_mirna <- mirna_counts %>% 
  dplyr::select(cd4_gene_meta$`Lab ID` %>% .[. %in% colnames(mirna_counts)]) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(`Lab ID`=rownames(.)) %>%
  inner_join(.,
             cd4_gene_meta %>% 
               dplyr::select(`Lab ID`,`aw_id`),
             by="Lab ID") %>%
  column_to_rownames("aw_id") %>% 
  dplyr::select(-`Lab ID`) %>% 
  t() %>% 
  as.data.frame()

cd4_mirna_fdata <- data.frame(gene=rownames(cd4_mirna))
rownames(cd4_mirna_fdata) <- cd4_mirna_fdata$gene

cd16_mirna_fdata <- data.frame(gene=rownames(cd16_mirna)) 
rownames(cd16_mirna_fdata) <- cd16_mirna_fdata$gene

# --- Protein ----------

prot <- readxl::read_excel("~/jhu/projects/exposomePilot/data/proteomics.xlsx")

prot_fdata <- prot[,1:2] %>% 
  as.data.frame() %>% 
  `row.names<-`(.$Protein.Group) |> 
  mutate(gene=Genes)
rownames(prot_fdata) <- prot_fdata$Protein.Group

prot_abd <- prot %>% dplyr::select(
  c("Protein.Group",colnames(prot)[grepl("E.._",colnames(prot))])
) %>% 
  column_to_rownames("Protein.Group")

prot_abd[prot_abd=="nd"] <- NA

prot_abd <- prot_abd %>% 
  mutate(across(all_of(colnames(.)), as.numeric)) %>% 
  dplyr::select(colnames(.)[grepl("Asth",colnames(.))]) %>% 
  `colnames<-`(gsub("_.*","",colnames(.)))

prot_meta <- read_excel("~/jhu/projects/exposomePilot/data/Exposomics_key.xlsx") %>% 
  filter(grepl("AW",Sample_ID)) %>% 
  mutate(Sample_ID=gsub("AW","",Sample_ID)) %>% 
  as.data.frame() %>% 
  left_join(.,
            aw_dataset_cv2_filt,
            by=c("Sample_ID"="id")) %>% 
  mutate(aw_id=paste0("s",Sample_ID)) %>% 
  `rownames<-`(.$Inj_ID) 

prot_abd <- prot_abd %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Inj_ID=rownames(.)) %>%
  inner_join(.,
             prot_meta %>% 
               dplyr::select(Inj_ID,aw_id),
             by="Inj_ID") %>%
  column_to_rownames("aw_id") %>% 
  dplyr::select(-Inj_ID) %>% 
  t() %>% 
  as.data.frame()

prot_abd_unnorm <- prot_abd |> 
  (\(x) 2^x)()

# --- Adductomics ------------------
adduct_long <- readRDS("~/jhu/projects/exposomePilot/data/Adductomics/FennaExposome_long.rds")

adduct <- adduct_long |> 
  dplyr::select(IonIntQuant_key,Subject,Raw_IonIntensity) |> 
  mutate(Subject=paste("s",Subject,sep="")) |> 
  filter(Subject %in% aw_dataset_cv2_filt$aw_id) |> 
  pivot_wider(
    names_from = Subject,
    values_from = Raw_IonIntensity
  ) |> 
  column_to_rownames("IonIntQuant_key") 

adduct_fdata <- adduct_long |> 
  dplyr::select(IonIntQuant_key,Protein_ID,
                Gene_name,residue,mean_mass,
                annotation) |> 
  distinct() |> 
  mutate(gene=Gene_name)

rownames(adduct_fdata) <- adduct_fdata$IonIntQuant_key

# --- Save Data ---------------------
save(
  # information and meta data
  des,
  aw_dataset_cv2_filt,
  
  # omics
  cd4_rna,
  cd16_rna,
  cd4_mirna,
  cd16_mirna,
  prot_abd_unnorm,
  adduct,
  
  # rowdata
  cd4_rna_fdata,
  cd16_rna_fdata,
  cd4_mirna_fdata,
  cd16_mirna_fdata,
  prot_fdata,
  adduct_fdata,
  
  file = "~/jhu/projects/tidyexposomics/data/expom.RData")


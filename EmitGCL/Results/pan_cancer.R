setwd("/users/PAS1475/anjunma/wxy/metastasis/")
#########breast
breast_LN5 <- read.csv("Breast_cancer5_LN.csv")
breast_P5 <- read.csv("Breast_cancer5_P.csv")

breast_LN9 <- read.csv("Breast_cancer9 LN.csv")
breast_P9 <- read.csv("Breast_cancer9 P.csv")

merged_df <- full_join(breast_LN5, breast_LN9, by = "index")
breast_LN <- merged_df %>%
  group_by(index) %>%
  summarise(Count = sum(Frequency.x, Frequency.y, na.rm = TRUE),
            Score = max(Attention.x, Attention.y, na.rm = TRUE))
breast_LN <- as.data.frame(breast_LN)
breast_LN <- breast_LN %>%
  filter(is.finite(Score))

merged_df <- full_join(breast_P5, breast_P9, by = "index")
breast_P <- merged_df %>%
  group_by(index) %>%
  summarise(Count = sum(Frequency.x, Frequency.y, na.rm = TRUE),
            Score = max(Attention.x, Attention.y, na.rm = TRUE))
breast_P <- as.data.frame(breast_P)
breast_P <- breast_P %>%
  filter(is.finite(Score))

breast_P <- breast_P[breast_P$Count>13,]
breast_LN <- breast_LN[breast_LN$Count>13,]
#########head and neck
HN_P <- read.csv("Head and neck cancer_P.csv")
HN_P <- HN_P %>%
  filter(is.finite(Attention))

HN_LN <- read.csv("Head and neck cancer_LN.csv")
HN_LN <- HN_LN %>%
  filter(is.finite(Attention))

HN_P <- HN_P[HN_P$Frequency>4,]
HN_LN <- HN_LN[HN_LN$Frequency>4,]
################HRA
HRA_P <- read.csv("HRA_P.csv")
HRA_P <- HRA %>%
  filter(is.finite(Attention))

HRA_LN <- read.csv("HRA_LN.csv")
HRA_LN <- HRA_LN %>%
  filter(is.finite(Attention))

HRA_LN <- HRA_LN[HRA_LN$Frequency>3,]
HRA_P <- HRA_P[HRA_P$Frequency>3,]
###############brain cancer
brain_P <- read.csv("Brain cancer_P.csv")
brain_P <- brain_P %>%
  filter(is.finite(Attention))

brain_LN <- read.csv("Brain cancer_LN.csv")
brain_LN <- brain_LN %>%
  filter(is.finite(Attention))

brain_LN <- brain_LN[brain_LN$Frequency>3,]
brain_P <- brain_P[brain_P$Frequency>3,]
########Pancreastic cancer
Pancreastic_P <- read.csv("Pancreastic_P.csv")
Pancreastic_P <- Pancreastic_P %>%
  filter(is.finite(Attention))

Pancreastic_LN <- read.csv("Pancreastic_LN.csv")
Pancreastic_LN <- Pancreastic_LN %>%
  filter(is.finite(Attention))

Pancreastic_LN <- Pancreastic_LN[Pancreastic_LN$Frequency>2,]
Pancreastic_P <- Pancreastic_P[Pancreastic_P$Frequency>2,]
########Papillary cancer
Papillary_P <- read.csv("Papillary thyroid cancer_P.csv")
Papillary_P <- Papillary_P %>%
  filter(is.finite(Attention))

Papillary_LN <- read.csv("Papillary thyroid cancer_LN.csv")
Papillary_LN <- Papillary_LN %>%
  filter(is.finite(Attention))

df_p <- data.frame(gene = c(breast_P$index,HN_P$index,brain_P$index,Pancreastic_P$index,Papillary_P$index),
                   att = c(breast_P$Score,HN_P$Attention,brain_P$Attention,Pancreastic_P$Attention,Papillary_P$Attention),
                   fre = c(breast_P$Count,HN_P$Frequency,brain_P$Frequency,Pancreastic_P$Frequency,Papillary_P$Frequency),
                   cancer = c(rep("Breast",length(breast_P$index)),rep("HN",length(HN_P$index)),rep("Brain",length(brain_P$index)),
                   rep("Pancreastic",length(Pancreastic_P$index)),rep("Papillary",length(Papillary_P$index))))
df_p <- df_p[-grep("^RP", df_p$gene),]
df_p <- df_p[df_p$att>0.1,]
library(dplyr)
library(purrr)

# 创建每个癌症类型的基因集合
gene_sets <- df_p %>% 
  group_by(cancer) %>%
  summarise(gene = list(unique(gene)), .groups = 'drop')

gene_list <- setNames(gene_sets$gene, gene_sets$cancer)
common_genes_p <- Reduce(intersect, gene_list)

# 计算每个癌症类型的特异基因
df_p[df_p$gene %in% names(table(df_p$gene)[table(df_p$gene)==1]),]
specific_genes <- map(gene_list, function(x) setdiff(x, common_genes))

dfp <- df_p[df_p$gene %in% c(names(table(df_p$gene)[table(df_p$gene)==1]),common_genes),]
dfp$cancer <- factor(dfp$cancer, levels = unique(dfp$cancer))

ggplot(dfp, aes(x= cancer ,y = gene ,color = att)) +
  geom_point()+scale_color_gradient(low = "#6DCADA", high = "#07A9CA")+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


df_ln <- data.frame(gene = c(breast_LN$index,HN_LN$index,brain_LN$index,Pancreastic_LN$index,Papillary_LN$index),
                   att = c(breast_LN$Score,HN_LN$Attention,brain_LN$Attention,Pancreastic_LN$Attention,Papillary_LN$Attention),
                   fre = c(breast_LN$Count,HN_LN$Frequency,brain_LN$Frequency,Pancreastic_LN$Frequency,Papillary_LN$Frequency),
                   cancer = c(rep("Breast",length(breast_LN$index)),rep("HN",length(HN_LN$index)),rep("Brain",length(brain_LN$index)),
                              rep("Pancreastic",length(Pancreastic_LN$index)),rep("Papillary",length(Papillary_LN$index))))
df_ln <- df_ln[-grep("^RP", df_ln$gene),]
df_ln <- df_ln[df_ln$att>0.1,]
library(dplyr)
library(purrr)

# 创建每个癌症类型的基因集合
gene_sets <- df_ln %>% 
  group_by(cancer) %>%
  summarise(gene = list(unique(gene)), .groups = 'drop')

gene_list <- setNames(gene_sets$gene, gene_sets$cancer)
common_genes_l <- Reduce(intersect, gene_list)

# 计算每个癌症类型的特异基因
#df_p[df_p$gene %in% names(table(df_p$gene)[table(df_p$gene)==1]),]
#specific_genes <- map(gene_list, function(x) setdiff(x, common_genes))

dfl <- df_ln[df_ln$gene %in% c(names(table(df_ln$gene)[table(df_ln$gene)==1]),common_genes),]
ggplot(dfl, aes(x= cancer ,y = gene ,color = att)) +
  geom_point()+scale_color_gradient(low = "#EBCAE0", high = "#E55CB7")+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dfp



# 转换为列表，方便后续处理
gene_list <- setNames(gene_sets$Genes, gene_sets$CancerType)
list_input <- list(
  Breast = df_p$gene[df_p$cancer=="Breast"],
  HN = df_p$gene[df_p$cancer=="HN"],
  Brain = df_p$gene[df_p$cancer=="Brain"],
  Pancreastic = df_p$gene[df_p$cancer=="Pancreastic"],
  Papillary = df_p$gene[df_p$cancer=="Papillary"]
)
upset(fromList(list_input), order.by = "freq")

## cell number 
for (sample in c("1","2","4","5","6","7","8","3")){
  print(sample)
  site_label_path <- paste0("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE167036/Patient_subsample/Patient",sample,"/label.csv")
  label <- read.csv(site_label_path)
  print(table(label$orig.ident))
}

for (sample in c("patientD2020_2","patientD2020_3")){
  print(sample)
  site_label_path <- paste0("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE180286/All_patient/",sample,"/D2020_cell_label.csv")
  site_label <- read.csv(site_label_path)
  label <- site_label$X0
  print(table(label))
}
for (sample in c("patientE2020_2","patientE2020_3")){
  print(sample)
  site_label_path <- paste0("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE180286/All_patient/",sample,"/E2020_cell_label.csv")
  site_label <- read.csv(site_label_path)
  label <- site_label$X0
  print(table(label))
}
for (sample in c("patientC2020_2","patientC2020_3")){
  print(sample)
  site_label_path <- paste0("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE180286/All_patient/",sample,"/C2020_cell_label.csv")
  site_label <- read.csv(site_label_path)
  label <- site_label$X0
  print(table(label))
}
for (sample in c("patientB2019_2","patientB2019_3")){
  print(sample)
  site_label_path <- paste0("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE180286/All_patient/",sample,"/B2019_cell_label.csv")
  site_label <- read.csv(site_label_path)
  label <- site_label$X0
  print(table(label))
}
for (sample in c("patientA2019_2","patientA2019_3")){
  print(sample)
  site_label_path <- paste0("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE180286/All_patient/",sample,"/A2019_cell_label.csv")
  site_label <- read.csv(site_label_path)
  label <- site_label$X0
  print(table(label))
}



for (sample in c("Patient1","Patient2")){
  site_label_path <- paste0("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE198099/",sample,"/Data/label.csv")
  site_label <- read.csv(site_label_path)
  print(table(site_label$Label))
}

for (sample in c("R1","R2","R3","R4")){
  site_label_path <- paste0("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE162631/All_patient/",sample,"/label.csv")
  site_label <- read.csv(site_label_path)
  print(table(site_label$Label))
}
for (sample in c("HN237","HN242","HN251","HN257","HN263","HN272","HN279")){

site_label_path <- paste0("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE188737/",sample,"/cell_label.csv")
site_label <- read.csv(site_label_path)
print(table(site_label$P_Mid))
}
for (sample in c("Patient1","Patient2","Patient3")){
  site_label_path <- paste0("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE197177/Patient_subsample/",sample,"/cell_label.csv")
  site_label <- read.csv(site_label_path)
  print(table(site_label))
}

label <- read.csv("/fs/ess/PAS1475/Maoteng/Metastasis_new/GSE241184/label.csv")
table(label$Label)

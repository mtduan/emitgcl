#Gene signature
library(dplyr)
library(purrr)
setwd("C:/Users/Wang Xiaoying/BMBL Dropbox/BMBL Project/Xiaoying-Wang_Projects/4-Metastasis/Results/breast cancer/early metastastic gene signature")
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

breast_P <- breast_P[-grep("^RP", breast_P$index),]
breast_LN <- breast_LN[-grep("^RP", breast_LN$index),]

breast_P <- breast_P[breast_P$Count>12,]
breast_LN <- breast_LN[breast_LN$Count>12,]
breast_P <- breast_P[breast_P$Score>0.1,]
breast_LN <- breast_LN[breast_LN$Score>0.1,]
breast <- rbind(breast_P,breast_LN)
breast$site <- c(rep("Primary",dim(breast_P)[1]),rep("Metastatic",dim(breast_LN)[1]))
#breast_P[breast_P$Count>9,]$index
library(ggplot2)
ggplot(breast, aes(x= site ,y = index ,size = Score)) +
  geom_point()+scale_color_gradient(low = "#EBCAE0", high = "#E55CB7")+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
###REGULON 
regulon <- function(E_name, R_name, L_name){
  Encode <- read.csv(paste0("C:/Users/Wang Xiaoying/BMBL Dropbox/BMBL Project/Xiaoying-Wang_Projects/4-Metastasis/Results/breast cancer/early metastastic gene signature/",E_name,".tsv"),header=T, sep="\t")
  ReMap <-  read.csv(paste0("C:/Users/Wang Xiaoying/BMBL Dropbox/BMBL Project/Xiaoying-Wang_Projects/4-Metastasis/Results/breast cancer/early metastastic gene signature/",R_name,".tsv"),header=T, sep="\t")
  Literature <-  read.csv(paste0("C:/Users/Wang Xiaoying/BMBL Dropbox/BMBL Project/Xiaoying-Wang_Projects/4-Metastasis/Results/breast cancer/early metastastic gene signature/",L_name,".tsv"),header=T, sep="\t")
  a <- Encode[Encode$FDR<0.001,]
  b <- ReMap[ReMap$FDR<0.001,]
  c <- Literature[Literature$FDR<0.001,]
  TF <- intersect(intersect(a$TF,b$TF),c$TF)
  regulon_l <- list()
  for (i in (1:nrow(a))){
    regulon_l[[a$TF[i]]] <- unique(unlist(strsplit(a$Overlapping_Genes[i],",")) )
  }
  for (i in (1:nrow(b))){
    if (!(b$TF[i] %in% names(regulon_l)))
      regulon_l[[b$TF[i]]] <- unique(unlist(strsplit(b$Overlapping_Genes[i],",")) )
    else{
      regulon_l[[b$TF[i]]]<-unique(c(regulon_l[[b$TF[i]]], regulon_l[[b$TF[i]]] ) )
    }
  }
  for (i in (1:nrow(c))){
    if (!(c$TF[i] %in% names(regulon_l)))
      regulon_l[[c$TF[i]]] <- unique(unlist(strsplit(c$Overlapping_Genes[i],",")) )
    else{
      regulon_l[[c$TF[i]]]<- unique(c(regulon_l[[c$TF[i]]], regulon_l[[c$TF[i]]] ) )
    }
  }
  return (regulon_l)
}
P_regulon <- regulon("ENCODE_ChIP-seq_P","Literature_ChIP-seq_P","ReMap_ChIP-seq_P")
L_regulon <- regulon("ENCODE_ChIP-seq_LN","Literature_ChIP-seq_LN","ReMap_ChIP-seq_LN")

#P_regulon <- regulon("ENCODE_early_RP","Literature_early_RP","ReMap_early_RP")
#L_regulon <- regulon("ENCODE_ChIP-seq_LN","Literature_ChIP-seq_LN","ReMap_ChIP-seq_LN")


df <- data.frame()
#(df) <- c('Regulator','Attention','Group') 
for (i in intersect(names(P_regulon),names(L_regulon))){
  #for (i in intersect(names(PM_regulon),names(LM_regulon))){
  att <- breast_P$Score[breast_P$index %in% P_regulon[[i]]]
  #att <- c(att,rep(0, (length(P257_regulon[[i]])-length(att))))
  df1 <-data.frame(cbind(cbind(i,as.numeric(att)),'P'))
  att2 <- breast_LN$Score[breast_LN$index %in% L_regulon[[i]]]
  #att2 <- c(att,rep(0, (length(P257_regulon[[i]])-length(att2))))
  df2 <-data.frame(cbind(cbind(i,as.numeric(att2)),'L'))
  names(df1) <- c('Regulator','Attention','Group')
  names(df2) <- c('Regulator','Attention','Group')
  df <- rbind(df, df1)
  df <- rbind(df, df2)
  #TF_att_272[[i]] <- mean(P272$Interest.Mean.Attention.Value[P272$Gene %in% P257_regulon[[i]]])
}
df$Attention <- as.numeric(df$Attention)
library(dplyr)
df1 <- data.frame(df %>%
                    group_by(Regulator,Group) %>%
                    summarise(MeanValue = mean(Attention, na.rm = TRUE)))

DF <- df1[df1$Regulator %in% c("YY1","STAT1","NFE2L2","JUND","JUNB","JUN","FOXP1","FOS","CEBPD","CEBPB","BHLHE40","ATF3"),]
library(ggplot2)
# 创建箱形图
#ggplot(df1, aes(x= Group ,y = Regulator ,color = MeanValue, size = )) +
#  geom_point()+scale_color_gradient(low = "orange", high = "red")+theme_minimal()+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#df1[df1$Group=='L',][order(df1[df1$Group=='L',]$MeanValue),]

merged_df <- merge(df1[df1$Group=='L',], df1[df1$Group=='P',], by = "Regulator", suffixes = c(".df1", ".df2"))
merged_df$score_diff <- merged_df$MeanValue.df1 - merged_df$MeanValue.df2
merged_df$score_ratio <- merged_df$MeanValue.df2 / merged_df$MeanValue.df1




dfa <- df1[df1$Group=='L',]
dfb <- df1[df1$Group=='P',]
dfa$rank_df1 <- rank(-dfa$MeanValue)  
dfb$rank_df2 <- rank(-dfb$MeanValue)   
merged_df <- merge(dfa[, c("Regulator", "rank_df1")], dfb[, c("Regulator", "rank_df2")], by = "Regulator")

print(merged_df)


#merged_df[order(merged_df$score_ratio),]
#merged_df[order(merged_df$score_diff),]  
merged_df$rank_change <- merged_df$rank_df2 - merged_df$rank_df1
ggplot(merged_df, aes(x = Regulator, y = rank_change, fill = Regulator)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "score rank variation", y = "rank variation", x = "gene")


########attention heatmap
data_normalized <- df1 %>%
  group_by(Group) %>%
  mutate(Value = (MeanValue - min(MeanValue)) / (max(MeanValue) - min(MeanValue))) %>%
  ungroup()
ggplot(data_normalized, aes(x = Group, y = Regulator, fill = Value)) +
  geom_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(fill = "Value", title = "Heatmap Example", x = "Column", y = "Row")
###############attention total

dfa <- df1[df1$Group=='L',]
dfb <- df1[df1$Group=='P',]
dfa$rank_df1 <- rank(-dfa$MeanValue)  
dfb$rank_df2 <- rank(-dfb$MeanValue)  
merged_df <- merge(dfa[, c("Regulator", "rank_df1")], dfb[, c("Regulator", "rank_df2")], by = "Regulator")
dfa$total <- dfa$MeanValue+dfb$MeanValue
DF <- dfa[dfa$Regulator %in% c("YY1","STAT1","NFE2L2","JUND","JUNB","JUN","FOXP1","FOS","CEBPD","CEBPB","BHLHE40","ATF3"),]
library(RColorBrewer)
color_palette <- colorRampPalette(brewer.pal(9, "Greens"))(100) 
ggplot(DF, aes(x = Group, y = Regulator, fill = total,)) +
  geom_tile() +  
  scale_fill_gradientn(colors = color_palette,limit = c(0.9, 1.1),oob = scales::oob_squish )+
  #scale_fill_gradient2(low = "blue", high = "red",midpoint =0.8) +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(fill = "Value", title = "Heatmap Example", x = "Column", y = "Row")
  
net <- function(TF){
  example_list_intersect <- P_regulon[ c("YY1","STAT1","NFE2L2","JUND","JUNB","JUN","FOXP1","FOS","CEBPD","CEBPB","BHLHE40","ATF3")]
  
  df_list <- lapply(names(example_list_intersect), function(name) {
    data.frame(source = name, target = example_list_intersect[[name]], stringsAsFactors = FALSE)
  })
  df_MP <- do.call(rbind, df_list)
  df_MP$group <- rep("MP",nrow(df_MP))
  example_list_intersect <- L_regulon[ c("YY1","STAT1","NFE2L2","JUND","JUNB","JUN","FOXP1","FOS","CEBPD","CEBPB","BHLHE40","ATF3")]
  df_list <- lapply(names(example_list_intersect), function(name) {
    data.frame(source = name, target = example_list_intersect[[name]], stringsAsFactors = FALSE)
  })
  df_LP <- do.call(rbind, df_list)
  df_LP$group <- rep("ML",nrow(df_LP))
  
  df <- rbind(df_MP,df_LP)
  #df <- df[-grep("^RP", df$target),]
  
  library(dplyr)
  get_group <- function(col1, col2, col3) {
    ifelse(duplicated(cbind(col1, col2)) | duplicated(cbind(col1, col2), fromLast = TRUE), "new", col3)
  }
  
  df$group <- get_group(df$source, df$target, df$group)
  
  df <- df[!duplicated(df[c("source", "target","group")]), ]
  
  library(igraph)
  library(gridExtra)
  #df <- df1 
  #df <- df[df$source %in% c('HIF1A','HSF1','RUNX1','FOXP2','NFKB1','STAT1','MAFB','SOX11'),]
  network <- graph_from_data_frame(df, directed = TRUE)
  print(df)
  #all_nodes <- unique(c(df$source, df$target))
  #node_types <- ifelse(all_nodes %in% df$source, "source", "target")
  #V(network)$type <- ifelse(V(network)$name %in% df$source, TRUE, FALSE)
  color_map <- setNames(c("#DCDCDC", "#f4b9a7", "#b0d1e0"), c("new", "MP", "ML"))
  E(network)$color <- sapply(E(network)$group, function(l) color_map[l])
  #(network)$color <- ifelse(V(network)$name %in% df$source, "#f4b9a7", E(network)$color)
  V(network)[ends(network, E(network))[,2]]$color <- E(network)$color
  V(network)[ends(network, E(network))[,1]]$color <- "#DCDCDC"
  return(network)
}

par(mfrow = c(3, 4),mar=rep(0.1,4), oma=rep(4,4))


library(igraph)
network <- net(c(1))
plot(network, edge.arrow.size=0.8,
     vertex.size=12, vertex.label.color="black", vertex.label.cex=0.8,vertex.label.font = 2,
     edge.color=E(network)$color,vertex.color=V(network)$color,
     vertex.frame.color="#DCDCDC", main = )
network <- net()
plot(network, edge.arrow.size=0.1,
     vertex.size=12, vertex.label.color="black", vertex.label.cex=0.8,vertex.label.font = 2,
     edge.color=E(network)$color,vertex.color=V(network)$color,
     vertex.frame.color="#DCDCDC")



#Fig6.b
example_list_intersect <- P_regulon[ c("YY1","STAT1","NFE2L2","JUND","JUNB","JUN","FOXP1","FOS","CEBPD","CEBPB","BHLHE40","ATF3")]

df_list <- lapply(names(example_list_intersect), function(name) {
  data.frame(source = name, target = example_list_intersect[[name]], stringsAsFactors = FALSE)
})
df_MP <- do.call(rbind, df_list)
df_MP$group <- rep("MP",nrow(df_MP))
example_list_intersect <- L_regulon[ c("YY1","STAT1","NFE2L2","JUND","JUNB","JUN","FOXP1","FOS","CEBPD","CEBPB","BHLHE40","ATF3")]
df_list <- lapply(names(example_list_intersect), function(name) {
  data.frame(source = name, target = example_list_intersect[[name]], stringsAsFactors = FALSE)
})
df_LP <- do.call(rbind, df_list)
df_LP$group <- rep("ML",nrow(df_LP))

df <- rbind(df_MP,df_LP)
#df <- df[-grep("^RP", df$target),]

library(dplyr)
get_group <- function(col1, col2, col3) {
  ifelse(duplicated(cbind(col1, col2)) | duplicated(cbind(col1, col2), fromLast = TRUE), "new", col3)
}

df$group <- get_group(df$source, df$target, df$group)

df <- df[!duplicated(df[c("source", "target","group")]), ]
E(network)$color <- ifelse(E(network)$group == "new", "gray",ifelse(E(network)$group == "MP", "#6DCADA", "#EBCAE0") )

V(network)$type <- bipartite_mapping(network)$type
g <- network

inner_radius <- 1
outer_radius <- 2


source_nodes <- unique(df$source)
target_nodes <- unique(df$target)

layout_circle <- matrix(nrow = vcount(g), ncol = 2)
names <- V(g)$name
source_angles <- seq(0, 2 * pi, length.out = length(source_nodes) + 1)
target_angles <- seq(0, 2 * pi, length.out = length(target_nodes) + 1)


for (i in 1:length(names)) {
  if (names[i] %in% source_nodes) {
    angle_index <- which(source_nodes == names[i])
    layout_circle[i, ] <- c(inner_radius * cos(source_angles[angle_index]), 
                            inner_radius * sin(source_angles[angle_index]))
  } else {
    angle_index <- which(target_nodes == names[i])
    layout_circle[i, ] <- c(outer_radius * cos(target_angles[angle_index]), 
                            outer_radius * sin(target_angles[angle_index]))
  }
}

windowsFonts(A=windowsFont("Arial"))
par(family="A")
plot(g, layout = layout_circle, vertex.color = ifelse(V(g)$name %in% source_nodes, "#96CCCB", "#CFEAF1"), 
     vertex.frame.color = "gray", vertex.size = 20, 
     vertex.label.color = "black", vertex.label.cex = 0.6, 
     edge.color = E(g)$color, edge.lty = 2,edge.arrow.size = 0.2)


#Fig6.c
DF$degree <- degree(g)[c("YY1","STAT1","NFE2L2","JUND","JUNB","JUN","FOXP1","FOS","CEBPD","CEBPB","BHLHE40","ATF3")]
DF$closeness <- closeness(g)[c("YY1","STAT1","NFE2L2","JUND","JUNB","JUN","FOXP1","FOS","CEBPD","CEBPB","BHLHE40","ATF3")]


color_palette <- colorRampPalette(brewer.pal(9, "Oranges"))(100) 
ggplot(DF, aes(x = Group, y = Regulator, fill = degree,)) +
  geom_tile() +  
  scale_fill_gradientn(colors = color_palette,limit = c(8, 17),oob = scales::oob_squish )+
  #scale_fill_gradient2(low = "blue", high = "red",midpoint =0.8) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  labs(fill = "Value", title = "Heatmap Example", x = "Column", y = "Row")

color_palette <- colorRampPalette(brewer.pal(9, "Reds"))(100) 
ggplot(DF, aes(x = Group, y = Regulator, fill = closeness,)) +
  geom_tile() + 
  scale_fill_gradientn(colors = color_palette,limit = c(0, 0.1),oob = scales::oob_squish )+
  #scale_fill_gradient2(low = "blue", high = "red",midpoint =0.8) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(fill = "Value", title = "Heatmap Example", x = "Column", y = "Row")
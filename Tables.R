# FUNCTIONS
desat=function(cols, sat=0.5) {
  X=diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

using=function(...) {
  libs=unlist(list(...))
  req=unlist(lapply(libs,require,character.only=TRUE))
  need=libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

get_percent_col=function(table_melt) {
  table_melt["percent"]=0
  for (sample in unique(table_melt$Sample)){
    total=sum(table_melt[table_melt$Sample== sample, 3])
    for (genus in unique(table_melt$Taxa)){
      table_melt[which(table_melt$Sample== sample & table_melt$Taxa== genus), 4]=100 * table_melt[which(table_melt$Sample== sample & table_melt$Taxa== genus), 3] / total
    }
  }
  return(data.frame(table_melt))
}

percent_table=function(table_melt, manyX=0.05, table_clean) {
  for (sample in unique(table_melt$Sample)) {
    value_sum=0
    threshold=sum(table_melt[table_melt$Sample==sample,3])*manyX
    for (genus in unique(table_melt$Taxa)) {
      value=table_melt[which(table_melt$Sample== sample & table_melt$Taxa== genus), 3]
      if (length(value) > 0){
        if (value >= threshold) {
          table_clean[genus , sample]=value
        }
        else {
          value_sum=value_sum + value
          table_clean[genus , sample]=0
        }
      }
    }
    table_clean["0thers", sample]=value_sum
  } 
  return(as.matrix(table_clean))
}

#PACKAGES
using("ggplot2","RColorBrewer", "randomcoloR", "reshape2", "dplyr", "readr", "stringr")



# STACKED BARPLOT Phyla
## input
V3V4_D2=read_delim("V3V4_MeanDepth.csv", 
                      delim="\t", escape_double=FALSE, 
                      col_names=FALSE, col_types=cols(X2=col_number(), 
                                                          X3=col_number(), X4=col_number(), 
                                                          X5=col_number(), X6=col_number()), 
                      na="0", trim_ws=TRUE)
colnames(V3V4_D2)=c("Taxa","N1","N2","T1","T2")
##  melt data
V3V4_D2$Taxa=str_replace(V3V4_D2$Taxa, "^.+;D_3__", "") 
V3V4_D2$Taxa=str_replace_all(V3V4_D2$Taxa, "D_.__", "")
melt_V3V4_D2=melt(V3V4_D2, id.vars="Taxa", variable.name="Sample")
unique_taxa=unique(melt_V3V4_D2)
## How many X
manyX=0.01 
## Table needed
melt_V3V4_D2=melt_V3V4_D2 %>% replace(is.na(.), 0) #replace NA by 0
empty_D2=matrix(nrow=length(unique_taxa)+1, ncol=length(unique(melt_V3V4_D2$Sample))) #Empty matrix
colnames(empty_D2)=unique(melt_V3V4_D2$Sample)
rownames(empty_D2)=c(unique_taxa,"0thers") #Add "Others" as Taxa
clean_V3V4_D2=percent_table(melt_V3V4_D2, manyX, empty_D2) #Get only > X % of sample
clean_V3V4_D2[clean_V3V4_D2==0]=NA #0 replaced by NA
R2=rownames(clean_V3V4_D2[rowSums(is.na(clean_V3V4_D2)) != ncol(clean_V3V4_D2), ]) #Get Taxa with not only NA as value for every sample
clean_V3V4_D2=clean_V3V4_D2 %>% replace(is.na(.), 0) #replace NA by 0
meltC_V3V4_D2=melt(clean_V3V4_D2, id.vars="Taxa", variable.name="Sample") #Melt
## Color
qual_col_pals=brewer.pal.info[brewer.pal.info$category== 'qual',] #Colors for ggplot
col_vector=unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols_ggplot=desat(c(brewer.pal(n=length(unique_taxa)-1), name="Set3"),"White"),0.4) #Colors with "Others" + Desat
## Plot
ggplot_V3V4=ggplot(data=meltC_V3V4_D2, aes(x=Var2, y=value, fill=Var1)) + #V3V4 ggplot2
  coord_flip() +
  geom_bar(position="fill", stat="identity", na.rm=TRUE, colour="darkgrey", lwd=0.1) +
  theme_minimal()+
  labs(x="Samples", y="Relative Abundance (%)", title=expression(paste('Staked barplot showing the relative abundane of the bacteria (Phylum;Class) in the different samples of ', italic("Haslea ostrearia"), ' identified by the 16S region V3V4')) ,
       caption=paste("* Only Phylum with an abundance >=", manyX*100, "% of the sample are representated")) + 
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values=cols_ggplot, name="")+
  theme(legend.position="bottom", legend.text=element_text(size=14), legend.title=element_text(size=14), axis.text.x=element_text(size=14), axis.text.y=element_text(face="bold", size=14))+
  guides(fill=guide_legend(ncol=2, reverse=TRUE))
ggplot_V3V4


# BUBBLE CHART
## Input
V3V4=read_delim("V3V4_taxons_uncultured.csv", 
                   delim="\t", escape_double=FALSE, 
                   col_names=TRUE, col_types=cols(Taxa=col_character(), N1=col_number(), 
                                                      N2=col_number(), T1=col_number(), 
                                                      T2=col_number()), 
                   trim_ws=TRUE)
## Melt
V3V4$Taxa=str_replace(V3V4$Taxa, ";D_4__.*", "")
V3V4$Taxa=str_replace(V3V4$Taxa, ".*;D_2__", "")
V3V4$Taxa=str_replace(V3V4$Taxa, "D_3__", "")
V3V4$Taxa[grepl("Bacteria", V3V4$Taxa)]=c("Z Others")
V3V4_agg=aggregate(cbind(V3V4$N1, V3V4$N2, V3V4$T1, V3V4$T2), by=list(Taxa=V3V4$Taxa), FUN=sum)
colnames(V3V4_agg)=c("Taxa","N1","N2","T1","T2")
melt_V3V4=melt(V3V4_agg[apply(V3V4_agg[, 2:5], 1, any),], id.vars="Taxa", variable.name="Sample")
percent_melt_V3V4=get_percent_col(melt_V3V4)

## Plot
bubble_V3V4=ggplot(percent_melt_V3V4, aes(x=Sample, y=Taxa, size=percent, color=value)) +
  geom_point(alpha=0.8) +
  scale_size(range=c(.1, 24), name="Relative Abundance (%)") + 
  labs(y="Class;Order", x="Samples", title="Clade found by the analysis of V3V4 region", color="Mean Depth") +
  scale_y_discrete(limits=rev)+
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12), axis.text.x=element_text(face="bold",size=14), axis.text.y=element_text(face="italic",size=14))
bubble_V3V4


# KHI DEUX TEST
chisq.test(V3V4_agg$T1, V3V4_agg$T2)


# MOCK BARPLOT
## DATA
V1V3=read_delim("HOB9-16S-V1V3-Mock-June2022_contigsNR.meandepth.taxa.tsv", 
                   delim="\t", col_names=FALSE, col_types=cols(X2=col_character(), X1=col_number()))
V3V4=read_delim("HOB9-16S-V3V4-Mock-June2022_contigsNR.meandepth.taxa.tsv", 
                   delim="\t", col_names=FALSE, col_types=cols(X2=col_character(), X1=col_number()))

## V1V3 melt
colnames(V1V3)=c("MeanDepth","Taxa")
V1V3$Taxa=str_replace(V1V3$Taxa, ";D_5__.*", "")
V1V3$Taxa=str_replace(V1V3$Taxa, ".*;D_3__", "")
V1V3$Taxa=str_replace(V1V3$Taxa, "D_4__", "")
V1V3$Taxa[grepl("Unassigned", V1V3$Taxa)]=c("Z Unassigned")
V1V3_agg=aggregate(V1V3$MeanDepth, by=list(Taxa=V1V3$Taxa), FUN=sum)
colnames(V1V3_agg)=c("Taxa","Mock_V1V3")
melt_V1V3=melt(V1V3_agg, id.vars="Taxa", variable.name="Sample")
percent_melt_V1V3=get_percent_col(melt_V1V3)
## V3V4 melt
colnames(V3V4)=c("MeanDepth","Taxa")
V3V4$Taxa=str_replace(V3V4$Taxa, ";D_5__.*", "")
V3V4$Taxa=str_replace(V3V4$Taxa, ".*;D_3__", "")
V3V4$Taxa=str_replace(V3V4$Taxa, "D_4__", "")
V3V4$Taxa[grepl("Unassigned", V3V4$Taxa)]=c("Z Unassigned")
V3V4_agg=aggregate(V3V4$MeanDepth, by=list(Taxa=V3V4$Taxa), FUN=sum)
colnames(V3V4_agg)=c("Taxa","Mock_V3V4")
melt_V3V4=melt(V3V4_agg, id.vars="Taxa", variable.name="Sample")
percent_melt_V3V4=get_percent_col(melt_V3V4)
## ALL melt
percent_all=rbind(percent_melt_V1V3,percent_melt_V3V4)
percent_all$Taxa=str_replace(percent_all$Taxa, ".*;", "")
## Barplot
barplot_all=ggplot(data=percent_all, aes(x=Taxa, y=percent, fill=Sample)) + 
  geom_bar(position="dodge",stat="identity", na.rm=TRUE, colour="darkgrey", lwd=0.1) +
  theme_minimal()+
  scale_y_continuous(labels=scales::percent_format(scale=1))+
  labs(x="Order or Family", y="Relative Abundance (%)", caption=paste("* Chloroplast is the only order which wasn't associated to any family")) + 
  scale_fill_manual(values=cols_ggplot, name="")+
  theme(legend.position="top", legend.text=element_text(size=16), axis.text.x=element_text(size=18, angle=20, hjust=1), axis.text.y=element_text(size=18))+
  geom_hline(yintercept=16.5, linetype="dashed", color="darkgrey", size=1)
barplot_all


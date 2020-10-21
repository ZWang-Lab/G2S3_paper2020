args = commandArgs(TRUE)
files_id= args[1]
dat_path<-list.files("./Data/")[as.numeric(files_id)]
method<-gsub(".*[_]([^.]+)[.].*", "\\1", dat_path)
library(Seurat)
library(monocle)
library(TSCAN)
library(SingleCellExperiment)
#Input gene expression matrix
dat<-readRDS(paste0("/gpfs/loomis/scratch60/xiting_yan/yl883/081020/Traj/Data/",dat_path))
dat<-UpdateSeuratObject(dat)

res_mat<-dat@assays$RNA@counts
print(dim(res_mat))
tot_counts <- Matrix::colSums(res_mat)
summary(tot_counts)
sce <- SingleCellExperiment(assays = list(counts = res_mat))
sce$cell_type <-dat$Cell_stage

pd <- data.frame(cell_id = colnames(dat), cell_type = dat$Cell_stage,row.names = colnames(dat))
pd <- new("AnnotatedDataFrame", data = pd)
fd <- data.frame(gene_id = rownames(dat), gene_short_name = rownames(dat),row.names = row.names(dat))
fd <- new("AnnotatedDataFrame", data = fd)


if(min(counts(sce))<0){
  cds <- newCellDataSet(counts(sce), phenoData = pd, featureData = fd,expressionFamily = uninormal())
  cds <- estimateSizeFactors(cds)
  ordering_genes <-data.frame(gene_id=row.names(counts(sce)),mean_expression=Matrix::rowMeans(counts(sce)))
}else{
  cds <- newCellDataSet(counts(sce), phenoData = pd, featureData = fd)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  disp_table <- dispersionTable(cds)
  ordering_genes <-subset(disp_table, mean_expression >= 0.1)
}
print(dim(ordering_genes))
print(str(ordering_genes))

cds_orig<-cds
od_gene_orig<-ordering_genes
ordering_genes<-od_gene_orig$gene_id[order(od_gene_orig$mean_expression,decreasing = T)][1:16935]

#cds <- setOrderingFilter(cds, ordering_genes$gene_id)
cds <- setOrderingFilter(cds, ordering_genes)
if(min(counts(sce))<0){
  cds <- reduceDimension(cds,norm_method = "none")
}else{
  cds <- reduceDimension(cds)
}
cds <- orderCells(cds)
saveRDS(cds,paste0("./monocle2_",method,"_full_order.rds"))
pseudo<-cds$Pseudotime
cellLabels=factor(dat$Cell_stage,levels = c("E3","E4","E5","E6","E7" ))
cor.kendall = cor(pseudo, as.numeric(cellLabels), method = "kendall", use = "complete.obs")
subpopulation <- data.frame(cell = colnames(dat), sub = as.numeric(cellLabels)-1)
pseudo.df<-data.frame(sample_name=c(1:length(pseudo)),State=as.numeric(cellLabels),Pseudotime=order(unlist(pseudo)))
POS <- orderscore(subpopulation, pseudo.df)
print(POS)
print(cor.kendall)



library(monocle)
library(ggplot2)
library(gridExtra)
library(grid)
res_path<-"~/Desktop/081020/Traj/"
res_path<-"~/yl883/scratch60/081020/Traj/"
methods<-c("alra","g2s3","knnsmooth","raw","sctssr","saver","saucie","enimpute","magic","scimpute","dca","viper")
method_title<-c("ALRA","G2S3","kNN-smoothing","Raw","scTSSR","SAVER","SAUCIE","EnImpute","MAGIC","scImpute","DCA","VIPER")
dat_file<-readRDS(paste0(res_path,"/Data/Petro_seurat_raw.rds"))
cellLabels=factor(dat_file$Cell_stage,levels = c("E3","E4","E5","E6","E7" ))
traj_score<-matrix(NA,nrow=length(methods),ncol=2)
figure_traj<-list()

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
for(i in 1:12){
  file_name<-paste0(res_path,"monocle2_",methods[i],"_full_order.rds")
  cds_file<-readRDS(file_name)
  pseudo<-cds_file$Pseudotime
  cor.kendall = cor(pseudo, as.numeric(cellLabels), method = "kendall", use = "complete.obs")
  subpopulation <- data.frame(cell = colnames(dat_file), sub = as.numeric(cellLabels)-1)
  pseudo.df<-data.frame(sample_name=c(1:length(pseudo)),State=as.numeric(cellLabels),Pseudotime=order(unlist(pseudo)))
  POS <- orderscore(subpopulation, pseudo.df)
  traj_score[i,]<-c(POS[3],cor.kendall)
  if(i%in%c(4,12)){
    figure_traj[[i]]<-plot_cell_trajectory(cds_file, color_by = "cell_type", cell_size = 1) +
      ggtitle(paste0(method_title[i],"\n","(POS: ",round(abs(POS[3]),2),", Cor: ",round(abs(cor.kendall),2),")"))+
      theme(plot.title = element_text(hjust = 0.5),legend.position = "right",legend.justification="center",
            legend.margin=margin(0,0,0,-100),
            legend.box.margin=margin(-5,-5,-5,-5),legend.text = element_text(size=14),legend.title = element_blank())+ guides(colour = guide_legend(override.aes = list(size=5)))
    legend <- get_legend(figure_traj[[i]])
    figure_traj[[i]]<-plot_cell_trajectory(cds_file, color_by = "cell_type", cell_size = 1,show_branch_points = F) +
      ggtitle(method_title[i],subtitle = paste0("(POS: ",round(abs(POS[3]),2),", Cor: ",round(abs(cor.kendall),2),")"))+
      theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"),plot.subtitle= element_text(hjust = 0.5,size=14),legend.position = "none",
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12))+
      xlab("Dim 1") + ylab("Dim 2")
  }else{
    figure_traj[[i]]<-plot_cell_trajectory(cds_file, color_by = "cell_type", cell_size = 1,show_branch_points = F) +
      ggtitle(method_title[i],subtitle = paste0("(POS: ",round(abs(POS[3]),2),", Cor: ",round(abs(cor.kendall),2),")"))+
      theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"),plot.subtitle= element_text(hjust = 0.5,size=14),legend.position = "none",axis.title.x = element_text(size=12),axis.title.y = element_blank())+
      xlab("Dim 1")
  }

}
colnames(traj_score)<-c("POS","Kendall")
figure_traj[[13]]<-legend
tiff("~/Desktop/Trajectory_MethodsComparison.tif", width = 18, height = 6, units = 'in', res = 300)
grid.arrange(grobs=figure_traj[c(4,2,6,3,9,10,1,5,11,7,8,12)],
             layout_matrix = cbind(c(13,1,1,13),rbind(c(2:6,12),c(2:6,12),c(7:11,12), c(7:11,12))))
dev.off()

pdf("~/Desktop/Trajectory_MethodsComparison.pdf", width = 18, height = 6)
grid.arrange(grobs=figure_traj[c(4,2,6,3,9,10,12,1,5,11,7,8,13)],
             layout_matrix = rbind(c(1:6,13),c(7:12,13)))
dev.off()

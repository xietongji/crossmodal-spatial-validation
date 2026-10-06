# crossmodal-spatial-validation
# 
# R scripts for cross-modal spatial validation of tumor differentiation
# signatures using serial-section spatial transcriptomics (ST) and
# multiplex immunofluorescence (mIF).
# 
# ## Two methodological contributions
# 
# 1. **Dual-perspective neighborhood analysis**
# - Tumor-centered: infiltration count within 20 um
# - Immune-centered: infiltration depth + differentiation score of
# the nearest border tumor cell
# 
# 2. **Grid-based cross-modal validation**
# - Fixed-boundary normalization of ST and mIF to a shared grid
# - binname-based pairing avoids single-cell alignment
# - Multi-resolution scan (10-500 bins) tests robustness
# 
# ## Requirements
# 
# - R >= 4.2.1
# - Packages: randomForest, ggplot2, tidyverse, ...
# 
# ## Data
# 
# - `data/mIF.Rdata`: single-cell mIF data
# - `data/ST.Rdata`: ST pseudotime data
# - Coordinate ranges: mIF X [3125, 9625] um, Y [6250, 12750] um;
# ST array_col [0,127], array_row [0,77]

# The RData files used in this code have been deposited at Zenodo:
# https://doi.org/10.5281/zenodo.23190640
# (includes mIF cell matrix and ST spot matrix)

#library packages
library(scales)
library(data.table)
library(ggplot2)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(ggnewscale)
library(ggpubr)
library(randomForest)
#load data (cell matrix from mIF)
load(file="C:\\Users\\Lenovo\\Desktop\\Github\\mIF.Rdata")
#identification of CellType 
mIF$CellType=ifelse(mIF$Parent=="Annotation (Tumor)" & is.na(mIF$Classification),"Tumor",NA)
mIF$CellType=ifelse(mIF$Classification %in% c("Opal520"),"Macrophage",mIF$CellType)
mIF$CellType=ifelse(mIF$Classification %in% c("Opal570"),"NK",mIF$CellType)
mIF$CellType=ifelse(mIF$Classification %in% c("Opal620"),"CD8T_cell",mIF$CellType)
mIF$CellType=ifelse(mIF$Classification %in% c("Opal690"),"CD4T_cell",mIF$CellType)
mIF$CellType=ifelse(mIF$Classification %in% c("Opal780"),"B_cell",mIF$CellType)
mIF$CellType=ifelse(is.na(mIF$CellType),"Other",mIF$CellType)
#selection of trainingregion
trainingregion=data.frame(cbind(Xmin=c(4000,4000,5200,4650,7300,3600),
                                Xmax=c(5000,5000,6000,5350,8000,4000),
                                Ymin=c(9700,6700,11500,8100,7000,11000),
                                Ymax=c(10500,7500,12500,9600,8000,12000),
                                Type=c(1,0,1,1,0,0)))
trainingset=data.frame()
for (i in 1:nrow(trainingregion)) {
  trainingset=rbind(trainingset,cbind(mIF[(mIF$Centroid_Y<trainingregion[i,"Ymax"]) & 
                                            (mIF$Centroid_Y>trainingregion[i,"Ymin"]) &
                                            (mIF$Centroid_X<trainingregion[i,"Xmax"]) &
                                            (mIF$Centroid_X>trainingregion[i,"Xmin"]),],
                                      Type=trainingregion[i,"Type"]))
}

trainingset=trainingset[trainingset$CellType=="Tumor",]
trainingset=trainingset[,c(5:115,117)]
#training of randomForest model
set.seed(1)
result=randomForest(Type~.,data=trainingset,importance=TRUE)
#calculation of model score
Score=as.data.frame(predict(result,mIF))
mIF=cbind(mIF,
          Score=Score[rownames(mIF),])
mIF$Score=round(mIF$Score,12)
mIF[mIF$CellType!="Tumor","Score"]=NA

# ============================================================
# Innovation 1: Dual-perspective neighborhood analysis
# Core idea: Quantify tumor-immune spatial relationships from
# two reference frames:
#   (A) Tumor-centered: immune infiltration count
#   (B) Immune-centered: infiltration depth (distance to nearest
#       border tumor cell, plus that tumor cell's differentiation score)
# ============================================================

# --- Mode A: Tumor-centered (infiltration count) ---
# For each tumor cell, count immune cells within a 20-um radius.
radius=20
cellcount=0
for (cell in 1:nrow(mIF)) {
  cellcount=cellcount+1
  if(cellcount%%500==0){
    print(paste0(cellcount,"_",Sys.time()))
  }
  if(mIF[cell,"CellType"]=="Tumor"){
    xposition=mIF[cell,"Centroid_X"]
    yposition=mIF[cell,"Centroid_Y"]
    submIF=mIF[mIF$Centroid_X<(xposition+radius) &
                 mIF$Centroid_X>(xposition-radius) &
                 mIF$Centroid_Y<(yposition+radius) &
                 mIF$Centroid_Y>(yposition-radius),]
    submIF$distance=((submIF$Centroid_X-xposition)^2+(submIF$Centroid_Y-yposition)^2)^(0.5)
    submIF=submIF[submIF$distance<radius,]
    mIF[cell,"Totalcell"]=nrow(submIF)-1
    mIF[cell,"Nontumorcell"]=nrow(submIF[submIF$CellType!="Tumor",])
    mIF[cell,"Macrophage"]=nrow(submIF[submIF$CellType=="Macrophage",])
    mIF[cell,"NK"]=nrow(submIF[submIF$CellType=="NK",])
    mIF[cell,"CD8T_cell"]=nrow(submIF[submIF$CellType=="CD8T_cell",])
    mIF[cell,"CD4T_cell"]=nrow(submIF[submIF$CellType=="CD4T_cell",])
    mIF[cell,"B_cell"]=nrow(submIF[submIF$CellType=="B_cell",])
    mIF[cell,"Other"]=nrow(submIF[submIF$CellType=="Other",])
    if(mIF[cell,"Nontumorcell"]!=0){
      mIF[rownames(submIF[submIF$CellType!="Tumor",]),"AroundTumor"]="Yes"
    }
  } else {
    mIF[cell,"Totalcell"]=NA
    mIF[cell,"Nontumorcell"]=NA
    mIF[cell,"Macrophage"]=NA
    mIF[cell,"NK"]=NA
    mIF[cell,"CD8T_cell"]=NA
    mIF[cell,"CD4T_cell"]=NA
    mIF[cell,"B_cell"]=NA
    mIF[cell,"Other"]=NA
  }
}
mIF$index=(mIF$Nontumorcell)/(mIF$Totalcell)
mIF$Border=ifelse(mIF$CellType!="Tumor",NA,
                  ifelse(mIF$Totalcell==0,"Island",
                         ifelse(mIF$index==0,"Core","Border")))
mIF$Macrophage_percent=(mIF$Macrophage)/(mIF$Totalcell)
mIF$NK_percent=(mIF$NK)/(mIF$Totalcell)
mIF$CD8T_cell_percent=(mIF$CD8T_cell)/(mIF$Totalcell)
mIF$CD4T_cell_percent=(mIF$CD4T_cell)/(mIF$Totalcell)
mIF$B_cell_percent=(mIF$B_cell)/(mIF$Totalcell)
# --- Mode B: Immune-centered (infiltration depth) ---
# For each immune cell, find the nearest border tumor cell and
# record (1) the distance and (2) that tumor cell's differentiation score.
cellcount=0
for (cell in 1:nrow(mIF)) {
  cellcount=cellcount+1
  if(cellcount%%1000==0){
    print(paste0(cellcount,"_",Sys.time()))
  }
  if(!(mIF[cell,"CellType"] %in% c("Tumor","Other"))){
    xposition=mIF[cell,"Centroid_X"]
    yposition=mIF[cell,"Centroid_Y"]
    for (extend in 1:4) {
      newradius=10^(extend)
      submIF=mIF[mIF$Centroid_X<(xposition+newradius) &
                   mIF$Centroid_X>(xposition-newradius) &
                   mIF$Centroid_Y<(yposition+newradius) &
                   mIF$Centroid_Y>(yposition-newradius),]
      submIF=submIF[submIF$Border %in% c("Border"),]
      if(nrow(submIF)==0){next}
      submIF$distance=((submIF$Centroid_X-xposition)^2+(submIF$Centroid_Y-yposition)^2)^(0.5)
      submIF=submIF[submIF$distance<newradius,]
      if(nrow(submIF)==0){next}
      submIF=submIF[order(submIF$distance),]
      mIF[cell,"mindistance"]=submIF[1,"distance"]
      mIF[cell,"mindistanceBorderScore"]=submIF[1,"Score"]
      break
    }
  } else {
    mIF[cell,"mindistance"]=NA
    mIF[cell,"mindistanceBorderScore"]=NA
  }
}

rm(list=setdiff(ls(),c("mIF","trainingregion")))
gc()
#overlapped part
mIFSTpart=mIF[(mIF$Centroid_Y<12750) & 
                (mIF$Centroid_Y>6250) &
                (mIF$Centroid_X<9625) &
                (mIF$Centroid_X>3125),]

trainingset=data.frame()
for (i in 1:nrow(trainingregion)) {
  trainingset=rbind(trainingset,cbind(mIF[(mIF$Centroid_Y<trainingregion[i,"Ymax"]) & 
                                            (mIF$Centroid_Y>trainingregion[i,"Ymin"]) &
                                            (mIF$Centroid_X<trainingregion[i,"Xmax"]) &
                                            (mIF$Centroid_X>trainingregion[i,"Xmin"]),],
                                      Type=trainingregion[i,"Type"]))
}
trainingset=trainingset[,c("CellType","Score","Type")]
trainingset=trainingset[trainingset$CellType=="Tumor",]

mIFSTpartTraining=mIFSTpart[rownames(trainingset),]
mIFSTpartTest=mIFSTpart[setdiff(rownames(mIFSTpart),rownames(trainingset)),]
#load data (spot matrix from ST)
load("C:\\Users\\Lenovo\\Desktop\\Github\\ST.Rdata")
rm(list=setdiff(ls(),c("mIF","mIFSTpart","mIFSTpartTraining","mIFSTpartTest","trainingregion","spatialpathwaydata")))
gc() 

# ============================================================
# Innovation 2: Grid-based cross-modal validation
# Core idea: Normalize ST and mIF to a shared grid coordinate
# system using fixed physical boundaries, then pair bins by
# binname string intersection. No single-cell alignment needed.
# ============================================================

# --- Step 1: Build virtual cell grid for ST (500 x 500) ---
# Each virtual cell inherits the pseudotime of its nearest ST spot
# within a radius of 5 grid units.
spatialpathwaydata=na.omit(spatialpathwaydata)
spatialpathwaydata$binx=spatialpathwaydata$array_col/127
spatialpathwaydata$biny=spatialpathwaydata$array_row/77
spatialpathwaydata$binx=spatialpathwaydata$binx*500
spatialpathwaydata$biny=spatialpathwaydata$biny*500

newST=data.frame()
cellcount=0
radius=5
for (pseudocellx in 0:500) {
  for (pseudocelly in 0:500) {
    cellcount=cellcount+1
    if(cellcount%%10000==0){
      print(paste0(cellcount,"_",Sys.time()))
    }
    
    subspatialpathwaydata=spatialpathwaydata[spatialpathwaydata$binx<(pseudocellx+radius) &
                                               spatialpathwaydata$binx>(pseudocellx-radius) &
                                               spatialpathwaydata$biny<(pseudocelly+radius) &
                                               spatialpathwaydata$biny>(pseudocelly-radius),]
    subspatialpathwaydata$distance=((subspatialpathwaydata$binx-pseudocellx)^2+(subspatialpathwaydata$biny-pseudocelly)^2)^(0.5)
    subspatialpathwaydata=subspatialpathwaydata[subspatialpathwaydata$distance<radius,]
    
    if(nrow(subspatialpathwaydata)==0){
      newST=rbind(newST,cbind(pseudocellx=pseudocellx,
                              pseudocelly=pseudocelly,
                              Score=NA))
    } else {
      subspatialpathwaydata=subspatialpathwaydata[order(subspatialpathwaydata$distance),]
      newST=rbind(newST,cbind(pseudocellx=pseudocellx,
                              pseudocelly=pseudocelly,
                              Score=subspatialpathwaydata[1,"PseudotimeF"]))
    }
  }
}

newST=na.omit(newST)
rm(list=setdiff(ls(),c("mIF","mIFSTpart","mIFSTpartTraining","mIFSTpartTest","trainingregion","spatialpathwaydata","newST")))
save.image(file="C:\\Users\\Lenovo\\Desktop\\Github\\STmIF.RData")
gc()

# --- Step 2: Multi-resolution scan (10 to 500, step 10) ---
binresult=data.frame()
plotbin=c(10,20,50,100,500,1000)
for (bin in (c(1:150)*10)) {
  #show procession
  print(paste0(bin,"_",Sys.time()))
  #load data
  load(file="C:\\Users\\Lenovo\\Desktop\\Github\\STmIF.RData")
  # mIF side: fixed physical boundaries -> bin index
  # mIF coordinate range: X [3125, 9625] um, Y [6250, 12750] um
  #Training part of mIF (ST overlapped)
  mIFSTpartTraining$binx=(mIFSTpartTraining$Centroid_X-3125)/(9625-3125)
  mIFSTpartTraining$binx=mIFSTpartTraining$binx*bin
  mIFSTpartTraining$binx=as.integer(mIFSTpartTraining$binx)
  mIFSTpartTraining$binx=mIFSTpartTraining$binx+1
  mIFSTpartTraining$binx=ifelse(mIFSTpartTraining$binx==(bin+1),
                                bin,mIFSTpartTraining$binx)
  
  mIFSTpartTraining$biny=(mIFSTpartTraining$Centroid_Y-6250)/(12750-6250)
  mIFSTpartTraining$biny=mIFSTpartTraining$biny*bin
  mIFSTpartTraining$biny=as.integer(mIFSTpartTraining$biny)
  mIFSTpartTraining$biny=mIFSTpartTraining$biny+1
  mIFSTpartTraining$biny=ifelse(mIFSTpartTraining$biny==(bin+1),
                                bin,mIFSTpartTraining$biny)
  mIFSTpartTraining=mIFSTpartTraining[mIFSTpartTraining$CellType=="Tumor",]
  
  mIFSTpartTraining=aggregate(mIFSTpartTraining$Score,by=list(mIFSTpartTraining$binx,mIFSTpartTraining$biny),"mean")
  
  mIFSTpartTraining$binname=paste0("X",(bin-mIFSTpartTraining$Group.2),"_",
                                   "Y",(bin-mIFSTpartTraining$Group.1))
  mIFSTpartTraining$binxnew=bin-mIFSTpartTraining$Group.2
  mIFSTpartTraining$binynew=bin-mIFSTpartTraining$Group.1
  mIFSTpartTraining=mIFSTpartTraining[,c("binxnew","binynew","x","binname")]
  colnames(mIFSTpartTraining)=c("binx","biny","Score","binname")
  #Test part of mIF (ST overlapped)
  mIFSTpartTest$binx=(mIFSTpartTest$Centroid_X-3125)/(9625-3125)
  mIFSTpartTest$binx=mIFSTpartTest$binx*bin
  mIFSTpartTest$binx=as.integer(mIFSTpartTest$binx)
  mIFSTpartTest$binx=mIFSTpartTest$binx+1
  mIFSTpartTest$binx=ifelse(mIFSTpartTest$binx==(bin+1),
                            bin,mIFSTpartTest$binx)
  
  mIFSTpartTest$biny=(mIFSTpartTest$Centroid_Y-6250)/(12750-6250)
  mIFSTpartTest$biny=mIFSTpartTest$biny*bin
  mIFSTpartTest$biny=as.integer(mIFSTpartTest$biny)
  mIFSTpartTest$biny=mIFSTpartTest$biny+1
  mIFSTpartTest$biny=ifelse(mIFSTpartTest$biny==(bin+1),
                            bin,mIFSTpartTest$biny)
  mIFSTpartTest=mIFSTpartTest[mIFSTpartTest$CellType=="Tumor",]
  
  mIFSTpartTest=aggregate(mIFSTpartTest$Score,by=list(mIFSTpartTest$binx,mIFSTpartTest$biny),"mean")
  
  mIFSTpartTest$binname=paste0("X",(bin-mIFSTpartTest$Group.2),"_",
                               "Y",(bin-mIFSTpartTest$Group.1))
  mIFSTpartTest$binxnew=bin-mIFSTpartTest$Group.2
  mIFSTpartTest$binynew=bin-mIFSTpartTest$Group.1
  mIFSTpartTest=mIFSTpartTest[,c("binxnew","binynew","x","binname")]
  colnames(mIFSTpartTest)=c("binx","biny","Score","binname")
  #All of mIF (ST overlapped)
  mIFSTpart$binx=(mIFSTpart$Centroid_X-3125)/(9625-3125)
  mIFSTpart$binx=mIFSTpart$binx*bin
  mIFSTpart$binx=as.integer(mIFSTpart$binx)
  mIFSTpart$binx=mIFSTpart$binx+1
  mIFSTpart$binx=ifelse(mIFSTpart$binx==(bin+1),
                        bin,mIFSTpart$binx)
  
  mIFSTpart$biny=(mIFSTpart$Centroid_Y-6250)/(12750-6250)
  mIFSTpart$biny=mIFSTpart$biny*bin
  mIFSTpart$biny=as.integer(mIFSTpart$biny)
  mIFSTpart$biny=mIFSTpart$biny+1
  mIFSTpart$biny=ifelse(mIFSTpart$biny==(bin+1),
                        bin,mIFSTpart$biny)
  mIFSTpart=mIFSTpart[mIFSTpart$CellType=="Tumor",]
  
  mIFSTpart=aggregate(mIFSTpart$Score,by=list(mIFSTpart$binx,mIFSTpart$biny),"mean")
  # Coordinate flip to correct for serial-section scanning direction
  mIFSTpart$binname=paste0("X",(bin-mIFSTpart$Group.2),"_",
                           "Y",(bin-mIFSTpart$Group.1))
  mIFSTpart$binxnew=bin-mIFSTpart$Group.2
  mIFSTpart$binynew=bin-mIFSTpart$Group.1
  mIFSTpart=mIFSTpart[,c("binxnew","binynew","x","binname")]
  colnames(mIFSTpart)=c("binx","biny","Score","binname")
  
  # ST side: virtual cell grid -> bin index (fixed boundary 0-500)
  newST$binx=newST$pseudocellx/500
  newST$binx=newST$binx*bin
  newST$binx=as.integer(newST$binx)
  newST$binx=newST$binx+1
  newST$binx=ifelse(newST$binx==(bin+1),
                    bin,newST$binx)
  
  newST$biny=newST$pseudocelly/500
  newST$biny=newST$biny*bin
  newST$biny=as.integer(newST$biny)
  newST$biny=newST$biny+1
  newST$biny=ifelse(newST$biny==(bin+1),
                    bin,newST$biny)
  
  newST=aggregate(newST$Score,by=list(newST$binx,newST$biny),"mean")
  
  newST$binname=paste0("X",newST$Group.1,"_","Y",newST$Group.2)
  colnames(newST)=c("binx","biny","Score","binname")
  newST$Score=(newST$Score-min(newST$Score))/(max(newST$Score)-min(newST$Score))
  #Plot(before Intersection)
  if(bin %in% plotbin){
    pdf(paste0("C:\\Users\\Lenovo\\Desktop\\Github\\",
               bin,"_All_BI.pdf"),width=16,height=8)
    plotdata=rbind(cbind(mIFSTpart,Group="mIF"),
                   cbind(newST,Group="ST"))
    plot=ggplot()+
      geom_tile(data=plotdata,
                aes(x=binx,y=biny,fill=Score))+
      scale_fill_gradientn(name="Score",
                           values=seq(0,1,0.2),
                           colors=c("#C7E9B4",
                                    "#60C2C0",
                                    "#1D91C0",
                                    "#24499E",
                                    "#081D58"))+
      facet_grid(.~Group)+
      ylab("")+xlab("")+xlim(0,bin)+ylim(0,bin)+
      coord_fixed()+
      theme_bw()+
      theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            panel.grid=element_blank())
    print(plot)
    dev.off()
  }
  # --- Step 3: Pair bins by binname intersection ---
  intesctbinTraining=intersect(newST$binname,mIFSTpartTraining$binname)
  intesctbinTest=intersect(newST$binname,mIFSTpartTest$binname)
  intesctbin=intersect(newST$binname,mIFSTpart$binname)
  
  rownames(newST)=newST$binname
  rownames(mIFSTpart)=mIFSTpart$binname
  rownames(mIFSTpartTraining)=mIFSTpartTraining$binname
  rownames(mIFSTpartTest)=mIFSTpartTest$binname
  # --- Step 4: Spearman correlation ---
  corresultTraining=cor.test(as.numeric(newST[intesctbinTraining,"Score"]),
                             as.numeric(mIFSTpartTraining[intesctbinTraining,"Score"]),
                             method="spearman")
  corresultTest=cor.test(as.numeric(newST[intesctbinTest,"Score"]),
                         as.numeric(mIFSTpartTest[intesctbinTest,"Score"]),
                         method="spearman")
  corresult=cor.test(as.numeric(newST[intesctbin,"Score"]),
                     as.numeric(mIFSTpart[intesctbin,"Score"]),
                     method="spearman")
  #Plot(after Intersection)
  if(bin %in% plotbin){
    #All
    pdf(paste0("C:\\Users\\Lenovo\\Desktop\\Github\\",
               bin,"_All_AI.pdf"),width=16,height=8)
    plotdata=rbind(cbind(mIFSTpart[intesctbin,],Group="mIF"),
                   cbind(newST[intesctbin,],Group="ST"))
    plot=ggplot()+
      geom_tile(data=plotdata,
                aes(x=binx,y=biny,fill=Score))+
      scale_fill_gradientn(name="Score",
                           values=seq(0,1,0.2),
                           colors=c("#C7E9B4",
                                    "#60C2C0",
                                    "#1D91C0",
                                    "#24499E",
                                    "#081D58"))+
      facet_grid(.~Group)+
      ylab("")+xlab("")+xlim(0,bin)+ylim(0,bin)+
      coord_fixed()+
      theme_bw()+
      theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            panel.grid=element_blank())
    print(plot)
    dev.off()
    #Training
    pdf(paste0("C:\\Users\\Lenovo\\Desktop\\Github\\",
               bin,"_Training_AI.pdf"),width=16,height=8)
    plotdata=rbind(cbind(mIFSTpartTraining[intesctbinTraining,],Group="mIF"),
                   cbind(newST[intesctbinTraining,],Group="ST"))
    plot=ggplot()+
      geom_tile(data=plotdata,
                aes(x=binx,y=biny,fill=Score))+
      scale_fill_gradientn(name="Score",
                           values=seq(0,1,0.2),
                           colors=c("#C7E9B4",
                                    "#60C2C0",
                                    "#1D91C0",
                                    "#24499E",
                                    "#081D58"))+
      facet_grid(.~Group)+
      ylab("")+xlab("")+xlim(0,bin)+ylim(0,bin)+
      coord_fixed()+
      theme_bw()+
      theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            panel.grid=element_blank())
    print(plot)
    dev.off()
    #Test
    pdf(paste0("C:\\Users\\Lenovo\\Desktop\\Github\\",
               bin,"_Test_AI.pdf"),width=16,height=8)
    plotdata=rbind(cbind(mIFSTpartTest[intesctbinTest,],Group="mIF"),
                   cbind(newST[intesctbinTest,],Group="ST"))
    plot=ggplot()+
      geom_tile(data=plotdata,
                aes(x=binx,y=biny,fill=Score))+
      scale_fill_gradientn(name="Score",
                           values=seq(0,1,0.2),
                           colors=c("#C7E9B4",
                                    "#60C2C0",
                                    "#1D91C0",
                                    "#24499E",
                                    "#081D58"))+
      facet_grid(.~Group)+
      ylab("")+xlab("")+xlim(0,bin)+ylim(0,bin)+
      coord_fixed()+
      theme_bw()+
      theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            panel.grid=element_blank())
    print(plot)
    dev.off()
  }
  #Results
  binresult=rbind(binresult,cbind(bin=bin,
                                  rTraining=corresultTraining$estimate,
                                  pTraining=corresultTraining$p.value,
                                  numTraining=length(intesctbinTraining),
                                  rTest=corresultTest$estimate,
                                  pTest=corresultTest$p.value,
                                  numTest=length(intesctbinTest),
                                  r=corresult$estimate,
                                  p=corresult$p.value,
                                  num=length(intesctbin),
                                  STnum=nrow(newST),
                                  mIFnumTraining=nrow(mIFSTpartTraining),
                                  mIFnumTest=nrow(mIFSTpartTest),
                                  mIFnum=nrow(mIFSTpart)))
  rm(list=setdiff(ls(),c("bin","binresult","plotbin")))
  gc() 
}

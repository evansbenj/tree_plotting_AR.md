```
setwd('/projects/2018_Hymenochirus_Pipa_mello/mtDNA_genomez')
library(ape)
# must export as nexus tree with FigTree before reading
tree<-read.nexus("All_mtDNA_genomez_new_and_Gb_align5_allpipids_NODLOOP.nexus_boot.con.nex.tre")
DMW_data<-read.csv("DMW_state.txt")
DMW_matrix<-as.matrix(DMW_data)
DMWW<-setNames(DMW_data[,1],rownames(DMW_data))
#DMWW
# prune rhino so we have a rooted tree
rooted_tree<-drop.tip(tree, 1, trim.internal = TRUE, root.edge = 0, rooted = is.rooted(tree))
#plot(rooted_tree, layout=1)
#tree$tip.label
# this roots the tree using the first taxon (1) as the outgroup
#rooted_tree <- root(tree, 1, resolve.root = TRUE )
#DMW<-c(0,0,0,0,0,NA,0,2,2,2,2,1,1,1,2,1,1,NA,1,2,2,2,2,2,1,2,1,2,2,2,2)
#names(DMW)<-rooted_tree$tip.label
# list by columns with zeros in the diagonal
transitions <- matrix (c(0, 0, 0, 1, 0, 2, 1, 2, 0), nrow=3)
# this transition matrix has a separate rate between fem-specific and non-specific
#transitions_separate <- matrix (c(0, 0, 0, 1, 0, 3, 1, 2, 0), nrow=3)
# This allows DMW to orignate once to female specific or not female specific state
# and then for transitions between female specific and notfemale specific to occur 
# with a different rate from the rate that DMW originates. Once it originates, an
# absence of DMW is coded as not female specific
AR <- ace(DMW_data$DMW_NA, rooted_tree, type="discrete", model=transitions)
#AR_separate <- ace(DMW, rooted_tree, type="discrete", model=transitions_separate)
#plotTree(rooted_tree,type="phylogram",fsize=0.4,ftype="i",lwd=1)
#cols<-setNames(c("black","red","blue"),levels(as.factor(DMWW)))

#DMWWW<-as.factor(DMWW)
#tiplabels(pie=as.matrix(DMW_data$species[rooted_tree$tip.label], levels(DMWWW)),piecol=cols,cex=0.3)


#plotTree(rooted_tree,type="phylogram",fsize=0.4,ftype="i",lwd=1, showTipLabel = FALSE)
#nodelabels(node=1:rooted_tree$Nnode+Ntip(rooted_tree),
#           pie=AR_separate$lik.anc,piecol=cols,cex=0.4)
#tiplabels(tip=1:31,pch = 19,col=DMW_data$color,bg=DMW_data$color,cex=0.8)
#edgelabels(DMW_data$species,adj = c(-0.8),frame = "none",
#           pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
#           col = "black", offset=3, cex = 0.6)

cols<-c("black","red","light blue","blue","pink","white")
names<-DMW_data$name

# GOOD
pdf("./phylogeny.pdf",w=8, h=4.5, version="1.4", bg="transparent")
# plot tree with no labels
plot(rooted_tree, use.edge.length = TRUE, font = 1, label.offset = 3,
     x.lim = 1, no.margin = TRUE, show.tip.label = FALSE)
# add scale bar
add.scale.bar(x=0.03,y=20, length = 0.05, ask = FALSE, cex=0.7,
              lwd = 1, lcol = "black")
# plot the ancestral node reconstructions using pie charts
nodelabels(node=1:rooted_tree$Nnode+Ntip(rooted_tree),
           pie=AR$lik.anc,piecol=cols,cex=0.35,frame="none")
# add support values to individual nodes
nodelabels(text="73", node=32,adj = c(-0.6),bg=NULL,frame="none", cex=0.5)
nodelabels(text="71", node=38,adj = c(-0.6),bg=NULL,frame="none", cex=0.5)
nodelabels(text="83", node=43,adj = c(-0.6),bg=NULL,frame="none", cex=0.5)
nodelabels(text="92", node=52,adj = c(-0.6),bg=NULL,frame="none", cex=0.5)
nodelabels(text="92", node=52,adj = c(-0.6),bg=NULL,frame="none", cex=0.5)
nodelabels(text="49", node=47,adj = c(-0.6),bg=NULL,frame="none", cex=0.5)

# add symbols to tips
tiplabels(tip=1:30,pch = 21,bg=as.matrix(DMW_data$color),col="black",cex=1,lwd=1)
# add nice names to tips
tiplabels(as.vector(names), font = 3,cex=0.5,adj = c(0),tip=1:30,bg=NULL,frame="none",offset=0.01)
#add.simmap.legend(colors=cols,prompt=FALSE,vertical=TRUE,shape="circle",x=0.1,y=20,)
#add.simmap.legend(colors=sapply(setNames(cols,c("No DM-W","DM-W fixed only in females","DM-W not detected in one female","DM-W not fixed only in females","DM-W present in one female, but unknown specificity","DM-W not present in one male")),make.transparent
  #                              ),fsize=0.6,vertical=TRUE,shape="circle",cex=0.1,lwd=0.5,prompt=FALSE,x=0.03,y=28)
legend(
       "topleft", 
       legend=c("Never had DM-W","DM-W fixed only in females","DM-W not fixed only in females","DM-W not detected in one female","DM-W present in one female, but unknown sex-specificity","DM-W not present in one male"), 
       col="black",
       #lty = c(NA,NA,NA,NA,NA),
       #pch = c(21,21,21,21,21),
       bty = "n", 
       fill=cols, 
       horiz=FALSE, 
       cex=0.6, 
       pt.cex = 1,
       inset=.02
      )
dev.off()

```

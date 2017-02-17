#websites referred https://sites.google.com/site/daishizuka/toolkits/sna/plotting-networks-pt-2

#preprocessing data using qiime

#packages required
require(igraph)
require(Hmisc)

###########################root 90% OTUs###################################

#this table contains both rhizosphere and endosphere otus that were clustered at 90% and then jacknifed
merge_table<-read.table("~/Dropbox/Co-occurence/otu_table/merge_bac_fung_90_v3_R.txt",header=TRUE,row.names=1)
#transforming the table for correlation 
merge_table<-t(merge_table)

#removing columns/otus with sum of 1 or singletons or others
#merge_table<-merge_table[,!(colSums(abs(merge_table)) == 1)]
merge_table<-merge_table[,!(colSums(abs(merge_table)) < 5)]


#correlation between each OTUs
cor.merge.otu<-rcorr(merge_table,type="spearman")

# correlation values spearman's
r<-cor.merge.otu$r

# p-values
p<-cor.merge.otu$P 

#subset OTUs that have p-values less than 0.5
cor.merge.otu<-ifelse(p < .05, r,0)

#assigns 0 to lower triangle and diagonal.
cor.merge.otu[lower.tri(cor.merge.otu, diag=TRUE) ]<- 0

#assigns 0 to correlation lower than 0.6
cor.merge.otu[ abs(cor.merge.otu) < 0.6]<- 0

# change to graph format
graph <- graph.adjacency(cor.bac.otu > 0.6, mode="upper",add.rownames="OTUs") 

#assigning spearman r value to Edge
E(graph)$weight<-t(cor.merge.otu)[abs(t(cor.merge.otu))>0.6] 
E(graph)[weight>0.6 & weight <0.8]$color <- "green" #assigning colors
E(graph)[weight>=0.8 ]$color <- "red"

# size of vertices based on total sequence in that OTU
V(graph)$size<-log2(colSums(merge_table)) 

# remove vertices with no connections
new.graph <- delete.vertices(graph, which(degree(graph) < 1))

#######assigning colors for fungal and bacterial otus
#reading the file with otu affiliation
bac_fung<-read.table("~/Dropbox/Co-occurence/others/bac_fungi_col.txt",header=TRUE)

#assigning the otus to bacfung variable
V(new.graph)$bacfung=as.character(bac_fung$taxa[match(V(new.graph)$name,bac_fung$OTU.ID)])

#reassigning bacfung to color variable
V(new.graph)$color=V(new.graph)$bacfung

# assigning actual color to bacteria and fungi
V(new.graph)$color=gsub("Bacteria","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("Fungi","red",V(new.graph)$color)

#plot the graph
plot(new.graph, layout=layout.fruchterman.reingold, vertex.label=NA )#structure of the graph (circle) 

##############################################################################

###########################rhizosphere 90% OTUs###################################
#this table contains rhizosphere otus
rhizo_table<-read.table("~/Dropbox/Co-occurence/otu_table/merge_bac_fung_90_v3_Rhizosphere_R.txt",header=TRUE,row.names=1)

#transforming the table for correlation 
rhizo_table<-t(rhizo_table)

#removing columns/otus with sum of 1 or singletons or others
#rhizo_table<-rhizo_table[,!(colSums(abs(rhizo_table)) == 1)]
rhizo_table<-rhizo_table[,!(colSums(abs(rhizo_table)) < 5)]


#correlation between each OTUs
cor.rhizo.otu<-rcorr(rhizo_table,type="spearman")

# correlation values spearman's
r<-cor.rhizo.otu$r

# p-values
p<-cor.rhizo.otu$P 

#subset OTUs that have p-values less than 0.5
cor.rhizo.otu<-ifelse(p < .05, r,0)

#assigns 0 to lower triangle and diagonal.
cor.rhizo.otu[lower.tri(cor.rhizo.otu, diag=TRUE) ]<- 0

#assigns 0 to correlation lower than 0.6
cor.rhizo.otu[ abs(cor.rhizo.otu) < 0.6]<- 0

# change to graph format
graph <- graph.adjacency(cor.rhizo.otu > 0.6, mode="upper",add.rownames="OTUs") 

#assigning spearman r value to Edge
E(graph)$weight<-t(cor.rhizo.otu)[abs(t(cor.rhizo.otu))>0.6] 
E(graph)[weight>0.6 & weight <0.8]$color <- "green" #assigning colors
E(graph)[weight>=0.8 ]$color <- "red"

# size of vertices based on total sequence in that OTU
V(graph)$size<-log2(colSums(rhizo_table)) 

# remove vertices with no connections
new.graph <- delete.vertices(graph, which(degree(graph) < 1))

#######assigning colors for fungal and bacterial otus
#reading the file with otu affiliation
bac_fung<-read.table("~/Dropbox/Co-occurence/others/bac_fungi_col.txt",header=TRUE)

#assigning the otus to bacfung variable
V(new.graph)$bacfung=as.character(bac_fung$taxa[match(V(new.graph)$name,bac_fung$OTU.ID)])

#reassigning bacfung to color variable
V(new.graph)$color=V(new.graph)$bacfung

# assigning actual color to bacteria and fungi
V(new.graph)$color=gsub("Bacteria","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("Fungi","red",V(new.graph)$color)

#plot the graph
plot(new.graph, layout=layout.fruchterman.reingold, vertex.label=NA )#structure of the graph (circle) 

##############################################################################


###########################endosphere 90% OTUs###################################
#this table contains endosphere otus
endo_table<-read.table("~/Dropbox/Co-occurence/otu_table/merge_bac_fung_90_v3_Endosphere_R.txt",header=TRUE,row.names=1)

#transforming the table for correlation 
endo_table<-t(endo_table)

#removing columns/otus with sum of 1 or singletons or others
endo_table<-endo_table[,!(colSums(abs(endo_table)) < 5)]

#correlation between each OTUs
cor.endo.otu<-rcorr(endo_table,type="spearman")

# correlation values spearman's
r<-cor.endo.otu$r

# p-values
p<-cor.endo.otu$P 

#subset OTUs that have p-values less than 0.5
cor.endo.otu<-ifelse(p < .05, r,0)

#assigns 0 to lower triangle and diagonal.
cor.endo.otu[lower.tri(cor.endo.otu, diag=TRUE) ]<- 0

#assigns 0 to correlation lower than 0.6
cor.endo.otu[ abs(cor.endo.otu) < 0.6]<- 0

# change to graph format
graph <- graph.adjacency(cor.endo.otu > 0.6, mode="upper",add.rownames="OTUs") 

#assigning spearman r value to Edge
E(graph)$weight<-t(cor.endo.otu)[abs(t(cor.endo.otu))>0.6] 
E(graph)[weight>0.6 & weight <0.8]$color <- "green" #assigning colors
E(graph)[weight>=0.8 ]$color <- "red"

# size of vertices based on total sequence in that OTU
V(graph)$size<-log2(colSums(endo_table)) 

# remove vertices with no connections
new.graph <- delete.vertices(graph, which(degree(graph) < 1))

#######assigning colors for fungal and bacterial otus
#reading the file with otu affiliation
bac_fung<-read.table("~/Dropbox/Co-occurence/others/bac_fungi_col.txt",header=TRUE)

#assigning the otus to bacfung variable
V(new.graph)$bacfung=as.character(bac_fung$taxa[match(V(new.graph)$name,bac_fung$OTU.ID)])

#reassigning bacfung to color variable
V(new.graph)$color=V(new.graph)$bacfung

# assigning actual color to bacteria and fungi
V(new.graph)$color=gsub("Bacteria","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("Fungi","red",V(new.graph)$color)

#plot the graph
plot(new.graph, layout=layout.fruchterman.reingold, vertex.label=NA )#structure of the graph (circle) 

##############################################################################


###########################root 97% OTUs###################################

#this table contains both rhizosphere and endosphere otus that were clustered at 90% and then jacknifed
merge_table<-read.table("~/Dropbox/Co-occurence/otu_table/merge_bac_fung_97_v3_R.txt",header=TRUE,row.names=1)
#transforming the table for correlation 
merge_table<-t(merge_table)

#removing columns/otus with sum of 1 or singletons or others
#merge_table<-merge_table[,!(colSums(abs(merge_table)) == 1)]
merge_table<-merge_table[,!(colSums(abs(merge_table)) < 20)]


#correlation between each OTUs
cor.merge.otu<-rcorr(merge_table,type="spearman")

# correlation values spearman's
r<-cor.merge.otu$r

# p-values
p<-cor.merge.otu$P 

#subset OTUs that have p-values less than 0.5
cor.merge.otu<-ifelse(p < .05, r,0)

#assigns 0 to lower triangle and diagonal.
cor.merge.otu[lower.tri(cor.merge.otu, diag=TRUE) ]<- 0

#assigns 0 to correlation lower than 0.6
cor.merge.otu[ abs(cor.merge.otu) < 0.6]<- 0

# change to graph format
graph <- graph.adjacency(cor.merge.otu > 0.6, mode="upper",add.rownames="OTUs") 

#assigning spearman r value to Edge
E(graph)$weight<-t(cor.merge.otu)[abs(t(cor.merge.otu))>0.6] 
E(graph)[weight>0.6 & weight <0.8]$color <- "green" #assigning colors
E(graph)[weight>=0.8 ]$color <- "red"

# size of vertices based on total sequence in that OTU
V(graph)$size<-log2(colSums(merge_table)) 

# remove vertices with no connections
new.graph <- delete.vertices(graph, which(degree(graph) < 1))

#######assigning colors for fungal and bacterial otus
#reading the file with otu affiliation for bacteria and fungi, 100,000+ are fungi
bac_fung<-read.table("~/Dropbox/Co-occurence/others/bac_fungi_col_97.txt",header=TRUE)

#assigning the otus to bacfung variable
V(new.graph)$bacfung=as.character(bac_fung$taxa[match(V(new.graph)$name,bac_fung$OTU.ID)])

#reassigning bacfung to color variable
V(new.graph)$color=V(new.graph)$bacfung

# assigning actual color to bacteria and fungi
V(new.graph)$color=gsub("Bacteria","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("Fungi","red",V(new.graph)$color)

#plot the graph
plot(new.graph, layout=layout.fruchterman.reingold, vertex.label=NA )#structure of the graph (circle) 


##############################################################################

###########################rhizosphere 97% OTUs###################################
#this table contains rhizosphere otus
rhizo_table<-read.table("~/Dropbox/Co-occurence/otu_table/merge_bac_fung_97_v3_Rhizosphere_R.txt",header=TRUE,row.names=1)

#transforming the table for correlation 
rhizo_table<-t(rhizo_table)

#removing columns/otus with sum of 1 or singletons or others
#rhizo_table<-rhizo_table[,!(colSums(abs(rhizo_table)) == 1)]
rhizo_table<-rhizo_table[,!(colSums(abs(rhizo_table)) < 5)]


#correlation between each OTUs
cor.rhizo.otu<-rcorr(rhizo_table,type="spearman")

# correlation values spearman's
r<-cor.rhizo.otu$r

# p-values
p<-cor.rhizo.otu$P 

#subset OTUs that have p-values less than 0.5
cor.rhizo.otu<-ifelse(p < .05, r,0)

#assigns 0 to lower triangle and diagonal.
cor.rhizo.otu[lower.tri(cor.rhizo.otu, diag=TRUE) ]<- 0

#assigns 0 to correlation lower than 0.6
cor.rhizo.otu[ abs(cor.rhizo.otu) < 0.6]<- 0

# change to graph format
graph <- graph.adjacency(cor.rhizo.otu > 0.6, mode="upper",add.rownames="OTUs") 

#assigning spearman r value to Edge
E(graph)$weight<-t(cor.rhizo.otu)[abs(t(cor.rhizo.otu))>0.6] 
E(graph)[weight>0.6 & weight <0.8]$color <- "green" #assigning colors
E(graph)[weight>=0.8 ]$color <- "red"

# size of vertices based on total sequence in that OTU
V(graph)$size<-log2(colSums(rhizo_table)) 

# remove vertices with no connections
new.graph <- delete.vertices(graph, which(degree(graph) < 1))

#######assigning colors for fungal and bacterial otus
#reading the file with otu affiliation
bac_fung<-read.table("~/Dropbox/Co-occurence/others/bac_fungi_col_97.txt",header=TRUE)

#assigning the otus to bacfung variable
V(new.graph)$bacfung=as.character(bac_fung$taxa[match(V(new.graph)$name,bac_fung$OTU.ID)])

#reassigning bacfung to color variable
V(new.graph)$color=V(new.graph)$bacfung

# assigning actual color to bacteria and fungi
V(new.graph)$color=gsub("Bacteria","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("Fungi","red",V(new.graph)$color)

#plot the graph
plot(new.graph, layout=layout.fruchterman.reingold, vertex.label=NA )#structure of the graph (circle) 

##############################################################################

###########################endosphere 97% OTUs###################################
#this table contains endosphere otus
endo_table<-read.table("~/Dropbox/Co-occurence/otu_table/merge_bac_fung_97_v3_Endosphere_R.txt",header=TRUE,row.names=1)

#transforming the table for correlation 
endo_table<-t(endo_table)

#removing columns/otus with sum of 1 or singletons or others
endo_table<-endo_table[,!(colSums(abs(endo_table)) < 5)]

#correlation between each OTUs
cor.endo.otu<-rcorr(endo_table,type="spearman")

# correlation values spearman's
r<-cor.endo.otu$r

# p-values
p<-cor.endo.otu$P 

#subset OTUs that have p-values less than 0.5
cor.endo.otu<-ifelse(p < .05, r,0)

#assigns 0 to lower triangle and diagonal.
cor.endo.otu[lower.tri(cor.endo.otu, diag=TRUE) ]<- 0

#assigns 0 to correlation lower than 0.6
cor.endo.otu[ abs(cor.endo.otu) < 0.6]<- 0

# change to graph format
graph <- graph.adjacency(cor.endo.otu > 0.6, mode="upper",add.rownames="OTUs") 

#assigning spearman r value to Edge
E(graph)$weight<-t(cor.endo.otu)[abs(t(cor.endo.otu))>0.6] 
E(graph)[weight>0.6 & weight <0.8]$color <- "green" #assigning colors
E(graph)[weight>=0.8 ]$color <- "red"

# size of vertices based on total sequence in that OTU
V(graph)$size<-log2(colSums(endo_table)) 

# remove vertices with no connections
new.graph <- delete.vertices(graph, which(degree(graph) < 1))

#######assigning colors for fungal and bacterial otus
#reading the file with otu affiliation
bac_fung<-read.table(~/Dropbox/Co-occurence/others/bac_fungi_col_97.txt",header=TRUE)

#assigning the otus to bacfung variable
V(new.graph)$bacfung=as.character(bac_fung$taxa[match(V(new.graph)$name,bac_fung$OTU.ID)])

#reassigning bacfung to color variable
V(new.graph)$color=V(new.graph)$bacfung

# assigning actual color to bacteria and fungi
V(new.graph)$color=gsub("Bacteria","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("Fungi","red",V(new.graph)$color)

#plot the graph
plot(new.graph, layout=layout.fruchterman.reingold, vertex.label=NA )#structure of the graph (circle) 











#rhizosphere otu table from 2010 data set
bac_table<-read.table("~/Desktop/PMI/PMI_Bacteria/FINAL_ANALYSIS_V1.5/IMPFILES/A.RHIZO_1000/final_otu_table_rhizo_1000_R.txt",header=TRUE,row.names=1)

#endosphere otu table from 2010 data set
bac_table<-read.table("~/Desktop/PMI/PMI_Bacteria/FINAL_ANALYSIS_V1.5/IMPFILES/B.ENDO_1000/final_otu_table_1000_Endosphere_R.txt",header=TRUE,row.names=1)

#reading the attribute files for endosphere taxa, remove row.names=1
endo.taxa<-read.table("~/Desktop/PMI/PMI_COOCCURENCE/endo_taxa.txt",header=TRUE)
bac_fung<-read.table("~/Desktop/bac_fungi_col.txt",header=TRUE)

###########################reading the table###################################

#removing columns with sum of 1 or singletons or others
bac_table<-bac_table[,!(colSums(abs(bac_table)) == 1)]
bac_table<-bac_table[,!(colSums(abs(bac_table)) < 5)]


#transforming the table
bac_table<-t(bac_table)

#converting to matrix
bac_table<-as.matrix(bac_table)

#correlation between each OTUs
cor.bac.otu<-rcorr(bac_table,type="spearman")

# correlation values
r<-cor.bac.otu$r

# p-values
p<-cor.bac.otu$P 

#subset OTUs that have p=values less than 0.5
cor.bac.otu<-ifelse(p < .05, r,0)

#assigns 0 to lower triangle and diagonal.
cor.bac.otu[lower.tri(cor.bac.otu, diag=TRUE) ]<- 0

#assigns 0 to correlation lower than 0.6
cor.bac.otu[ abs(cor.bac.otu) < 0.6]<- 0

# change to graph format
graph <- graph.adjacency(cor.bac.otu > 0.6, mode="upper",add.rownames="OTUs") 

#assigning spearman r value to E
E(graph)$weight<-t(cor.bac.otu)[abs(t(cor.bac.otu))>0.6] 

E(graph)[weight>0.6 & weight <0.8]$color <- "green" #assigning colors
E(graph)[weight>=0.8 ]$color <- "red"

#V(graph)$size<-(degree(graph)+1)*2 # size of the vertices based on # of connection

# size of vertices based on total sequence in that OTU
V(graph)$size<-log2(colSums(bac_table)) 

#assign OTU number as vertics name, order is maintained
#V(graph)$label <- colnames(cor.bac.otu)

# remove vertices with no connections
new.graph <- delete.vertices(graph, which(degree(graph) < 1))

#assign taxonomic information to vertices for endophyte
V(new.graph)$taxonomy=as.character(endo.taxa$taxonomy[match(V(new.graph)$name,endo.taxa$XOTU.ID)])
V(new.graph)$bacfung=as.character(bac_fung$taxa[match(V(new.graph)$name,bac_fung$OTU.ID)])



#assign taxonomy to color attribute
V(new.graph)$color=V(new.graph)$taxonomy
V(new.graph)$color=V(new.graph)$bacfung

#assigning colors to each category
V(new.graph)$color=gsub("Proteobacteria","red",V(new.graph)$color)
V(new.graph)$color=gsub("Bacteria","green",V(new.graph)$color)
V(new.graph)$color=gsub("Acidobacteria","black",V(new.graph)$color)
V(new.graph)$color=gsub("Actinobacteria","brown",V(new.graph)$color)
V(new.graph)$color=gsub("Chloroflexi","gray",V(new.graph)$color)
V(new.graph)$color=gsub("WS3","blue",V(new.graph)$color)
V(new.graph)$color=gsub("Verrucomicrobia","darkorange",V(new.graph)$color)
V(new.graph)$color=gsub("Planctomycetes","lightblue4",V(new.graph)$color)
V(new.graph)$color=gsub("Firmicutes","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("Gemmatimonadetes","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("TM7","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("Elusimicrobia","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("Bacteria","midnightblue",V(new.graph)$color)
V(new.graph)$color=gsub("Fungi","red",V(new.graph)$color)


#plot the graph
plot(new.graph, layout=layout.fruchterman.reingold, vertex.label=NA )#structure of the graph (circle) 
     
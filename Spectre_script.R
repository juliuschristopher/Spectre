## ## Spectre Workflow ####
# Example Data: Old Cyp11a1 Phenutyping 1 - 25/03/21

####Install Spectre from GitHub
if(!require('devtools')) {install.packages('devtools')}
library('devtools')

options(timeout=6000) # this changes the 'timeout' limit for downloading the package
devtools::install_github("immunedynamics/spectre") # this will download and install Spectre

####Load libraries####
library(Spectre)
Spectre::package.check() # Check that all required packages are installed
Spectre::package.load() # Load required packages

####Set working directory####
setwd("~/Desktop/Spectre")
PrimaryDirectory <- getwd()

####Set input directory####
setwd("~/Desktop/Spectre/data/") #where the data can be found
InputDirectory <- getwd()

####Set metadata directory####
setwd("~/Desktop/Spectre/metadata/")
MetaDirectory <- getwd()

####Set output directory####
setwd(PrimaryDirectory)
dir.create("Output_spectre", showWarnings = FALSE)
setwd("Output_spectre")
OutputDirectory <- getwd()
setwd(PrimaryDirectory) #set to primary directory!

####Import data####
##Load the data
setwd(InputDirectory)
list.files(InputDirectory, ".csv")

data.list <- Spectre::read.files(file.loc = InputDirectory,
                                 file.type = ".csv",
                                 do.embed.file.names = TRUE)
##Check the data
check <- do.list.summary(data.list)

check$name.table # Review column names and their subsequent values
check$ncol.check # Review number of columns (features, markers) in each sample
check$nrow.check # Review number of rows (cells) in each sample

data.list[1] #First six rows (cells)

##Merge data
cell.dat <- Spectre::do.merge.files(dat = data.list)
cell.dat

##Read in metadata
setwd(MetaDirectory)
meta.dat <- fread("sample.details.csv")
meta.dat

####Arcsinh transformation####
#Not required when CSV channal vaules have been exported from FlowJo

####Add sample metadata and set preferences####
##Select metadata columns to be added
setwd(MetaDirectory)

sample.info <- meta.dat[,c(1:4)] #Only first four columns will be added
sample.info

##Add metadata to cell.dat
cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "Filename", rmv.ext = TRUE)
cell.dat
as.matrix(names(cell.dat)) #Check columns

##Specifiy columns which represent cellular features
cellular.cols <- names(cell.dat)[c(1:15)]
as.matrix(cellular.cols)


##Specifiy columns for tSNE/UMAP clustering
cluster.cols <- names(cell.dat)[c(1:15)]
as.matrix(cluster.cols)

#Specifiy sample, group and batch columns
exp.name <- "Old Cyp11a1 Phenotyping (1)"
sample.col <- "Sample"
group.col <- "Group"
batch.col <- "Batch"


##Determine number for downsampling for dimensionality reduction
data.frame(table(cell.dat[[group.col]])) # Check number of cells per sample

sub.targets <- c(2000, 2000, 2000) # target subsample numbers from each group
sub.targets


####Clustering and dimensionality reduction####
setwd(OutputDirectory)
dir.create("Output - clustering")
setwd("Output - clustering")

##Clustering
cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = "auto") #auto clusters
fwrite(cell.dat, "clustered.data.csv")

##Dimensionality reduction
cell.sub <- do.subsample(cell.dat, sub.targets, group.col)
cell.sub <- run.umap(cell.sub, cluster.cols) #umap choosen

fwrite(cell.sub, "clustered.data.DR.csv")

##Visualise dimensionality reduction (DR) plots - metaclusters
make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", cellular.cols)
make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", group.col, col.type = 'factor')
make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_cluster", col.type = 'factor', plot.width = 20, plot.height = 20)
make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "Group", col.type = 'factor', plot.width = 10, plot.height = 10)
make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "Sample", col.type = 'factor', plot.width = 10, plot.height = 10)


##Expression heatmaps
exp <- do.aggregate(cell.dat, cellular.cols, by = "FlowSOM_metacluster")
make.pheatmap(exp, "FlowSOM_metacluster", cellular.cols)

exp1 <- do.aggregate(cell.dat, cellular.cols, by = "Group")
make.pheatmap(exp1, "Group", cellular.cols)

exp2 <- do.aggregate(cell.dat, cellular.cols, by = "Sample")
make.pheatmap(exp2, "Sample", cellular.cols)

exp3 <- do.aggregate(cell.dat, cellular.cols, by = "FlowSOM_cluster")
make.pheatmap(exp3, "FlowSOM_cluster", cellular.cols)

cell.sub %>%
  group_by(FlowSOM_metacluster) %>%
  summarise(n())

####Statistical analysis####
setwd(OutputDirectory)
dir.create("Output 4 - summary data")
setwd("Output 4 - summary data")

##Setup
variance.test <- 'kruskal.test'
pairwise.test <- "wilcox.test"

comparisons <- list(c("WT", "Mb1Cyp11a1KO", "Mb1_E1020KCyp11a1KO"))
comparisons

grp.order <- c("WT", "Mb1Cyp11a1KO", "Mb1E1020KCyp11a1KO")
grp.order

##Specifiy columns to measure MFI on
as.matrix(cellular.cols)
dyn.cols <- cellular.cols[c(5,13)] #CXCR4 and CD95
dyn.cols

##Create a summary table
sum.dat <- create.sumtable(dat = cell.dat,
                           sample.col = sample.col,
                           pop.col = "FlowSOM_metacluster",
                           use.cols = dyn.cols,
                           annot.cols = c(group.col, batch.col))
##Review summary data
sum.dat
as.matrix(names(sum.dat))

##Specifiy which columns to plot
plot.cols <- names(sum.dat)[c(4:15)]
plot.cols

##Re-order summary data
sum.dat <- do.reorder(sum.dat, group.col, grp.order)
sum.dat[,c(1:3)]

###Autographs
for(i in plot.cols){
  
  measure <- gsub("\\ --.*", "", i)
  measure
  
  pop <- gsub("^[^--]*.-- ", "", i)
  pop
  
  make.autograph(sum.dat,
                 x.axis = group.col,
                 y.axis = i,
                 y.axis.label = measure,
                 
                 grp.order = grp.order,
                 my_comparisons = comparisons,
                 
                 Variance_test = variance.test,
                 Pairwise_test = pairwise.test,
                 
                 title = pop,
                 subtitle = measure,
                 filename = paste0(i, '.pdf'))
}


####Create a fold change heat map####
##Z-score calculation
sum.dat.z <- do.zscore(sum.dat, plot.cols)

##Group
t.first <- match(grp.order, sum.dat.z[[group.col]])
t.first <- t.first -1
t.first

##Make heatmap
make.pheatmap(sum.dat.z,
              sample.col = sample.col,
              plot.cols = paste0(plot.cols, '_zscore'),
              is.fold = TRUE,
              plot.title = 'Z-score',
              dendrograms = 'column',
              row.sep = t.first,
              cutree_cols = 3)

####Output session info####
### Session info and metadata
setwd(OutputDirectory)
dir.create("Output - info", showWarnings = FALSE)
setwd("Output - info")

sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
session_info()
sink()

write(cellular.cols, "cellular.cols.txt")
write(cluster.cols, "cluster.cols.txt")



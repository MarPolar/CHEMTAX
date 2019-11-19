#CHEMTAX for HPLC data
setwd("~/HPLC_Chemtax")
install.packages("BCE")
install.packages("limSolve")
# load packages required
library(limSolve) # for CHEMTAX analyses

##Import data from MSM80 CUSCO Cruise  December 2018- January 2019 in Peruvian waters.##

#1.Pigment composition of the various algal groups (based on DiTullio et al 2005 Table 1 transects 8-13 specific for Peruvian waters)

Ratio <- read.table(file = "Ratio.txt", # needs to be a matrix
                    header = TRUE) # be careful here that the group names are not in the matrix themselves...
str(Ratio)

ratio <- as.matrix(Ratio) # ensure this remains as a matrix! This is important for further lines of script below.

subset(ratio[,1:16])#to check how the matrix looks like.

#2.Pigment composition of the samples measured with HPLC. Ratio pigment:Chla calculated
#First try just with the experiments 1 and 2 data

Field <- read.table(file = "Field.txt", # needs to be converted or extracted as a vector
                    header = TRUE)
str(Field)

Field.n <- Field[,-1]# removes the first column of names
Field <- Field.n[1,] # select the first sample
Field.v <- as.vector(unlist(Field)) # extracts data for 1st sample as vector called Field.v

str(Field.v)

# now create list of Ratio and Field data for sample 1. This is called CHEMTAX.
CHEMTAX <- list(Ratio = ratio, Field = Field.v)

# 1. Graphical representation of chemtax example input data
mp     <- apply(CHEMTAX$Ratio, MARGIN = 2, max) 
pstars <- rbind(t(t(CHEMTAX$Ratio)/mp),
                sample = CHEMTAX$Field/max(CHEMTAX$Field))

# this plot seems to distract a little but kept it here for good measure (graphical representation)
stars(pstars, len = 0.9, key.loc = c(7.2, 1.7),scale=FALSE,ncol=4,
      main = "CHEMTAX pigment composition", draw.segments = TRUE,
      flip.labels=FALSE)

# 2. Estimating the algal composition of the field sample
Nx     <- nrow(CHEMTAX$Ratio) #  int 10

# equations that have to be met exactly Ex=f: 
# sum of all fraction must be equal to 1.
EE <- rep(1, Nx) ## must equal 1 as full fraction # num [1:10] 
FF <- 1  #  num 1

# inequalities, Gx>=h:
# all fractions must be positive numbers
GG <- diag(nrow = Nx) # num [1:10, 1:10] 
HH <- rep(0, Nx) # num [1:10]

# equations that must be reproduced as close as possible, Ax ~ b
# = the field data; the input ratio matrix and field data are rescaled
AA     <- CHEMTAX$Ratio/rowSums(CHEMTAX$Ratio) # num [1:10, 1:16]  
BB     <- CHEMTAX$Field/sum(CHEMTAX$Field) # Named num [1:12] --> should be 10! or 16? #mine is 1:16

# now for limSolve! linear regression? This is the calculation step
# 1. Solve with lsei method
X <- lsei(t(AA), BB, EE, FF, GG, HH)$X # t is for transposition of matrix AA
(Sample <- data.frame(Algae = rownames(CHEMTAX$Ratio),
                      fraction = X))

# plot results
barplot(X, names = rownames(CHEMTAX$Ratio), col = heat.colors(8),
        cex.names = 0.8, main = "CHEMTAX analysis",las=2,
        ylab = 'Contribution to Chla')

#### For multiple samples --> how to loop it using functions etc #### by Allanah Paul

# this function is designed to extract pigment ratios from each sample as vectors
extract.field.dat <- function(i) {
  Field <- Field.n[i,] # select the first sample 
  Field.v <- as.vector(unlist(Field)) # extracts data for 1st sample as vector called Field.v
}

# this line produces a list 'CHEMTAX.Field.dat' of samples and corresponding pigment ratios to Chla
CHEMTAX.Field.dat <- lapply(seq_along(Field[,1]), extract.field.dat) # apply function along sample list

# now create list of Ratio and Field data for each sample within list CHEMTAX.l. 
#This is called CHEMTAX.l where l = list.

CHEMTAX.l <- list() # create list

# create function 'list.CHEMTAX'
list.CHEMTAX <- function(i) {
  Field.temp = CHEMTAX.Field.dat[[i]]
  CHEMTAX.l[[i]] <- list(ratios = Ratio, field = Field.temp)
  # names(CHEMTAX.l[[i]]) <- Field[[,i]]
}

# apply function 'list.CHEMTAX' along sample list to create a list of all samples with two components
# 1) ratios as a matrix and 2) pigment concs as a vector
CHEMTAX.list <- lapply(seq_along(CHEMTAX.Field.dat), list.CHEMTAX) 

# renames the list items with the sample number for easy identification
names(CHEMTAX.list) <- Field[,1] 

##### now the data is ready, time for running the CHEMTAX analysis part! ####

# create new function called 'CHEMTAX.analysis'
# note: Ratio --> ratios and Field --> field to distinguish between individual analysis and looping
CHEMTAX.output <- list() # creates list for output
CHEMTAX.analysis <- function(i){
  mp     <- apply(CHEMTAX.list[[i]]$ratios, MARGIN = 2, max) 
  pstars <- rbind(t(t(CHEMTAX.list[[i]]$ratios)/mp),
                  sample = CHEMTAX.list[[i]]$field/max(CHEMTAX.list[[i]]$field))
  # this plot seems to distract a little but kept it here for good measure (graphical representation)
  # stars(pstars, len = 0.9, key.loc = c(7.2, 1.7),scale=FALSE,ncol=4,
  #  main = "CHEMTAX pigment composition", draw.segments = TRUE,
  # flip.labels=FALSE)
  # 2. Estimating the algal composition of the field sample
  Nx     <- nrow(CHEMTAX.list[[i]]$ratios) #  int 10
  
  # equations that have to be met exactly Ex=f: 
  # sum of all fraction must be equal to 1.
  EE <- rep(1, Nx) ## must equal 1 as full fraction # num [1:10] 
  FF <- 1  #  num 1
  
  # inequalities, Gx>=h:
  # all fractions must be positive numbers
  GG <- diag(nrow = Nx) # num [1:10, 1:10] 
  HH <- rep(0, Nx) # num [1:10]
  
  # equations that must be reproduced as close as possible, Ax ~ b
  # = the field data; the input ratio matrix and field data are rescaled
  AA     <- CHEMTAX.list[[i]]$ratios/rowSums(CHEMTAX.list[[i]]$ratios) # num [1:10, 1:16]  
  BB     <- CHEMTAX.list[[i]]$field/sum(CHEMTAX.list[[i]]$field) # Named num [1:12] --> should be 10! or 16?
  
  # now for limSolve! linear regression? This is the calculation step
  # 1. Solve with lsei method
  X <- lsei(t(AA), BB, EE, FF, GG, HH)$X # t is for transposition of matrix AA
  Sample <- data.frame(Algae = rownames(CHEMTAX.list[[i]]$ratios),
                       fraction = X)
}

# apply 'CHEMTAX.analysis' function along the list 'CHEMTAX.list' 
# to add to a new list created above called 'CHEMTAX.output'
CHEMTAX.output <- lapply(seq_along(CHEMTAX.list), CHEMTAX.analysis) 
# renames the list items with the sample number for easy identification
names(CHEMTAX.output) <- Field[,1] 

# this output can then be used for plotting :) 
# a suggestion would be to create another function that then plots a stacked bar chart
# by station (x-axis) and relative contribution to Chla (y-axis)

# this function selects just the fraction data for extraction into a new df in following steps
output.extract <- function(i) {
  Fraction <- unlist(CHEMTAX.output[[i]]$fraction) 
}

# using the function, creates new list of output
CHEMTAX.output.list <- lapply(seq_along(CHEMTAX.output), output.extract)

# create data frame from output
CHEMTAX.output.cbind.df <- as.data.frame(do.call(cbind, CHEMTAX.output.list)) # can also make matrix

# data wrangling :) ....
# define new col and row names in df
phyto.names <- CHEMTAX.output[[47]]$Algae
rownames(CHEMTAX.output.cbind.df) <- phyto.names
colnames(CHEMTAX.output.cbind.df) <- Field[,1] 
CHEMTAX.output.cbind.df$Phyto.group<-CHEMTAX.output[[47]]$Algae #adds phyto group names as a column
str(CHEMTAX.output.cbind.df)

# convert data to usable format for plotting a stacked bar graph (long format required)
library(reshape2)
melted <- melt(CHEMTAX.output.cbind.df)
str(melted)

## now plotting with ggplot 
# variable = sample no.
# value = fraction of chla
stacked.bar <- ggplot(melted, aes(x = variable, y = value, fill = Phyto.group)) + 
  geom_bar(stat = 'identity') + 
  ggtitle('CHEMTAX analyses') +
  scale_x_discrete(name = 'Sample no.') +
  scale_y_continuous(name = 'Contribution to Chl a', #limits = c(0,1),
                     expand = c(0,0)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3))

stacked.bar 

####### END ############
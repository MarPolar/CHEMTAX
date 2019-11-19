#CHEMTAX for HPLC data
setwd("~/HPLC_Chemtax")
install.packages("BCE")
install.packages("limSolve")

#Import data from MSM80 CUSCO Cruise  December 2018- January 2019 in Peruvian waters.
install.packages("readxl")
library(readxl)
#1.Pigment composition of the various algal groups (based on DiTullio et al 2005 Table 1 transects 8-13 specific for Peruvian waters)

ratio<-read_excel("Ratio.xlsx")
str(ratio)
#2.Pigment composition of the samples measured with HPLC. Ratio pigment:Chla calculated
#First try just with the experiments 1 and 2 data

field<-read_excel("Field.xlsx")
str(field)


# 1. Graphical representation of the chemtax example input data
palette(rainbow(12, s = 0.6, v = 0.75))

mp     <- apply(Chemtax$Ratio, MARGIN = 2, max)
pstars <- rbind(t(t(Chemtax$Ratio)/mp) ,
                sample = Chemtax$Field/max(Chemtax$Field))
stars(pstars, len = 0.9, key.loc = c(7.2, 1.7),scale=FALSE,ncol=4,
      main = "CHEMTAX pigment composition", draw.segments = TRUE,
      flip.labels=FALSE)

# 2. Estimating the algal composition of the field sample
#nrow returns the number of rows present in x.
Nx     <-nrow(ratio)

# equations that have to be met exactly Ex=f: 
# sum of all fraction must be equal to 1.
EE <- rep(1, Nx)
FF <- 1

# inequalities, Gx>=h:
# all fractions must be positive numbers
GG <- diag(nrow = Nx)
HH <- rep(0, Nx)

# equations that must be reproduced as close as possible, Ax ~ b
# = the field data; the input ratio matrix and field data are rescaled
AA     <- ratio/rowSums(ratio)
BB     <- field/sum(field)

# 1. Solve with lsei method
X <- lsei(t(AA), BB, EE, FF, GG, HH)$X
(Sample <- data.frame(Algae = rownames(Chemtax$Ratio),
                      fraction = X))

# plot results
barplot(X, names = rownames(Chemtax$Ratio), col = heat.colors(8),
        cex.names = 0.8, main = "Chemtax example solved with lsei")

# 2. Bayesian sampling; 
# The standard deviation on the field data is assumed to be 0.01
# jump length not too large or NO solutions are found!
xs <- xsample(t(AA), BB, EE, FF, GG, HH, sdB = 0.01, jmp = 0.025)$X
pairs(xs, main= "Chemtax, Bayesian sample")
# }

#Author(s)
#Karline Soetaert <karline.soetaert@nioz.nl>.
#References
#Mackey MD, Mackey DJ, Higgins HW, Wright SW, 1996. CHEMTAX - A program for estimating class abundances from chemical markers: Application to HPLC measurements of phytoplankton. Marine Ecology-Progress Series 144 (1-3): 265-283.
#Van den Meersche, K., Soetaert, K., Middelburg, J., 2008. A Bayesian compositional estimator for microbial taxonomy based on biomarkers. Limnology and Oceanography Methods, 6, 190-199.
#R-package BCE See Also
#lsei, the function to solve for the algal composition of the field sample.

#Another example from: https://github.com/low-decarie/Useful-R-functions/blob/master/Chemtax/chemtax.R
chemtax<-function(sample.conc, pigment.matrix, baysian=F){
  
  
  sample.a<-data.matrix(sample.conc[,-1])
  dimnames(sample.a)[[1]]<-sample.conc[,1]
  
  pigment.a<-data.matrix(pigment.matrix[,-1])
  dimnames(pigment.a)[[1]]<-pigment.matrix[,1]
  
  #match names
  sample.a<-sample.a[ ,dimnames(sample.a)[[2]] %in% dimnames(pigment.a)[[2]]]
  pigment.a<-pigment.a[ ,dimnames(pigment.a)[[2]] %in% dimnames(sample.a)[[2]]]
  sample.a<-sample.a[,match(dimnames(pigment.a)[[2]],dimnames(sample.a)[[2]])]
  
  #Check match
  print(match(dimnames(pigment.a)[[2]],dimnames(sample.a)[[2]]))
  
  #Remove unidentifiable groups
  pigment.a<-pigment.a[rowSums(pigment.a)!=0,]
  
  
  pdf("./Plots/baysian baysian.pdf")
  
  X<-apply(X=sample.a, MARGIN=1,FUN=function(sample.a){
    #From Chemtax in limsolve
    # 2. Estimating the algal composition of the field sample.a
    Nx     <-nrow(pigment.a)
    
    # equations that have to be met exactly Ex=f: 
    # sum of all fraction must be equal to 1.
    EE <- rep(1,Nx)
    FF <- 1
    
    # inequalities, Gx>=h:
    # all fractions must be positive numbers
    GG <- diag(nrow=Nx)
    HH <- rep(0,Nx)
    
    # equations that must be reproduced as close as possible, Ax ~ b
    # = the field data; the input ratio matrix and field data are rescaled
    AA     <- pigment.a/rowSums(pigment.a)
    BB     <- sample.a/sum(sample.a)
    
    # 1. Solve with lsei method
    X <-lsei(t(AA),BB,EE,FF,GG,HH)$X
    
    
    #   # 2. Bayesian baysian; 
    #   # The standard deviation on the field data is assumed to be 0.01
    #   # jump length not too large or NO solutions aer found!
    xs <- xsample(t(AA),BB,EE,FF,GG,HH, sdB=0.01, jmp=0.025)$X
    pairs(xs, main= "Chemtax, Bayesian sample")
    
    
    return(X)
  })
  
  
  dev.off()
  
  X<-data.frame(t(X))
  
  X$Sample<-rownames(X)
  
  return(X)
}
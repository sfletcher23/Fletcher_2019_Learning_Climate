

#######################################
# SETUP
#######################################

# Enable running on SLURM scheduler on cluster

# Get environment var from Slurm
JOBID = Sys.getenv("SLURM_JOB_ID")

# Check if on Slurm. If not, set working directory
if (JOBID != "") {
  cat("Slurm Job ID Found")
} else {
  setwd("/Users/sarahfletcher/Dropbox (MIT)/Mombasa_Climate/BMA_code")
  cat("Slurm Job ID Not Found")
}

source("./REA.Gibbs.r")
install.packages("foreach", repos = "http://cran.us.r-project.org")
library(foreach)

# Need to use this package for the dopar commands to work
install.packages("doParallel", repos = "http://cran.us.r-project.org")
library(doParallel)

# Set up the parallelization. If running on slurm, use the number of cores available in the job. If running on desktop, use 2. 
if (JOBID != "") {
  registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
} else {
  registerDoParallel(cores=2)
}

# Get date for file save name
currentDate = Sys.Date()


#######################################
# RUN BAYESIAN ANALYSIS
#######################################

args(REA.Gibbs)


#  Loop over each possible climate change observation in each time period
year = c(1990,2010,2030,2050,2070,2090)
foreach (scen_ii = 1:48) %dopar%{
  foreach (ii = 1:5) %dopar%{
    
    # Set time periods
    yearX = year[1] # Historical time period
    yearY = year[ii+1] # Future time period
    
    # Load climate model projections for historical time period
    str1 = sprintf("Input/X_%d.csv",yearX)
    tmp = read.csv(str1,header = FALSE)
    X = as.matrix(tmp)
    
    # Load climate model projections for future time period
    str2 = sprintf("Input/X_%d.csv",yearY)
    tmp = read.csv(str2,header = FALSE)
    Y = as.matrix(tmp)
    
    # Load samples of lambda prior
    str3 = sprintf("Input/lambda0_%d.csv",yearX)
    tmp = read.csv(str3,header = FALSE)
    lambda0 = as.matrix(tmp)

    # Load virtual climate change observation for this time period
    tmp = read.csv("Input/X0PU.csv",header = FALSE)
    X0P = tmp[c(ii),c(scen_ii)]
   
    # Run Bayesian model 
    REA.Gibbs(X[c(2),],X0P[c(1)],lambda0[c(2)],Y[c(2),],N=1000)->rg0
    
    # Save output
    mu0str = paste(sprintf("Output/muUP_%d_scen%d",yearY,scen_ii),"job", JOBID, currentDate, ".csv", sep="_")
    nu0str = paste(sprintf("Output/nuUP_%d_scen%d",yearY,scen_ii),"job", JOBID, currentDate, ".csv", sep="_")
    lambdastr = paste(sprintf("Output/lambdaUP_%d_scen%d",yearY,scen_ii),"job", JOBID, currentDate, ".csv", sep="_")
    write.table(rg0$lambda[,c(1:21),c(1:1000)], file = lambdastr,row.names=FALSE, col.names=FALSE,sep=",")
    write.table(rg0$mu, file = mu0str,row.names=FALSE, col.names=FALSE)
    write.table(rg0$nu, file = nu0str,row.names=FALSE, col.names=FALSE)
  }
}


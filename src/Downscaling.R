path_matlabRuntime="" #Where Matlab runtime is installed
path=""#Path to your working directory
Freq=12#Frequency of your data. For example if its is monthly data,frequency is 12. If it is daily data, frequency is 365.

Downscaling=function(path_matlabRuntime,path,Freq){
  # Check if the package is installed; if not, install it
  # Define the package name you want to load
  package_name <- "data.table"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  # Define the package name you want to load
  package_name <- "akima"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  # Define the package name you want to load
  package_name <- "parallel"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  # Define the package name you want to load
  package_name <- "doParallel"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  # Define the package name you want to load
  package_name <- "foreach"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  # Define the package name you want to load
  package_name <- "pryr"
  
  # Check if the package is installed; if not, install it
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    install.packages(package_name, repos = "http://cran.us.r-project.org")
    library(package_name, character.only = TRUE)
  }
  
  ##################################################################################################################################################
  # Load data from CSV files
  Model=as.matrix(data.table::fread(paste0(path,"/Model_data.csv")))
  Obs=as.matrix(data.table::fread(paste0(path,"/Obs_data.csv")))
  
  T_o=ncol(Obs[,-c(1:2)])#Length of the observation period
  T_all=ncol(Model[,-c(1:2)])#Length of the future downscaling period
  N=nrow(Obs)#Total number of fine resolution locations
  
  cat("Data pre-processing starts!\n")
  source(paste0(path,"/src/Stand_downscaling.R"))
  Stand_SD=Stand_downscaling(Obs,Model,Freq)
  Stand_SD=cbind(Obs[,1:2],Stand_SD)
  fwrite(Stand_SD,paste0(path,"Stand_Downscaled_training.csv"))
  cat("Data pre-processing ends!\n")
  
  cat("Mean estimation starts!\n")
  # Initialize matrices to store results
  Mu=matrix(NA, nrow = N, ncol = T_all)
  for(t in 1:T_all){
      if(t==1){
        t_seq=c(t+Freq,t+1)
        Mu[,t]=rowMeans(Stand_SD[,t_seq+2]) 
      }else if(t>1 & t<Freq-1){
        t_seq=c(t+Freq,t-1)
        Mu[,t]=rowMeans(Stand_SD[,t_seq+2]) 
      }else{
        t_seq=c(t-12,t-1)
        Mu[,t]=rowMeans(Stand_SD[,t_seq+2]) 
    }}
  cat("Mean estimation ends!\n")
    
  # Write the results to CSV files
  fwrite(Mu,paste0(path,"Mu1.csv"))
  fwrite(Mu,paste0(path,"Mu2.csv"))

  # Create a command to run an external script
  items=c(paste0(path,"src/run_BGL.sh"),path_matlabRuntime,path)
    
  # Wrap each item with double quotes and collapse them into a single string
  command=paste0("\"", paste(items, collapse = "\" \""), "\"")
    
  # Execute the external script using the 'system' function
  system(print(command))
  
}

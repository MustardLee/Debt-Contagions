library(pracma)
library(readr)
library(dplyr)

# 导入数据集
crossholding <- read_csv("crossholding.csv")
gdp <- read_csv("gdp21.csv")
gdp19 <- read_csv("gdp19.csv")

# 共有25个国家
n = 25
# 存入国家名称
names = t(t(gdp$Country))

# normalized country GDP(nominal)
gdp_21 = matrix(gdp$gdp_normalized)
gdp_19 = matrix(gdp19$gdp_normalized)

# save crossholding data as dataframe
C_raw <- data.frame(crossholding)
# convert to matrix
C_raw_m <- data.matrix(C_raw[2:26])
# sum every country column
sum_C <- t(data.matrix(colSums(C_raw[,2:26])))

# Here we estimate the ratio of total debt held outside the issuing country by 1/3.
c <- .33

# get C bar and save as diag matrix
C_row <- as.vector(matrix(rep(1,25), ncol = 25)/sum_C)
C_diag <- diag(C_row)

# normalizes columns to sum to 1
# imposes that only c of a country is held outside
C <- c * C_raw_m %*% C_diag

# the remainder is held by private shareholders, i.e. the citizens
Chat <- (1-c)*diag(n) 

# get dependency matrix A
A <- Chat%*%solve(diag(n)-C)
A_df <- data.frame(A)

# name A_df column and row to country, then save to csv
names(A_df) <- t(gdp$Country)
row.names(A_df) <- t(gdp$Country)
write_csv(A_df, "A.csv")

# p as predicted gdp
p <- gdp_21
pcurrent <- as.vector(p)
beta <- 0.5

# function for wave trend by theta
contagion <- function(theta){
  # set wave number
  wave <<- 1
  first_failed_num <- 0
  
  # define theta and threshold
  v_threshold <- theta*(A %*% gdp_19)
  
  # set failure indicator and failed country saver
  all_failure_indicator <- matrix(rep(0),25,1)
  new_failed_countries <- c(1)
  all_failure_list <- list()
  
  # loop for waves
  while(!isempty(new_failed_countries)){
    
    a <- as.matrix((as.numeric(A %*% pcurrent < v_threshold)))
    b <- as.matrix(as.numeric(all_failure_indicator == 0))
    
    new_failure_indicator <-  a * b
    
    # find pairwise maximum
    all_failure_indicator <- pmax(all_failure_indicator, new_failure_indicator)
    
    # cut the failed country's gpd by half of threshold
    pcurrent <- pcurrent - new_failure_indicator * (v_threshold * beta)
    
    # find new_failed_countries 
    new_failed_countries <- which(new_failure_indicator == 1)
    # get fialed country names
    new_failed_names = matrix(names[new_failed_countries], nrow = 1)
    
    if(wave == 1) {
      first_failed_num <<- length(new_failed_countries)
    }
    
    if (!isempty(new_failed_countries)){
      # cat('Wave ', wave, ' failures are ')
      # cat(new_failed_names)
      all_failure_list[[length(all_failure_list)+1]] <- list(new_failed_names)
      
    }
    else{
      break
    }
    # cat("\n","end wave", wave,"\n")
    wave <<- wave + 1
  }
  all_failed_num <<- length(which(all_failure_indicator==1))
  print(all_failure_list)
}

theta = c()
wave_length = c()
first_fail = c()
all_failed = c()

i <- 1
# for loop of theta from .9 to 1.2
for (t in seq(0.9, 1.2, by = 0.01)){
  cat("start loop", i, "\n")
  contagion(t)
  theta <- append(theta, t)
  wave_length <- append(wave_length, wave-1)
  first_fail <- append(first_fail, first_failed_num)
  all_failed <- append(all_failed, all_failed_num)
  i <- i+1
}

wave.df<- data.frame(
  Theta = theta,
  Total_Wave = wave_length,
  First_Fail_Num = first_fail,
  All_Fail_Num = all_failed
)

write_csv(wave.df, "wave.csv")
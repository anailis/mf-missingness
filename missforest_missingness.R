# Methods: Madley-Dowd, P. et al. (2019) 'The proportion of missing data should not be used to guide decisions on multiple imputation', 
# Journal of Clinical Epidemiology. Elsevier USA, 110, pp. 63-73. 
# doi: 10.1016/j.jclinepi.2019.02.016.

pacman::p_load(tidyverse, MASS, missForest, DBI, here)

create_data <- function(seed) {
  # all variables other than Y are not correlated with each other 
  # Y and X have a correlation of 0.6
  # Y and Z1:Z2 have a correlation of 0.4
  # Y and Z3:Z7 have a correlation of 0.2
  # Y and Z8:Z11 have a correlation of 0.1
  variable_names <- c('Y', 'X', 'Z1', 'Z2', 'Z3', 'Z4', 'Z5', 'Z6', 'Z7', 'Z8', 'Z9', 'Z10', 'Z11')
  
  covar <- matrix(data = c(1, 0.6, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1,
                         0.6, 1, rep.int(0, 11),
                         0.4, 0, 1, rep.int(0, 10),
                         0.4, rep.int(0, 2), 1, rep.int(0, 9),
                         0.2, rep.int(0, 3), 1, rep.int(0, 8),
                         0.2, rep.int(0, 4), 1, rep.int(0, 7),
                         0.2, rep.int(0, 5), 1, rep.int(0, 6),
                         0.2, rep.int(0, 6), 1, rep.int(0, 5),
                         0.2, rep.int(0, 7), 1, rep.int(0, 4),
                         0.1, rep.int(0, 8), 1, rep.int(0, 3),
                         0.1, rep.int(0, 9), 1, rep.int(0, 2),
                         0.1, rep.int(0, 10), 1, rep.int(0, 1),
                         0.1, rep.int(0, 11), 1), 
                ncol = 13,
                nrow = 13,
                byrow = TRUE,
                dimnames = list(variable_names, variable_names))

  set.seed(seed)
  sim_data <- mvrnorm(n = 1000,
                    mu = rep.int(1, 13),
                    Sigma = covar)
  return(sim_data)
}

# generate NAs in Y missing completely at random (mcar)
mcar <- function(dataMatrix, percentMissing) {
  numRowsToRemove <- (nrow(dataMatrix)/100) * percentMissing
  dataMatrix[1:numRowsToRemove,"Y"] <- NA
  return(dataMatrix)
}

#
mar <- function(dataMatrix, percentMissing) {
  # logistic regression
}

# impute data with the missforest algorithm 
impute <- function(dataMatrix, variables) {
  dataMatrix <- dataMatrix[,variables]
  results <- missForest(dataMatrix)
  return(results)
}

analyse <- function(dataMatrix, varset) {
  dataDf <- as.data.frame(dataMatrix)[,c('Y', varset)]
  model <- lm(Y ~ ., data = dataDf)
  return(model)
}

# store results in a database 
store_results <- function(results, runID, percentMissing, impMethod, impVars, db) {
  coef <- summary(results)$coefficients['X','Estimate']
  # because X has a correlation of 0.6 with Y
  # only considering the exposure coefficient, not the constant coefficient here 
  bias <- 0.6 - coef
  # coefficient se
  se <- summary(results)$coefficients['X','Std. Error']
  impVars <- paste(impVars,collapse=" ")
  
  dbSendQuery(conn = db,
            "INSERT INTO results
            VALUES (?, ?, ?, ?, ?, ?, ?)", 
            c(runID, percentMissing, impMethod, impVars, coef, bias, se))
}

# SET UP DATABASE

resultsDb <- dbConnect(RSQLite::SQLite(), 
                       here("Learning", "Blog", "missforest", "missForest.sqlite"))

dbSendQuery(conn = resultsDb,
            "CREATE TABLE results
            (runID INTEGER,
             percentMissing INTEGER,
             imputationMethod TEXT,
             imputationVars TEXT,
             regressCoef REAL,
             bias REAL,
             se REAL)")

# MISSING COMPLETELY AT RANDOM

# 1%, 5%, 10%, 20%, 40%, 60%, 80%, and 90% missing

# 1000 datasets per scenario
# generate 1000 random seeds
set.seed(5)
seeds <- sample(1:1000000000, 1000, replace=FALSE)
n = 0

# variables to be included in imputation: Y, X
# variables to be included in imputation: Y, X, Z1
# variables to be included in imputation: Y, X, Z3
# variables to be included in imputation: Y, X, Z1 - Z4
# variables to be included in imputation: Y, X, Z1 - Z11
for (percent in c(1, 5, 10, 20, 40, 60, 80, 90)) {
  for (seed in seeds) {
    sim_data <- create_data(seed)
    sim_data <- mcar(sim_data, percent)
  
    for (varset in list(c("X", "Z1"), 
                      c("X", "Z3"), 
                      c("X", "Z1", "Z2", "Z3", "Z4"),
                      c("X", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z8", "Z9", "Z10", "Z11"))) {
      results_cc <- analyse(sim_data, varset)
      n = n + 1
      store_results(results_cc, n, percent, 'complete cases', varset, resultsDb)
      
      imputed <- impute(sim_data, varset)
      results_imp <- analyse(imputed$ximp)
      n = n + 1
      store_results(results_imp, n, percent, 'missForest', varset, resultsDb)
    }
  }
}

# MISSING AT RANDOM

# 1%, 5%, 10%, 20%, 40%, 60%, 80%, and 90% missing
# 1000 datasets per scenario 

# variables to be included in imputation: Y, X
# variables to be included in imputation: Y, X, Z1
# variables to be included in imputation: Y, X, Z3
# variables to be included in imputation: Y, X, Z1 - Z4
# variables to be included in imputation: Y, X, Z1 - Z11

dbDisconnect(resultsDb)

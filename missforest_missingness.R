# Methods: Madley-Dowd, P. et al. (2019) 'The proportion of missing data should not be used to guide decisions on multiple imputation', 
# Journal of Clinical Epidemiology. Elsevier USA, 110, pp. 63-73. 
# doi: 10.1016/j.jclinepi.2019.02.016.

pacman::p_load(tidyverse, MASS, missForest, DBI, here)

create_data <- function() {
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

  sim_data <- mvrnorm(n = 1000,
                    mu = rep.int(1, 13),
                    Sigma = covar)
}

mcar <- function(dataMatrix, percentMissing) {
  numRowsToRemove <- (nrow(dataMatrix)/100) * percentMissing
  dataMatrix[1:numRowsToRemove,"Y"] <- NA
  return(dataMatrix)
}

mar <- function(dataMatrix, percentMissing) {
  # logistic regression
}

impute <- function(dataMatrix, variables) {
  dataMatrix <- dataMatrix[,variables]
  results <- missForest(dataMatrix)
  return(results)
}

analyse <- function(imputedMatrix) {
  return(lm(imputedMatrix[,"Y"] ~ imputedMatrix[,"X"]))
}

store_results <- function(results, runID, percentMissing, impMethod, impVars, db) {
  bias <- 0.6 - results$coefficients
  # empirical se
  empSE <- 0
  # fraction of missing information - see paper for how to calculate
  fmi <- 0
  
  dbSendQuery(conn = db,
            "INSERT INTO results
            VALUES (?, ?, ?, ?, ?, ?, ?)", 
            c(runID, percentMissing, impMethod, impVars, bias, empSE, fmi))
}

# SET UP DATABASE

resultsDb <- dbConnect(RSQLite::SQLite(), 
                       here("Learning", "Blog", "missforest", "missForest.sqlite"))

#dbSendQuery(conn = resultsDb,
            # "CREATE TABLE results
            # (runID INTEGER,
            #   percentMissing INTEGER,
              # imputationMethod TEXT,
              # imputationVars TEXT,
              # bias REAL,
              # se REAL,
              # fmi REAL)")

dbSendQuery(conn = resultsDb,
            "INSERT INTO results
            VALUES (1, 2, 'ignore', 'me', 1.12, 1.34, 5.6666)")

dbWriteTable()

# MISSING COMPLETELY AT RANDOM

#1%, 5%, 10%, 20%, 40%, 60%, 80%, and 90% missing
# 1000 samples per scenario

for (percent in c(1, 5, 10, 20, 40, 60, 80, 90)) {
  sim_data <- create_data()
  sim_data <- mcar(sim_data, percent)
  
  for (varset in list(c("Y", "X", "Z1"), 
                      c("Y", "X", "Z3"), 
                      c("Y", "X", "Z1", "Z2", "Z3", "Z4"),
                      c("Y", "X", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z8", "Z9", "Z10", "Z11"))) {
    imputed <- impute(sim_data, varset)
    results <- analyse(imputed$ximp)
    #store_results(results, resultsDb)
  }
}

# variables to be included in imputation: Y, X
# variables to be included in imputation: Y, X, Z1
# variables to be included in imputation: Y, X, Z3
# variables to be included in imputation: Y, X, Z1 - Z4
# variables to be included in imputation: Y, X, Z1 - Z11

for (run in 1:1000) {
  for (percent in c(1, 5, 10, 20, 40, 60, 80, 90)) { 
    sim_data <- create_data()
    sim_data <- mcar(sim_data, percent)
    results <- analyse(sim_data)
    store_results(results, resultsDb)
  }
}

# complete cases

# MISSING AT RANDOM

# 1%, 5%, 10%, 20%, 40%, 60%, 80%, and 90% missing
# 1000 samples per scenario 

# variables to be included in imputation: Y, X
# variables to be included in imputation: Y, X, Z1
# variables to be included in imputation: Y, X, Z3
# variables to be included in imputation: Y, X, Z1 - Z4
# variables to be included in imputation: Y, X, Z1 - Z11
# complete cases 

dbDisconnect(resultsDb)

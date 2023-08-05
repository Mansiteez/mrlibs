#' To Find Differentially Expressed Genes
#'
#' This function is used to find Differentially Expressed Genes from 2 or more sub groups of data.
#'
#' @param data_file Path to the data file containing gene expression levels with columns containing gene names.
#' @param dataID_file Path to the dataID file containg columns as SampleId and SubType.
#'
#' @return A vector containing the results of the function.
#'
#' @examples
#' finDef("data.csv", "dataID.csv")
#'
finDef <- function(data_file, dataID_file, fold) {
  library(readr)
  library(siggenes)
  library(pamr)

  # Read the data file
  data <- read.csv(data_file)

  # Read the dataID file
  dataID <- read_csv(dataID_file)

  # Extract the class labels
  cl <- c(dataID$Subtype)
  cl

  # Transpose the data matrix
  TData <- t(data)
  TData

  # Perform SAM analysis
  sam.out <- sam(TData, cl, rand = 123)
  View(sam.out)

  # Obtain the summary table of SAM results
  table <- summary(sam.out)
  table

  # Initialize a vector to store the results
  results <- c()

  # Convert the summary table to a data frame
  table_df <- as.data.frame(table@mat.fdr)

  # Iterate over each delta value in the summary table
  for (i in 1:nrow(table_df)) {
    delta_value <- table_df$Delta[i]

    # Select list of genes based on the delta value
    DEG <- list.siggenes(sam.out, delta_value)
    DEGdf <- subset(data, select = c(DEG))
    TDEGdf <- t(DEGdf)
    TrainData <- list(x = TDEGdf, y = cl)

    # Train the PAM model
    PAMtable <- pamr.train(TrainData, fold = fold)

    # Perform cross-validation
    CVtable <- pamr.cv(PAMtable, TrainData)

    # Plot the cross-validation results
    result <- pamr.plotcv(CVtable)
    results <- c(results, result)

    # Calculate the threshold with minimum error
    threshold <- CVtable$threshold[which(CVtable$error == min(CVtable$error))][1]

    # Calculate confusion matrix using the optimal threshold
    pamr.confusion(CVtable, threshold = threshold)

    # Generate a list of differentially expressed genes
    pamcen <- pamr.listgenes(PAMtable, TrainData, threshold = threshold)

    # Update column names of the gene list
    colnames(pamcen) <- c("genes", colnames(pamcen)[-1])

    # Print the dimensions of the gene list
    print(dim(pamcen))

    # Export pamcen to a CSV file
    csv_file <- paste0("pamcen_", delta_value, ".csv")
    write.csv(pamcen, file = csv_file, row.names = FALSE)
  }

  # Return results
  return(results)
}

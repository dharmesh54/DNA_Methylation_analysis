# Ensure the necessary libraries are loaded
library(dplyr)

setwd("C:/Users/dharm/Downloads/Pupil_bio/")

# Define the files and variable names
files <- c("cfdna_m1_r1.txt", "cfdna_m2_r1.txt", "cfdna_m3_r1.txt",
           "cfdna_m1_r2.txt", "cfdna_m2_r2.txt", "cfdna_m3_r2.txt",
           "islet_m1.txt", "islet_m2.txt", "islet_m3.txt")

variable_names <- c("cfdna_m1_r1", "cfdna_m2_r1", "cfdna_m3_r1",
                    "cfdna_m1_r2", "cfdna_m2_r2", "cfdna_m3_r2",
                    "islet_m1", "islet_m2", "islet_m3")

variables <- list(cfdna_m1_r1, cfdna_m2_r1, cfdna_m3_r1,
                  cfdna_m1_r2, cfdna_m2_r2, cfdna_m3_r2,
                  islet_m1, islet_m2, islet_m3)

# Check if lengths of 'files' and 'variable_names' match
if(length(files) != length(variable_names)) {
  stop("Mismatch between the number of files and variable names!")
}

# Function to summarize the data
summarize_data <- function(data) {
  total_count_column <- names(data)[3]
  methylation_column <- names(data)[4]
  print(total_count_column)
  print(methylation_column)
  
  summarized_data <- data %>%
    group_by(CpG_Coordinate) %>%
    summarize(
      Sample_ID = first(Sample_ID),
      Total_Count = sum(get(total_count_column), na.rm = TRUE),
      Methylation = sum(get(methylation_column), na.rm = TRUE)
    )
  
  return(summarized_data)
}

# Loop through each file, read it, process it, and store the result in the corresponding variable
for(i in 1:length(files)) {
  # Read the file
  data <- read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Summarize the data
  summarized_data <- summarize_data(data)
  
  # Store the summarized data in the corresponding variable
  assign(variable_names[i], summarized_data)
  
  # Print the head of the summarized data for verification
  #cat("Summary for", variable_names[i], ":\n")
  #print(head(summarized_data))
  #cat("\n")
}



variable_comparison_results <- list()






# Key-value mapping for renaming columns
rename_mapping <- c(
  "Sample_ID" = "chr",
  "CpG_Coordinate" = "pos",
  "Total_Count" = "N",
  "Methylation" = "X"
)

# Loop through each variable and rename columns
for (var_name in variable_names) {
  # Retrieve the data frame dynamically
  data <- get(var_name)
  
  # Rename columns using the mapping
  colnames(data) <- rename_mapping[colnames(data)]
  
  # Assign the updated data frame back to the variable
  assign(var_name, data)
}

# Print a sample to verify changes
print(cfdna_m1_r1)


# Start with the first dataset
common_positions <- get(variable_names[1]) %>% select(pos)

# Iteratively find common `pos` across all datasets
for (var_name in variable_names[-1]) {
  data <- get(var_name) %>% select(pos)
  common_positions <- inner_join(common_positions, data, by = "pos")
}

# Filter each dataset to keep only rows with common `pos`
for (var_name in variable_names) {
  data <- get(var_name)
  filtered_data <- data %>% filter(pos %in% common_positions$pos)
  assign(var_name, filtered_data)
}



#Tissue specific PMP's using "DSS"
#for m1 position
library("DSS")

BSobj_m1 <- makeBSseqData( list( cfdna_m1_r2, cfdna_m1_r1, islet_m1),
                        c("C1","C2", "N1") )

BSobj_m2 <- makeBSseqData( list( cfdna_m2_r2, cfdna_m2_r1, islet_m2),
                        c("C1","C2", "N1") )

BSobj_m3 <- makeBSseqData( list( cfdna_m3_r2, cfdna_m3_r1, islet_m3),
                           c("C1","C2", "N1") )

dmlTest_m1 <- DMLtest(BSobj_m1, group1=c("C1", "C2"), group2=c("N1"), equal.disp=TRUE, smoothing=TRUE)
dmlTest_m2 <- DMLtest(BSobj_m2, group1=c("C1", "C2"), group2=c("N1"), equal.disp=TRUE, smoothing=TRUE)
dmlTest_m3 <- DMLtest(BSobj_m3, group1=c("C1", "C2"), group2=c("N1"), equal.disp=TRUE, smoothing=TRUE)

write.csv(dmlTest_m1, "test_results_m1.csv")
write.csv(dmlTest_m2, "test_results_m2.csv")
write.csv(dmlTest_m3, "test_results_m3.csv")
head(dmlTest)
 




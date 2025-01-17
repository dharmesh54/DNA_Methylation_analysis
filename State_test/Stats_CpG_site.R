# Load dplyr package
library(dplyr)
# Load dplyr package

library(dplyr)
 
library(DSS)

setwd("C:/Users/dharm/Downloads/Pupil_bio/")

#Reading Input CSv file 
df = read.csv("PupilBioTest_PMP_revA.csv")

#getting information about the df
head(df)

str(df)

#counting frequency of column values  
table(cfDNA_data$strand)
table(Islet_data$Replicate)
table(df$Sample_ID)
table(df$Tissue)



# Create a new column for unmethylated values (column 'X.000')
df$unmethylated <- df$X.000

# Create a new column for methylated values (sum of columns 'X.001' to 'X.111')
df$methylated <- rowSums(df[, c("X.001", "X.010", "X.011", "X.100", "X.101", "X.110", "X.111")])

# Create a new column for total values (sum of unmethylated and methylated)
df$total_value <- df$unmethylated + df$methylated

# View the first few rows of the updated dataframe
head(df)


# Separate reverse strand
reverse_strand <- subset(df, strand == "r")

# Separate forward strand
forward_strand <- subset(df, strand == "f")




# Separate the dataframe into two based on Tissue
cfDNA_data <- reverse_strand %>% filter(Tissue == "cfDNA")
Islet_data <- reverse_strand %>% filter(Tissue == "Islet")

# View the separated data
head(cfDNA_data)
head(Islet_data)


# Extract the CpG_Coordinates columns
cfDNA_Coordinates <- cfDNA_data$CpG_Coordinates
Islet_Coordinates <- Islet_data$CpG_Coordinates

# Find the intersection of the two columns
common_coordinates <- intersect(cfDNA_Coordinates, Islet_Coordinates)


#all_same <- identical(cfDNA_data$CpG_Coordinates, Islet_data$CpG_Coordinates)



cfDNA_data_R1 <- cfDNA_data %>% filter(Replicate == "Rep1")
cfDNA_data_R2 <- cfDNA_data %>% filter(Replicate == "Rep2")



# Define the function to process data for cfDNA
process_data <- function(data_R1, data_R2) {
  # Step 1: Find common CpG_Coordinates
  common_coordinates <- intersect(data_R1$CpG_Coordinates, data_R2$CpG_Coordinates)
  
  # Step 2: Filter the dataframes to keep only common CpG_Coordinates
  data_R1_common <- data_R1 %>% filter(CpG_Coordinates %in% common_coordinates)
  data_R2_common <- data_R2 %>% filter(CpG_Coordinates %in% common_coordinates)
  
  # Remove duplicates by keeping only the first occurrence of each CpG_Coordinates
  data_R1_common <- data_R1_common %>%
    distinct(CpG_Coordinates, .keep_all = TRUE)
  
  data_R2_common <- data_R2_common %>%
    distinct(CpG_Coordinates, .keep_all = TRUE)
  
  # Merge the dataframes
  merged_data <- data_R1_common %>%
    full_join(data_R2_common, by = "CpG_Coordinates", suffix = c("_r1", "_r2"))
  
  return(merged_data)
}


cfDNA_data_merged <- process_data(cfDNA_data_R1, cfDNA_data_R2)


length(intersect(Islet_data$CpG_Coordinates, cfDNA_data_merged$CpG_Coordinates))

unique_islet = setdiff(Islet_data$CpG_Coordinates, cfDNA_data_merged$CpG_Coordinates)

unique_Islet_data <- Islet_data %>%
  filter(CpG_Coordinates %in% unique_islet)




# View the resulting dataframe
View(pmp_reverse)


cfDNA_data_merged <- cfDNA_data_merged %>%
  mutate(methylated_cfdna_avg = (methylated_r1 + methylated_r2) / 2)



#calculating median and cv for cfdna and islet in reverse strand  
# Function to calculate median and CV
calculate_median_cv <- function(df, column_name) {
  # Ensure the column exists in the dataframe
  if (!column_name %in% colnames(df)) {
    stop(paste("Column", column_name, "not found in the dataframe"))
  }
  
  # Calculate median and CV
  result <- df %>%
    summarise(
      median_value = median(.data[[column_name]], na.rm = TRUE),
      cv_value = sd(.data[[column_name]], na.rm = TRUE) / mean(.data[[column_name]], na.rm = TRUE) * 100
    )
  
  return(result)
}

calculate_median_cv(cfDNA_data_merged, "methylated_cfdna_avg")
calculate_median_cv(Islet_data, "methylated")

#output:Reverse Strand 
#cfdna
#median_value cv_value
#          6.5 358.9594

#islet
#median_value cv_value
#            2 225.8209

# Merge cfDNA_data_merged and Islet_data_merged on CpG_Coordinates
pmp_reverse <- cfDNA_data_merged %>%
  full_join(Islet_data, by = "CpG_Coordinates", suffix = c("_cfDNA", "_Islet"))

write.csv(pmp_reverse, "Phased_mutation_reveresed.csv")



# Function to process data for r1, r2, and islet

# Define the coordinates and suffixes
coordinates <- c("x", "y", "z")
suffixes <- c("r1", "r2", "")  # r1, r2, and no suffix for islet data
data <- read.csv("Phased_mutation_reveresed.csv")

# If there are missing values, remove rows with any missing values
if (any(is.na(data))) {
  data <- data[complete.cases(data), ]
}


# Set the directory where the files will be saved
output_dir <- "C:/Users/dharm/Downloads/Pupil_bio/"

# Function to process the data based on coordinate and suffix

process_data <- function(data, coordinate, suffix) {
  print(paste("Processing data for coordinate:", coordinate, "and suffix:", suffix))
  
  # Check column names in data
  print("Column names in data:")
  #print(names(data))
  
  # Determine the index based on the coordinate (x = 1, y = 2, z = 3)
  coord_index <- if (coordinate == "x") {
    1
  } else if (coordinate == "y") {
    2
  } else if (coordinate == "z") {
    3
  } else {
    stop("Invalid coordinate value. Must be 'x', 'y', or 'z'.")
  }
  
  # Extract the CpG coordinate for the specified coordinate (x, y, or z)
  data$CpG_Coordinate <- sapply(strsplit(data$CpG_Coordinates, ":"), function(x) x[coord_index])
  
  # Define the binary columns for methylation based on the suffix (r1, r2, or islet)
  if (suffix == "r1") {
    methylation_columns <- grep("^X\\.\\d{3}_r1", names(data), value = TRUE)
  } else if (suffix == "r2") {
    methylation_columns <- grep("^X\\.\\d{3}_r2", names(data), value = TRUE)
  } else {
    methylation_columns <- grep("^X\\.\\d{3}$", names(data), value = TRUE)  # For islet data, no suffix
  }
  print(methylation_columns)
  
  #Define patterns to exclude based on the coordinate
  if (coordinate == "x") {
    exclude_patterns <- c("X.000", "X.010", "X.001", "X.011")
  } else if (coordinate == "y") {
    exclude_patterns <- c("X.000", "X.101", "X.001", "X.100")
  } else if (coordinate == "z") {
    exclude_patterns <- c("X.000", "X.110", "X.100", "X.010")
  }
  
  # Filter columns by excluding patterns (match only the part before "_suffix")
  methylation_columns <- methylation_columns[!sapply(methylation_columns, function(col) {
    any(startsWith(col, exclude_patterns))
  })]
  
  
  # Check selected methylation columns
  print("Selected methylation columns:")
  print(methylation_columns)
  
  
  # Calculate methylation by summing the relevant columns
  data <- data %>%
    mutate(Methylation = rowSums(select(., all_of(methylation_columns))))
  
  # Extract total_value_r1 (for r1 and r2 files) or no suffix for islet data
  if (suffix == "r1" || suffix == "r2") {
    data$Total_Count <- data[[paste0("total_value_", suffix)]]
  } else {
    data$Total_Count <- data$total_value
  }
  
  # Select and rename the required columns
  # Select and rename the required columns
  result <- data %>%
    select(
      Sample_ID = ifelse(suffix == "", "Sample_ID", paste0("Sample_ID_", suffix)), 
      CpG_Coordinate, Total_Count, Methylation
    )
  
  # Check processed data
  print("Processed data (first few rows):")
  print(head(result))
  return(result)
}




# Function to save data to a text file
save_data <- function(data, file_name) {
  write.table(
    data, 
    file = file_name, 
    sep = "\t",          # Use tab-separated format
    row.names = FALSE,   # Do not include row numbers
    col.names = TRUE,    # Include column names
    quote = FALSE        # Do not include quotes around values
  )
}

# Process and save each dataset
cfdna_m1_r1 <- process_data(data, "x", "r1")
save_data(cfdna_m1_r1, "cfdna_m1_r1.txt")

cfdna_m2_r1 <- process_data(data, "y", "r1")
save_data(cfdna_m2_r1, "cfdna_m2_r1.txt")

cfdna_m3_r1 <- process_data(data, "z", "r1")
save_data(cfdna_m3_r1, "cfdna_m3_r1.txt")

cfdna_m1_r2 <- process_data(data, "x", "r2")
save_data(cfdna_m1_r2, "cfdna_m1_r2.txt")

cfdna_m2_r2 <- process_data(data, "y", "r2")
save_data(cfdna_m2_r2, "cfdna_m2_r2.txt")

cfdna_m3_r2 <- process_data(data, "z", "r2")
save_data(cfdna_m3_r2, "cfdna_m3_r2.txt")

islet_m1 <- process_data(data, "x", "")
save_data(islet_m1, "islet_m1.txt")

islet_m2 <- process_data(data, "y", "")
save_data(islet_m2, "islet_m2.txt")

islet_m3 <- process_data(data, "z", "")
save_data(islet_m3, "islet_m3.txt")


summarize_data(cfdna_m1_r1)



# Check if there are any duplicate rows
duplicates <-  cfdna_m2_r1 %>% filter(Total_Count == Methylation)
print(duplicates)





# This program calculates the 4 nucleotide frequencies across all samples of a FASTA file
# Test cases were conducted with two Primate Nucleotide Files
# The seqinr library was downloaded for this program, which enables the use of FASTA files as inputs

# FASTA PROCESSING
# Paste the FASTA file in the main script, replacing "path/to/your/input.fasta"
# Once the program has executed, the resulting samples and frequencies are stored as an output file in the same folder as the input file and R code


# BEGINNING OF SCRIPT
calculate_nucleotide_frequencies <- function(sequence) {

# Handle case variations by converting the sequence to uppercase
  sequence <- toupper(sequence)
  
# Counters initialized for each nucleotide

  A_count <- 0
  T_count <- 0
  C_count <- 0
  G_count <- 0
  
# Loop iterates through each character in the sequence
# If conditional evaluates to true, counter will increase by +=1

  for (base in strsplit(sequence, NULL)[[1]]) {
    if (base == "A") {
      A_count <- A_count + 1
    } else if (base == "T") {
      T_count <- T_count + 1
    } else if (base == "C") {
      C_count <- C_count + 1
    } else if (base == "G") {
      G_count <- G_count + 1
    }
# Irrelevant characters are bypassed
  }
  
# Frequencies are calculated

  total <- A_count + T_count + C_count + G_count
  freq_A <- A_count / total
  freq_T <- T_count / total
  freq_C <- C_count / total
  freq_G <- G_count / total
  
# Assign names to frequencies and return list
  return(list(A = freq_A, T = freq_T, C = freq_C, G = freq_G))
}


# Function to process a FASTA file and calculate nucleotide frequencies
process_fasta_file <- function(file_path) {
  
# Read the FASTA file
  fasta_data <- readLines(file_path)
  
# Initialize an empty list to store frequencies for each sequence
  sequence_frequencies <- list()
  
# Initialize counters for overall nucleotide frequencies
  overall_counts <- c(A = 0, T = 0, C = 0, G = 0)
  
# Each line in the FASTA file is processed
  current_sequence <- NULL
  for (line in fasta_data) {
    if (startsWith(line, ">")) {
      
# If a sequence is already being processed, calculate its frequencies
      if (!is.null(current_sequence)) {
        frequencies <- calculate_nucleotide_frequencies(current_sequence)
        sequence_frequencies[[current_sequence_header]] <- frequencies
        
# Update overall counts
        overall_counts <- overall_counts + as.numeric(unname(frequencies))
      }
      
# Initialize a new sequence
      current_sequence_header <- substr(line, 2, nchar(line))
      current_sequence <- ""
    } else {
      
# Concatenate lines to form the sequence
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
# Calculate frequencies for the last sequence
  frequencies <- calculate_nucleotide_frequencies(current_sequence)
  sequence_frequencies[[current_sequence_header]] <- frequencies
  overall_counts <- overall_counts + as.numeric(unname(frequencies))
  
# Calculate overall frequencies
  overall_frequencies <- overall_counts / sum(overall_counts)
  
# Return list containing sequence-wise and overall frequencies
  return(list(sequence = sequence_frequencies, overall = overall_frequencies))
}

}

# Function to deploy frequencies to an output file
write_frequencies_to_file <- function(output_file, sequence_frequencies, overall_frequencies) {
  # Extract the input file name without extension
  file_name <- tools::file_path_sans_ext(basename(output_file))
  
# Construct the output file path
  output_path <- file.path(dirname(output_file), paste0(file_name, "_output.txt"))
  
# Open the output file for writing
  con <- file(output_path, "w")
  
# Write header
  cat("Sequence_ID\tA\tT\tC\tG\n", file = con)
  
# Report sequence-wise frequencies
  for (sequence_id in names(sequence_frequencies)) {
    frequencies <- sequence_frequencies[[sequence_id]]
    cat(paste(sequence_id, frequencies$A, frequencies$T, frequencies$C, frequencies$G, sep = "\t"), file = con, append = TRUE)
    cat("\n", file = con, append = TRUE)
  }
  
# Report overall frequencies
  cat("Overall\t", paste(overall_frequencies, collapse = "\t"), "\n", file = con, append = TRUE)
  
# Close the output file
  close(con)
  
  # Return the path to the output file
  return(output_path)
}

# INPUT FASTA FILE BELOW (COPY OF FILE PATH)

input_fasta_file <- "path/to/your/input.fasta"

# Process the FASTA file, transfer frequencies to output file

result <- process_fasta_file(input_fasta_file)

output_file <- write_frequencies_to_file(input_fasta_file, result$sequence, result$overall)

message("Nucleotide frequencies have been calculated and written to:", output_file, "\n")


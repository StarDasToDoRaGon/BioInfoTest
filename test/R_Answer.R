#install.packages("dplyr")
library(dplyr)

# read .rds document
setwd('E:/CSYE6200/BioInfoTest/test/')
df_R <- readRDS("exgr_test.rds")

# Task 1: How many unique transcripts are there? 
# Transcript is distinguished by transcript_id. Use unique() to get all of the unique transcript_id, and use length() to get the amount of it.
num_unique_transcript_id <- length(unique(df_R$transcript_id))
print(paste("For Task 1, the number of unique transcript is:", num_unique_transcript_id))


# Task 2: How many unique Exons are there?
# Exon is distinguished by exon_id. Use unique() to get all of the unique exon_id, and use length() to get the amount of it.
num_unique_exon_id <- length(unique(df_R$exon_id))
print(paste("For Task 2, the number of unique transcript is:", num_unique_exon_id))


# Task 3: What is the average length of an exon? What is the median length?
# Average length of exon:
mean_width <- mean(df_R$width) 
print('For Task 3:')
print(paste("Mean of width is:", mean_width))

# Median length of exon:
median_width <- median(df_R$width) 
print(paste("Median of width is:", median_width))


# Task 4: Calculate the length of intron
start_time <- Sys.time()

df_R <- df_R %>%
  group_by(transcript_id) %>%
  mutate(intron_length = ifelse(strand == "+", start - lag(end), start - lead(end)))
#df_R <- df_R %>%
#  group_by(transcript_id) %>%
#  mutate(intron_length = ifelse(strand == "+" & row_number() == 1, 0, intron_length),
#         intron_length = ifelse(strand == "-" & row_number() == n(), 0, intron_length))

end_time <- Sys.time()
 
time_taken <- end_time - start_time
print(time_taken)

df_R

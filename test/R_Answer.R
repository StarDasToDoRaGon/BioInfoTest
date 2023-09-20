#install.packages("dplyr")
library(dplyr)
library(tidyr)

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
  group_by(transcript_id) %>%       # We don't know if the data set was grouped by transcript_id, so do it before operation
  mutate(intron_length = ifelse(strand == "+",      
                                start - lag(end),          # if strand is positive: intron_length = nth start - (n-1)th end
                                start - lead(end)))%>%     # if strand is negative: intron_length = nth start - (n+1)th end
  # since group by transcript_id, the first element of positive and last element of negative are none. Those elements are just the target that should be changed to 0 
  mutate(intron_length = ifelse(is.na(intron_length), 0, intron_length))

end_time <- Sys.time() 
time_taken <- end_time - start_time      # print the time used for task 4
print(time_taken)


df_R <- df_R %>%
  mutate(l1 = start - 100,        # init 
         l2 = start + 100,
         u1 = end - 100,
         u2 = end + 100) %>%
  group_by(transcript_id, strand) %>%
  # add 2 new column to indicate if the row is a head or tail, through this col, we can decide if we should update the begin and end coordinate of l and u
  mutate(first_row = row_number() == 1,                     
         last_row = row_number() == n()) %>%
  ungroup() %>%
  # update the beginning and end coordinate
  mutate(l1 = ifelse(first_row & strand == '+', 0, l1),
         u2 = ifelse(last_row & strand == '+', 0, u2),
         u2 = ifelse(first_row & strand == '-', 0, u2),
         l1 = ifelse(last_row & strand == '-', 0, l1),
         l2 = ifelse(width < 200, (start + end) / 2.0, l2),
         u1 = ifelse(width < 200, (start + end) / 2.0, u1))

# For positive strand intron_length adjustment
df_R <- df_R %>%
  group_by(transcript_id, strand) %>%
  mutate(l1 = ifelse(strand == '+' & intron_length > 0 & intron_length < 200, as.numeric(l1 - intron_length / 2.0), l1),
         u2 = ifelse(strand == '+' & lag(intron_length) > 0 & lag(intron_length) < 200, as.numeric(lag(l1) - lag(intron_length) / 2.0), u2)) %>%
  ungroup()

# For negative strand intron_length adjustment
df_R <- df_R %>%
  group_by(transcript_id, strand) %>%
  mutate(l1 = ifelse(strand == '-' & intron_length > 0 & intron_length < 200, as.numeric(l1 - intron_length / 2.0), l1),
         u2 = ifelse(strand == '-' & lead(intron_length) > 0 & lead(intron_length) < 200, as.numeric(lead(l1) - lead(intron_length) / 2.0), u2)) %>%
  ungroup()


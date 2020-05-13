#!/usr/bin/env Rscript

# load necessary packages
library(dplyr)
library(ggplot2)
library(ggExtra)

# find all the postfilt_all files to use
postfilt_files <- Sys.glob("/home/idies/workspace/covid19/sequencing_runs/*/artic-pipeline/5-post-filter/postfilt_all.txt")
postfilt_summary <- Sys.glob("/home/idies/workspace/covid19/sequencing_runs/*/artic-pipeline/5-post-filter/postfilt_summary.txt")
validation_files <- Sys.glob("/home/idies/workspace/covid19/sequencing_runs/*/artic-pipeline/validation/postfilt_all.txt")

# run names to exclude
runs_to_exclude <- c("20200410_2018_X4_FAN32204_327837a0","20200423_2210_X1_FAN25527_ae579cee")
runs_to_exclude <- paste(runs_to_exclude, collapse = "|")

# directory in which to put plots
outdir <- "/home/idies/workspace/Storage/swohl/persistent/covid/allele_freq_validation/"

# loop through each postfilt all file
# filter out rows without illumna data at threshold
# and concatenate to large dataframe

called_vars <- data.frame()

for (postfilt_file in postfilt_files){

  # skip problematic runs
  if (grepl(runs_to_exclude,postfilt_file)){
    next
  }
  
  # read in text file
  df <- read.table(postfilt_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # filter rows without illumina support
  df <- filter(df, illumina_support != ".")
  
  # only keep position and sample name columns
  df <- df %>% select(sample,pos)
  
  # combine all files into one data frame
  called_vars <- rbind(called_vars,df)
  
}

# repeat this process for the validation files
# we need both to distinguish between called and not called variants

af_data <- data.frame()

for (val_file in validation_files){

  # skip problematic runs
  if (grepl(runs_to_exclude,val_file)){
    next
  }
  
  # read in text file
  df <- read.table(val_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # filter rows without illumina support
  df <- filter(df, illumina_support != ".")
  
  # combine all files into one data frame
  af_data <- rbind(af_data,df)
  
}

# loop through the postfilt summary files
# to determine what samples are 'no' and should not be included

sample_status <- data.frame()

for (postfilt_sum in postfilt_summary){

  # skip problematic runs
  if (grepl(runs_to_exclude,postfilt_sum)){
    next
  }
  
  # read in text file
  df <- read.table(postfilt_sum, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # keep only the relevant columns
  df <- df %>% select(Sample,Status)
  df <- df %>% rename(sample = Sample, overall_status = Status)
  
  # combine all sample data into one data frame
  sample_status <- rbind(sample_status,df)
}

# clean up the data frame

# combine allele frequency data with status data
af_data <- merge(af_data,sample_status,by = "sample",all.x = TRUE)

# exclude samples marked as no
af_data <- filter(af_data, overall_status != "No")

# make sure the values in the allele frequency columns are numeric
af_data$ont_AF = as.numeric(af_data$ont_AF)
af_data$illumina_AF = as.numeric(af_data$illumina_AF)

# add a column for low confidence flags
maybe_flags = c("depth_flag","ntc_flag","sb_flag")
af_data$maybe_flags = apply(af_data[maybe_flags],1,function(x) any(x!="."))

# add a column for called or not called
called_vars$called <- TRUE
af_data <- merge(af_data,called_vars,by=c("sample","pos"),all.x = TRUE)
af_data$called[is.na(af_data$called)] <- FALSE

# summary column for confidence and called
af_data <- af_data %>% mutate(confidence=ifelse(called==FALSE,"uncalled",
                                                ifelse(maybe_flags==TRUE,"low","high")))

# add a column for the absolute value of the difference in allele frequencies
af_data$freq_diff = af_data$ont_AF - af_data$illumina_AF

# optional filter to remove uncalled variants
af_data_called <- filter(af_data,confidence!="uncalled")

# make plots

### COLOR PLOT BY CONFIDENCE FLAGS

# plot of all variants
p_all <- ggplot(af_data, aes(x=ont_AF, y=illumina_AF, color=confidence)) + 
  geom_point() +
  #annotate("rect", xmin=-Inf, xmax=0.15, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") +
  #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.15, alpha=0.2, fill="gray") +
  geom_abline(slope=1, intercept=0, linetype = "dashed") +
  scale_color_manual(labels = c("variant", "low-confidence variant", "uncalled varaint"),values=c("blue", "skyblue2","lightgray")) +
  labs(x = "ONT allele frequency", y = "Illumina allele frequency") +
  theme_classic() +
  theme(legend.position = c(0.75, 0.1),legend.title=element_blank(),legend.background = element_blank(),
        panel.grid.major = element_line(colour="black", size = (0.2), linetype = "dotted"),
        panel.grid.minor = element_line(colour="black", size = (0.1), linetype = "dotted"))

outfile <- paste(outdir,"af_plot.png",sep="")
ggsave(outfile, ggMarginal(p_all, type="density", color="darkgray", fill="gray", alpha=0.5, size=12),
       device = "png", width = 8, height = 7)

# plot of called variants only
p_called <- ggplot(af_data_called, aes(x=ont_AF, y=illumina_AF, color=confidence)) + 
  geom_point() +
  annotate("rect", xmin=-Inf, xmax=0.15, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.15, alpha=0.2, fill="gray") +
  geom_abline(slope=1, intercept=0, linetype = "dashed") +
  scale_color_manual(labels = c("variant", "low-confidence variant", "uncalled varaint"),values=c("blue", "skyblue2","lightgray")) +
  labs(x = "ONT allele frequency", y = "Illumina allele frequency") +
  theme_classic() +
  theme(legend.position = c(0.75, 0.1),legend.title=element_blank(),legend.background = element_blank(),
        panel.grid.major = element_line(colour="black", size = (0.2), linetype = "dotted"),
        panel.grid.minor = element_line(colour="black", size = (0.1), linetype = "dotted"))

outfile <- paste(outdir,"af_plot_called.png",sep="")
ggsave(outfile, ggMarginal(p_called, type="density", color="darkgray", fill="gray", alpha=0.5, size=12),
       device = "png", width = 8, height = 7)



## PLOT THE DIFFERENCE IN ALLELE FREQUENCIES BY POSITION

p3 <- ggplot(af_data_called %>% arrange(desc(confidence)), aes(x=pos, y=freq_diff, color=confidence)) + 
  geom_point() +
  geom_hline(yintercept = 0) + 
  geom_text(aes(label=pos),hjust=0, vjust=0) +
  scale_color_manual(labels = c("variant", "low-confidence variant", "uncalled varaint"),values=c("blue", "skyblue2","lightgray")) +
  labs(x = "position", y = "ONT allele freq - Illumina allele freq") +
  theme_classic() +
  theme(legend.position = "top", legend.title=element_blank(),legend.background = element_blank(),
        panel.grid.major = element_line(colour="black", size = (0.2), linetype = "dotted"),
        panel.grid.minor = element_line(colour="black", size = (0.1), linetype = "dotted"),
        axis.line.x=element_blank())

outfile <- paste(outdir,"af_pos_plot.png",sep="")
ggsave(outfile,device = "png", width = 8, height = 5)

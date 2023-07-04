#prepare the data from fastq with awk 'NR%4 == 2 { print length($1) }' ./path/to/fastq.fastq | sort | uniq -c > read_lenght_output.txt


library(ggplot2)


# Read the output files for both samples
data1 <- read.table("read_lenght_RNA.txt", header = FALSE)
data2 <- read.table("read_lenght_DNA.txt", header = FALSE)

# Rename the columns
colnames(data1) <- c("Read_Count", "Read_Length")
colnames(data2) <- c("Read_Count", "Read_Length")

# Add a column to indicate sample identity
data1$Sample <- "vRNA"
data2$Sample <- "cDNA"

# Combine the data
combined_data <- rbind(data1, data2)


# Create the plot
ggplot(combined_data, aes(x = Read_Length, y = Read_Count, color = Sample)) +
  geom_point() +
  scale_x_log10(limits = c(20, 10000), expand = c(0, 0)) +
  xlab("Read Length") +
  ylab("Read Count") +
  theme_minimal()



ggsave("combined_read_lenght.png", width = 8, height = 7)


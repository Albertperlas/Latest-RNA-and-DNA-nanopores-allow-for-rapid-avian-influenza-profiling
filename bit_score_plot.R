
# Install and load the required package
install.packages("ggplot2")
library(ggplot2)

# Read the data from the input file
data <- read.table("average_tfm_sin_raw.txt", header = TRUE, stringsAsFactors = FALSE)

##whole genome

data$downsampling <- factor(data$downsampling, levels = c("max", "med", "min", "env"))
data$technique <- factor(data$technique, levels = c("bcftools", "ivar", "irma"))

# Create the bar chart using ggplot2
chart <- ggplot(data, aes(x = technique, y = average_bit_score, fill = technique)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(nucleotide ~ downsampling) +
  labs(x = "Technique", y = "Bit Score") +
  theme_bw(base_size = 18)

# Display the chart
print(chart)

ggsave("whole_genome.png", width = 14, height = 7)


#por segmentos 



# Create the bar chart using ggplot2
chart <- ggplot(data, aes(x = technique, y = bit_score_in_percentage, fill = technique)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
  facet_grid(segment ~ nucleotide * downsampling, scales = "fixed") +
  labs(title = "Bit Score by Technique, Nucleotide, Segment, and Downsampling",
       x = "Technique", y = "Bit Score (in %)") +
  theme_bw() +
  theme(axis.text.x = element_blank())

# Display the chart
print(chart)


#DNA plot

data$downsampling <- factor(data$downsampling, levels = c("raw", "max", "med", "min", "env"))
data$technique <- factor(data$technique, levels = c("bcftools", "ivar", "irma"))

# Create the bar chart using ggplot2
chart <- ggplot(data, aes(x = technique, y = bit_score, fill = technique)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(segment ~ nucleotide * downsampling, scales = "free") +
  labs(x = "Technique", y = "Bit Score") +
  theme_bw() +
  theme(axis.text.x = element_blank())

# Display the chart
print(chart)

ggsave("DNA.png", width = 7, height = 7)


#RNA

data$downsampling <- factor(data$downsampling, levels = c("max", "med", "min", "env"))
data$technique <- factor(data$technique, levels = c("bcftools", "ivar", "irma"))

# Create the bar chart using ggplot2
chart <- ggplot(data, aes(x = technique, y = bit_score, fill = technique)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(segment ~ nucleotide * downsampling, scales = "free") +
  labs(x = "Technique", y = "Bit Score") +
  theme_bw() +
  theme(axis.text.x = element_blank())

# Display the chart
print(chart)

ggsave("RNA.png", width = 7, height = 7)

##whole genome in percentage

options(OutDec = ",")

data$downsampling <- factor(data$downsampling, levels = c("max", "med", "min", "env"))
data$technique <- factor(data$technique, levels = c("bcftools", "ivar", "irma"))

# Create the bar chart using ggplot2
chart <- ggplot(data, aes(x = technique, y = bit_score_in_percentage, fill = technique)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(nucleotide ~ downsampling) +
  labs(x = "Technique", y = "Bit Score") +
  theme_bw(base_size = 18)

# Display the chart
print(chart)

ggsave("whole_genome.png", width = 14, height = 7)


data$downsampling <- factor(data$downsampling, levels = c("max", "med", "min", "env"))
data$technique <- factor(data$technique, levels = c("bcftools", "ivar", "irma"))

# Create the bar chart using ggplot2
chart <- ggplot(data, aes(x = technique, y = average_bit_score, fill = technique)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(nucleotide ~ downsampling) +
  labs(x = "Technique", y = "Bit Score") +
  theme_bw(base_size = 18)

# Display the chart
print(chart)

ggsave("whole_genome.png", width = 14, height = 7)




library(ggplot2)

data.frame <-read.table("read_length_vRNA.txt")
data.frame_cDNA <-read.table("read_length_cDNA.txt")

ggplot(data.frame, aes(x = data.frame$V1)) +
       geom_histogram(binwidth = 0.02,fill = 'brown', color = "brown") +
       scale_x_log10(limits = c(60, 10000), expand = c(0, 0)) +
       xlab("Read length") + ylab("Number of reads") +
       scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
       theme(plot.background = element_rect(fill = "white"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.line = element_line(color = "black"),
                         axis.text = element_text(size = 20),
                         axis.title.y = element_text(size = 20, face = "bold"),
                         axis.title.x = element_text(size = 20, face = "bold"),
                         axis.text.x = element_text(size = 20, face = "bold"),
                         axis.text.y = element_text(size = 20, face = "bold"))


ggplot(data.frame_cDNA, aes(x = data.frame$V1)) +
  geom_histogram(binwidth = 0.02,fill = 'brown', color = "brown") +
  scale_x_log10(limits = c(60, 10000), expand = c(0, 0)) +
  xlab("Read length") + ylab("Number of reads") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme(plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"))


#combined plot

ggplot() +
  geom_histogram(data = data.frame, aes(x = V1), binwidth = 0.02, fill = 'brown', color = "brown") +
  geom_histogram(data = data.frame_cDNA, aes(x = V1), binwidth = 0.02, fill = 'blue', color = "blue", alpha = 0.2) +
  scale_x_log10(limits = c(60, 10000), expand = c(0, 0)) +
  xlab("Read length") +
  ylab("Number of reads") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme(plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"))

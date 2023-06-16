library(ggplot2)
library(dplyr)

datos <- read.delim("coverage_all.txt", sep = '\t', header = TRUE, colClasses = "character")
datos$coverage <- as.numeric(datos$coverage)
datos$sample <- as.character(datos$sample)

datos$gen <- factor(datos$gen, levels = c("PB1", "PB2", "PA", "HA", "NP", "NA", "M", "NS"), ordered = TRUE)

coverage <- subset(datos, select = c(1:4))
coverage$position <- as.numeric(coverage$position)
coverage$coverage <- as.numeric(coverage$coverage)
textcol <- "grey40"

coverage$sample <- as.factor(coverage$sample)
coverage$sample <- factor(coverage$sample, levels = c("vRNA_max", "cDNA_max", "cDNA_raw", "cDNA_med", "cDNA_min", "vRNA_med", "vRNA_min", "cDNA_env","vRNA_env"), order = TRUE,
                          labels = c("vRNA_max", "cDNA_max", "cDNA_raw", "cDNA_med", "cDNA_min", "vRNA_med", "vRNA_min", "cDNA_env","vRNA_env"))

color_group <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628")

COV <- ggplot(coverage, aes(x = position, y = coverage, colour = sample)) +
  geom_point(size = 0.1) +
  geom_line() +
  scale_y_log10(breaks = c(10, 50, 250, 1000, 2500, 10000), labels = c(10, 50, 250, "1,000", "2,500", "10,000")) +
  coord_cartesian(ylim = c(1, 10000)) +
  scale_colour_manual(values = color_group) +
  geom_hline(yintercept = 50, color = "black", linetype = "dashed", size = 0.5) +
  labs(y = "Depth (log scale)", title = "Coverage, depth per gen and position", colour = "Sequenced samples") +
  facet_grid(~gen, scales = "free_x", space = "free", switch = "x") +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, face = "bold", colour = textcol),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.margin = margin(grid::unit(0, "cm")),
        legend.text = element_text(colour = textcol, size = 10, face = "bold"),
        legend.key.height = grid::unit(0.8, "cm"),
        legend.key.width = grid::unit(0.2, "cm"),
        plot.title = element_blank()) +
  guides(colour = guide_legend(nrow = 7, title.position = "top", override.aes = list(size = 3))) +
  theme(strip.text.x = element_text(size = 10),
        strip.placement = "outside",
        strip.background.x = element_rect(color = NA, fill = NA),
        strip.background.y = element_rect(color = NA, fill = NA),
        strip.text.y = element_text(size = 10, face = "bold"))

print(COV)
ggsave("coverage_all.png", width = 14, height = 7)

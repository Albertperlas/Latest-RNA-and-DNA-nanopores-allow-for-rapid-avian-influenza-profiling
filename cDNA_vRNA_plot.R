#!/bin/R
library(ggplot2)
library(dplyr)


datos<-read.delim("PATH/TO/COVERAGE/DATA/graphicR_coveragecDNA.txt", sep='\t', header=TRUE, colClasses = "character")



datos$coverage<-as.numeric(datos$coverage)

table(datos$sample)


#le decimos que la columna sample no son numeros
datos$sample<-as.character(datos$sample)
table(datos$sample)

#datos$gen<-factor(datos$gen, levels=c(levels(datos$gen), "NA"))
#datos$gen[is.na(datos$gen)] <- "NA"

#Creamos un factor para poner en orden las facetas
datos$gen<-factor(datos$gen, levels=c("PB2","PB1","PA","HA","NP","NA","M","NS"), ordered=TRUE)


#seleccionamos los datos
coverage<-subset(datos, select=c(1:4))
coverage$position<-as.numeric(coverage$position)
coverage$coverage<-as.numeric(coverage$coverage)
textcol <- "grey40"



coverage$sample<-as.factor(coverage$sample)
table(coverage$sample)


coverage$sample<-factor(coverage$sample, levels=c("vRNA", "cDNA"), order=T,
                        labels=c("vRNA", "cDNA")
)




#Creamos una palet de colores
azul<-colorRampPalette(c("#AED6F1","#1B4F72"))
nazul<- sort(unique(coverage$sample))
azulgrafico<- azul(length(nazul))

verde<-colorRampPalette(c("#A2D9CE","#0B5345"))
nverde <- sort(unique(coverage$sample[coverage$sample=="vRNA0.1"]))
verdegrafico<- verde(length(nverde))

color_group <- c("#D35400",azulgrafico, verdegrafico)

COV<-ggplot(coverage,  aes(x = position , y = coverage, colour = sample)) +
  geom_point(size=0.1) +
  geom_line()+
  scale_y_log10(breaks=c(10,50,250,1000,2500,10000), labels=c(10,50,250,"1,000","2,500","10,000")) +
  coord_cartesian(ylim = c(1 ,10000)) +
  scale_colour_manual(values=color_group) +
  geom_hline(yintercept = 50, color = "black", linetype="dashed", size=0.5)+
  labs(y="Depth (log scale)", title ="Coverage, depth per gen and position", colour="Sequenced samples") + #colour = "Samples" - title is given the aes used.
  facet_grid(~gen, scales="free_x", space="free", switch="x") +
  theme_bw() +
  theme(panel.spacing = unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=20),
        axis.title.x= element_blank(),
        axis.title.y=element_text(size=20,face="bold", colour=textcol),
        axis.ticks.x = element_blank(),
        #legend.title=element_text(size=12,face="bold", colour=textcol),
        legend.title=element_blank(),
        legend.position = "right",
        legend.justification = "left", 
        legend.direction="horizontal",
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=20,face="bold"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        plot.title = element_blank())+
  guides(colour=guide_legend(nrow = 7, title.position="top", override.aes = list(size=3)))+
  theme(strip.text.x = element_text(size = 20),
        strip.placement = "outside",
        strip.background.x=element_rect(color = NA, fill=NA),
        strip.background.y=element_rect(color = NA, fill=NA),
        strip.text.y = element_text(size = 20, face = "bold"))
last_plot()

#ggsave("coverageTODO.png", width=16, height=7)
ggsave("coverageTODOcDNA.png", width=14, height=7)



#!/bin/R
library(ggplot2)
library(dplyr)


datos<-read.delim("graphicR_coveragecDNA.txt", sep='\t', header=TRUE, colClasses = "character")



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

#fins aqui

#GRAFICO INOCULO

coverage_inoculum<-coverage[coverage$Group == "Inoculum",]
colourorange<-c("#D35400")

COVI<-ggplot(coverage,  aes(x = position , y = coverage, colour = sample)) +
  geom_point(size=0.1) +
  geom_line()+
  scale_y_log10(breaks=c(10,50,250,1000,2500,10000), labels=c(10,50,250,"1,000","2,500","10,000")) +
  coord_cartesian(ylim = c(1 ,10000)) +
  scale_colour_manual(values=colourorange) +
  geom_hline(yintercept = 50, color = "black", linetype="dashed", size=0.5)+
  labs(y="Depth (log scale)", title ="Coverage, depth per gen and position", colour="Sequenced samples") + #colour = "Samples" - title is given the aes used.
  facet_grid(~gen, scales="free_x", space="free", switch="x") +
  theme_bw() +
  theme(panel.spacing = unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y=element_text(size=10,face="bold", colour=textcol),
        axis.ticks.x = element_blank(),
        #legend.title=element_text(size=12,face="bold", colour=textcol),
        legend.title=element_text(size=10, colour=textcol,face="bold"),
        legend.position = "right",
        legend.justification = "left", 
        legend.direction="horizontal",
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=10,face="bold"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        plot.title = element_blank())+
  guides(colour=guide_legend(nrow = 1, ncol=3,title.position="top", override.aes = list(size=3)))+
  theme(strip.text.x = element_text(size = 10),
        strip.placement = "outside",
        strip.background.x=element_rect(color = NA, fill=NA),
        strip.background.y=element_rect(color = NA, fill=NA),
        strip.text.y = element_text(size = 10, face = "bold"))
last_plot()


ggsave("coverageinoculum.png", width=11, height=5)

#GRAFICO VACUNADOS E INOCULO

coverage_vaccinated<-coverage[coverage$Group == "Vaccinated",]
colour2<-c(azulgrafico)

COVV<-ggplot(coverage_vaccinated,  aes(x = position , y = coverage, colour = sample)) +
  geom_point(size=0.1) +
  geom_line()+
  scale_y_log10(breaks=c(10,50,250,1000,2500,10000), labels=c(10,50,250,"1,000","2,500","10,000")) +
  coord_cartesian(ylim = c(1 ,10000)) +
  scale_colour_manual(values=colour2) +
  geom_hline(yintercept = 50, color = "black", linetype="dashed", size=0.5)+
  labs(y="Depth (log scale)", title ="Coverage, depth per gen and position", colour="Sequenced samples") + #colour = "Samples" - title is given the aes used.
  facet_grid(~gen, scales="free_x", space="free", switch="x") +
  theme_bw() +
  theme(panel.spacing = unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y=element_text(size=10,face="bold", colour=textcol),
        axis.ticks.x = element_blank(),
        #legend.title=element_text(size=12,face="bold", colour=textcol),
        legend.title=element_text(size=10, colour=textcol,face="bold"),
        legend.position = "right",
        legend.justification = "left", 
        legend.direction="horizontal",
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=10,face="bold"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        plot.title = element_blank())+
  guides(colour=guide_legend(ncol = 3, title.position="top", override.aes = list(size=3)))+
  theme(strip.text.x = element_text(size = 10),
        strip.placement = "outside",
        strip.background.x=element_rect(color = NA, fill=NA),
        strip.background.y=element_rect(color = NA, fill=NA),
        strip.text.y = element_text(size = 10, face = "bold"))
last_plot()


ggsave("coverageVAcinoculum.png", width=11, height=5)

#GRAFICO NO VACUNADOS

coverage_unvaccinated<-coverage[coverage$Group == "Nonvaccinated",]
colour3<-c(verdegrafico)

COVN<-ggplot(coverage_unvaccinated,  aes(x = position , y = coverage, colour = sample)) +
  geom_point(size=0.1) +
  geom_line()+
  scale_y_log10(breaks=c(10,50,250,1000,2500,10000), labels=c(10,50,250,"1,000","2,500","10,000")) +
  coord_cartesian(ylim = c(1 ,10000)) +
  scale_colour_manual(values=colour3) +
  geom_hline(yintercept = 50, color = "black", linetype="dashed", size=0.5)+
  labs(y="Depth (log scale)", title ="Coverage, depth per gen and position", colour="Sequenced samples") + #colour = "Samples" - title is given the aes used.
  facet_grid(~gen, scales="free_x", space="free", switch="x") +
  theme_bw() +
  theme(panel.spacing = unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y=element_text(size=10,face="bold", colour=textcol),
        axis.ticks.x = element_blank(),
        #legend.title=element_text(size=12,face="bold", colour=textcol),
        legend.title=element_text(size=10, colour=textcol,face="bold"),
        legend.position = "right",
        legend.justification = "left", 
        legend.direction="horizontal",
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=10,face="bold"),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        plot.title = element_blank())+
  guides(colour=guide_legend(ncol = 3, title.position="top", override.aes = list(size=3)))+
  theme(strip.text.x = element_text(size = 10),
        strip.placement = "outside",
        strip.background.x=element_rect(color = NA, fill=NA),
        strip.background.y=element_rect(color = NA, fill=NA),
        strip.text.y = element_text(size = 10, face = "bold"))
last_plot()


ggsave("coveragenovac.png", width=11, height=5)


library(ggpubr)
TODOCOV<-ggarrange(COVI, COVV, COVN , 
                   #labels = c("A.1", "A.2", "A.3"),
                   ncol = 1, nrow = 3, align = "v")

ggsave("coverageTODO.png", width=12, height=7)




#HeatMap Plot de profundidad media

#HeatMap Plot de profundidad media
datos$sample<-as.character(datos$sample)
datos$coverage<- as.numeric(datos$coverage)

fragment_media<-aggregate(coverage~sample+gen+Group, datos, median)


fragment_media_2 <- fragment_media %>%
  mutate (sample=factor(sample, levels=rev(sort(unique(sample))))) %>%
  mutate(coveragefactor=cut(coverage, breaks = c(-1,50,100,1000,max(coverage, na.rm=T)), labels=c("< 50","50-100","> 100","> 1,000"))) %>%
  mutate(coveragefactor=factor(as.character(coveragefactor), levels=rev(levels(coveragefactor))))

levels(fragment_media_2$sample)
#assign text colour

fragment_media_2$Group<-factor(fragment_media_2$Group, levels=c("Inoculum","Vaccinated", "Nonvaccinated"), labels=c("","Vaccinated", "Nonvaccinated"))


fragment_media_2$sample<-factor(fragment_media_2$sample, levels=c("Inoculum",
                                                                  "01BALF",
                                                                  "02BALF", 
                                                                  "04BALF",
                                                                  "053dpi" ,"054dpi", "055dpi", "05BALF",
                                                                  "064dpi", "06BALF", 
                                                                  "073dpi", "074dpi", 
                                                                  "084dpi",
                                                                  "09BALF",
                                                                  "10BALF", 
                                                                  "11BALF", 
                                                                  "123dpi", "124dpi", "125dpi", "12BALF",
                                                                  "133dpi", "134dpi", "135dpi", "13BALF",
                                                                  "143dpi", "144dpi", "14BALF", 
                                                                  "153dpi", "155dpi",
                                                                  "163dpi", "164dpi", "165dpi", "169dpi" ), order=T,
                                labels=c("Inoculum",
                                         "1 BALF",
                                         "2 BALF", 
                                         "4 BALF",
                                         "5 3 dpi" ,"5 4 dpi", "5 5 dpi", "5 BALF",
                                         "6 4 dpi", "6 BALF", 
                                         "7 3 dpi", "7 4 dpi", 
                                         "8 4 dpi",
                                         "9 BALF",
                                         "10 BALF", 
                                         "11 BALF", 
                                         "12 3 dpi", "12 4 dpi", "12 5 dpi", "12 BALF",
                                         "13 3 dpi", "13 4 dpi", "13 5 dpi", "13 BALF",
                                         "14 3 dpi", "14 4 dpi", "14 BALF", 
                                         "15 3 dpi", "15 5 dpi",
                                         "16 3 dpi", "16 4 dpi", "16 5 dpi", "16 9 dpi")
)



fragment_green<-ggplot(fragment_media_2,  aes(x = gen, y = sample, fill=coveragefactor)) +
  geom_tile(colour="black", size=0.1)+
  labs(x="",y="", title = "Median depth per sequenced fragment")+
  guides(fill=guide_legend(title="Median depth value"), colour=textcol)+
  scale_y_discrete(limits=rev)+
  scale_fill_manual(values=c("#145A32","#1E8449","#27AE60","#A9DFBF"),
                    na.value = "grey90")+
  theme_classic() +
  facet_grid(Group~., scales = "free_y", space = "free_y", switch = "y")+
  theme(plot.title = element_blank(),
        #plot.title = element_text(hjust = 0,1, size=12, face="bold",colour=textcol),
        axis.line=element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, face="bold"),
        axis.ticks = element_blank(),
        panel.spacing =unit(0.2, "lines"),
        strip.placement = "outside",
        strip.background.x = element_blank(),
        strip.background.y=element_rect(color = NA),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 13, face = "bold"),
        legend.title=element_text(size=10,face="bold", colour=textcol),
        legend.position = "right",
        legend.direction="vertical",
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=10,face="bold"),
        legend.key.height=grid::unit(0.6,"cm"),
        legend.key.width=grid::unit(0.2,"cm"))

ggsave("fragmentgreenmedian.png", width=7, height=7)

fragment_blue<-ggplot(fragment_media_2,  aes(x = gen, y = sample, fill=coveragefactor)) +
  geom_tile(colour="black", size=0.1)+
  labs(x="",y="", title = "Median depth per sequenced segment")+
  scale_y_discrete(limits=rev)+
  guides(fill=guide_legend(title="Median depth value"), colour=textcol)+
  scale_fill_manual(values=c("#1B4F72","#2874A6","#3498DB","#85C1E9"),
                    na.value = "grey90")+
  theme_classic() +
  facet_grid(Group~., scales = "free_y", space = "free_y", switch = "y")+
  theme(plot.title = element_blank(),
        #plot.title = element_text(hjust = 0,1, size=12, face="bold",colour=textcol),
        axis.line=element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, face="bold"),
        axis.ticks = element_blank(),
        panel.spacing =unit(0.2, "lines"),
        strip.placement = "outside",
        strip.background.x = element_blank(),
        strip.background.y=element_rect(color = NA),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 13, face = "bold"),
        legend.title=element_text(size=10,face="bold", colour=textcol),
        legend.position = "right",
        legend.direction="vertical",
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=10,face="bold"),
        legend.key.height=grid::unit(0.6,"cm"),
        legend.key.width=grid::unit(0.2,"cm"))

ggsave("fragmentbluemedian.png", width=7, height=7)

fragment_orange<-ggplot(fragment_media_2,  aes(x = gen, y = sample, fill=coveragefactor)) +
  geom_tile(colour="black", size=0.1)+
  labs(x="",y="", title = "Median depth per sequenced fragment")+
  guides(fill=guide_legend(title="Median depth value"), colour=textcol)+
  scale_fill_manual(values=c("#873600","#BA4A00","#DC7633","#EDBB99"),
                    na.value = "grey90")+
  theme_classic() +
  facet_grid(Group~., scales = "free_y", space = "free_y", switch = "y")+
  theme(plot.title = element_blank(),
        #plot.title = element_text(hjust = 0,1, size=12, face="bold",colour=textcol),
        axis.line=element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, face="bold"),
        axis.ticks = element_blank(),
        panel.spacing =unit(0.2, "lines"),
        strip.placement = "outside",
        strip.background.x = element_blank(),
        strip.background.y=element_rect(color = NA),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 13, face = "bold"),
        legend.title=element_text(size=10,face="bold", colour=textcol),
        legend.position = "right",
        legend.direction="vertical",
        legend.margin=margin(grid::unit(0,"cm")),
        legend.text=element_text(colour=textcol,size=10,face="bold"),
        legend.key.height=grid::unit(0.6,"cm"),
        legend.key.width=grid::unit(0.2,"cm"))

ggsave("fragmentredmedian.png", width=7, height=7)


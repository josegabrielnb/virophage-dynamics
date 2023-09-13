library(ggplot2)
library(gridExtra)

#plot example of time course of infection

setwd("~/Desktop/JoseGabriel/PLoS_CompBio/r_figures/")

#simulation random seed 1234

### Virophages

#list.files(".")

cells1_neutral <- read.csv("cvv_dynamics_cells1_neutral.csv")
cells1_inhibition <- read.csv("cvv_dynamics_cells1_inhibition.csv")
cells1_pcd <- read.csv("cvv_dynamics_cells1_pcd.csv")

#View(cells1_neutral)

colours <- c("Neutral virophages" = "#009E73", "Inhibitory virophages" = "#D55E00", "PCD" = "#56B4E9")

#Single-celled organisms
p1 <- ggplot()+
  geom_line(aes(x=time,y=virophages,colour="Neutral virophages"),data=cells1_neutral)+
  geom_line(aes(x=time,y=virophages,colour="Inhibitory virophages"),data=cells1_inhibition)+
  geom_line(aes(x=time,y=virophages,colour="PCD"),data=cells1_pcd)+
  ylab("Virophage population")+
  xlab("Time step")+
  ggtitle("Single-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = c(.7,.7),legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"))+
  scale_color_manual(values = colours)+
  labs(colour="Simulation")+
  ylim(c(0,135000))

#2-celled organisms

cells2_neutral <- read.csv("cvv_dynamics_cells2_neutral.csv")
cells2_inhibition <- read.csv("cvv_dynamics_cells2_inhibition.csv")
cells2_pcd <- read.csv("cvv_dynamics_cells2_pcd.csv")

p2 <- ggplot()+
  geom_line(aes(x=time,y=virophages,colour="Neutral virophages"),data=cells2_neutral)+
  geom_line(aes(x=time,y=virophages,colour="Inhibitory virophages"),data=cells2_inhibition)+
  geom_line(aes(x=time,y=virophages,colour="PCD"),data=cells2_pcd)+
  ylab("Virophage population")+
  xlab("Time step")+
  #ylim(c(0,45000))+
  ggtitle("2-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = c(.7,.7),legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"))+
  scale_color_manual(values = colours)+
  labs(colour="Simulation")+
  ylim(c(0,135000))

#4-celled organisms

cells4_neutral <- read.csv("cvv_dynamics_cells4_neutral.csv")
cells4_inhibition <- read.csv("cvv_dynamics_cells4_inhibition.csv")
cells4_pcd <- read.csv("cvv_dynamics_cells4_pcd.csv")

p3 <- ggplot()+
  geom_line(aes(x=time,y=virophages,colour="Neutral virophages"),data=cells4_neutral)+
  geom_line(aes(x=time,y=virophages,colour="Inhibitory virophages"),data=cells4_inhibition)+
  geom_line(aes(x=time,y=virophages,colour="PCD"),data=cells4_pcd)+
  ylab("Virophage population")+
  xlab("Time step")+
  ggtitle("4-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = c(.7,.7),legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"))+
  scale_color_manual(values = colours)+
  labs(colour="Simulation")+
  ylim(c(0,135000))

#8-celled organisms

cells8_neutral <- read.csv("cvv_dynamics_cells8_neutral.csv")
cells8_inhibition <- read.csv("cvv_dynamics_cells8_inhibition.csv")
cells8_pcd <- read.csv("cvv_dynamics_cells8_pcd.csv")

p4 <- ggplot()+
  geom_line(aes(x=time,y=virophages,colour="Neutral virophages"),data=cells8_neutral)+
  geom_line(aes(x=time,y=virophages,colour="Inhibitory virophages"),data=cells8_inhibition)+
  geom_line(aes(x=time,y=virophages,colour="PCD"),data=cells8_pcd)+
  ylab("Virophage population")+
  xlab("Time step")+
  ggtitle("8-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = c(.7,.7),legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"))+
  scale_color_manual(values = colours)+
  labs(colour="Simulation")+
  ylim(c(0,135000))

grid.arrange(p1,p2,p3,p4,nrow=2)

### Viruses

#Single-celled organisms
p5 <- ggplot()+
  geom_line(aes(x=time,y=viruses,colour="Neutral virophages"),data=cells1_neutral)+
  geom_line(aes(x=time,y=viruses,colour="Inhibitory virophages"),data=cells1_inhibition)+
  geom_line(aes(x=time,y=viruses,colour="PCD"),data=cells1_pcd)+
  ylab("Virus population")+
  xlab("Time step")+
  ggtitle("Single-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = c(.7,.7),legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"))+
  scale_color_manual(values = colours)+
  labs(colour="Simulation")+
  ylim(c(0,14000))

#2-celled organisms
p6 <- ggplot()+
  geom_line(aes(x=time,y=viruses,colour="Neutral virophages"),data=cells2_neutral)+
  geom_line(aes(x=time,y=viruses,colour="Inhibitory virophages"),data=cells2_inhibition)+
  geom_line(aes(x=time,y=viruses,colour="PCD"),data=cells2_pcd)+
  ylab("Virus population")+
  xlab("Time step")+
  ylim(c(0,4000))+
  ggtitle("2-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = c(.7,.7),legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"))+
  scale_color_manual(values = colours)+
  labs(colour="Simulation")+
  ylim(c(0,14000))

#4-celled organisms
p7 <- ggplot()+
  geom_line(aes(x=time,y=viruses,colour="Neutral virophages"),data=cells4_neutral)+
  geom_line(aes(x=time,y=viruses,colour="Inhibitory virophages"),data=cells4_inhibition)+
  geom_line(aes(x=time,y=viruses,colour="PCD"),data=cells4_pcd)+
  ylab("Virus population")+
  xlab("Time step")+
  ggtitle("4-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = c(.7,.7),legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"))+
  scale_color_manual(values = colours)+
  labs(colour="Simulation")+
  ylim(c(0,14000))

#8-celled organisms
p8 <- ggplot()+
  geom_line(aes(x=time,y=viruses,colour="Neutral virophages"),data=cells8_neutral)+
  geom_line(aes(x=time,y=viruses,colour="Inhibitory virophages"),data=cells8_inhibition)+
  geom_line(aes(x=time,y=viruses,colour="PCD"),data=cells8_pcd)+
  ylab("Virus population")+
  xlab("Time step")+
  ggtitle("8-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = c(.7,.7),legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"))+
  scale_color_manual(values = colours)+
  labs(colour="Simulation")+
  ylim(c(0,14000))

grid.arrange(p5,p6,p7,p8,nrow=2)

### Cells

#Single-celled organisms
p9 <- ggplot()+
  geom_line(aes(x=time,y=cells,colour="Neutral virophages"),data=cells1_neutral)+
  geom_line(aes(x=time,y=cells,colour="Inhibitory virophages"),data=cells1_inhibition)+
  geom_line(aes(x=time,y=cells,colour="PCD"),data=cells1_pcd)+
  ylab("Cell population")+
  xlab("Time step")+
  ylim(c(250,1025))+
  ggtitle("Single-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = "bottom",legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"), legend.margin=margin(t=-3))+
  scale_color_manual(values = colours)+
  labs(colour="")

#2-celled organisms
p10 <- ggplot()+
  geom_line(aes(x=time,y=cells,colour="Neutral virophages"),data=cells2_neutral)+
  geom_line(aes(x=time,y=cells,colour="Inhibitory virophages"),data=cells2_inhibition)+
  geom_line(aes(x=time,y=cells,colour="PCD"),data=cells2_pcd)+
  ylab("Cell population")+
  xlab("Time step")+
  ylim(c(250,1025))+
  ggtitle("2-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = "bottom",legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"), legend.margin=margin(t=-3))+
  scale_color_manual(values = colours)+
  labs(colour="")

#4-celled organisms
p11 <- ggplot()+
  geom_line(aes(x=time,y=cells,colour="Neutral virophages"),data=cells4_neutral)+
  geom_line(aes(x=time,y=cells,colour="Inhibitory virophages"),data=cells4_inhibition)+
  geom_line(aes(x=time,y=cells,colour="PCD"),data=cells4_pcd)+
  ylab("Cell population")+
  xlab("Time step")+
  ylim(c(250,1025))+
  ggtitle("4-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = "bottom",legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"), legend.margin=margin(t=-3))+
  scale_color_manual(values = colours)+
  labs(colour="")

#8-celled organisms
p12 <- ggplot()+
  geom_line(aes(x=time,y=cells,colour="Neutral virophages"),data=cells8_neutral)+
  geom_line(aes(x=time,y=cells,colour="Inhibitory virophages"),data=cells8_inhibition)+
  geom_line(aes(x=time,y=cells,colour="PCD"),data=cells8_pcd)+
  ylab("Cell population")+
  xlab("Time step")+
  ylim(c(250,1025))+
  ggtitle("8-celled organisms")+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5), 
        axis.title=element_text(size=12), axis.text = element_text(size=10))+
  theme(legend.position = "bottom",legend.text=element_text(size=10), legend.title=element_text(size=10,face="bold"), legend.margin=margin(t=-3))+
  scale_color_manual(values = colours)+
  labs(colour="")

grid.arrange(p9,p10,p11,p12,nrow=2)

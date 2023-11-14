## Combining all Species data
library(readxl)
library(readr)
data1 <- read_csv("aurea_variables.csv")
data2 <- read_csv("cyna_variables.csv")
data3 <- read_excel("PRMUND.xlsx")
data4 <- read_excel("PRREPE.xlsx")
Protea_Cyna <- read_csv("Protea_Cyna.csv")
Protea_Aurea <- read_csv("Protea_Aurea.csv")
data4 <- subset(data4, select = c(-Lats, -PRO))
data3 <- subset(data3, select = c(-Lats, -PRO))
data1 <- subset(data1,select  = c(-...1))
data2 <- subset(data2,select  = c(-...1))
data1 <- cbind(Protea_Aurea$OBS,data1)
colnames(data1) <- c("OBS","latitude" ,"longitude","july_min","jan_max","prec","evap","soil","wind" )
data2 <- cbind(Protea_Cyna$OBS,data2)
colnames(data2) <- c("OBS","latitude" ,"longitude","july_min","jan_max","prec","evap","soil","wind" )
data3 <- data3[,c(1,2,3,5,4,6,7,8,9)]
data4 <- data4[,c(1,2,3,5,4,6,7,8,9)]
colnames(data3) <- c("OBS","latitude" ,"longitude","july_min","jan_max","prec","evap","soil","wind" )
colnames(data4) <- c("OBS","latitude" ,"longitude","july_min","jan_max","prec","evap","soil","wind" )
data5 <- read_csv("Species_LDSG.csv")
data6 <- read_csv("Species_PRPUNC.csv")
colnames(data5) <- c("OBS","latitude" ,"longitude","july_min","jan_max","prec","evap","soil","wind" )
colnames(data6) <- c("OBS","latitude" ,"longitude","july_min","jan_max","prec","evap","soil","wind" )
Complete_data <- rbind(data1,data2,data3,data4,data5,data6)

### Plotting Graphs 
library(ggmap)
library(ggplot2)
library(viridis)
library(maps)
world = map_data('world')
south_africa = subset(world, region == "South Africa")

a <-  ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = Complete_data ,aes(x = longitude  , y = latitude, color = july_min ),size = 2) + 
  labs(title = "July (Winter)Minimum Temperature", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in"),
        panel.border = element_rect(size = 6, color = "black") # Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31)) +
  scale_fill_viridis_c(guide = guide_colorbar(title = NULL)) +
  theme_bw()+
  scale_color_viridis_c(
    guide = guide_colorbar(title = 'July Min'),
    limits = c(min(Complete_data$wind),quantile(Complete_data$wind,0.95)))

b <-  ggplot() +
  geom_polygon(data = south_africa,aes(x = long, y = lat, group = group), fill = "white", color = 'black', size = 1) + 
  geom_point(data = Complete_data, aes(x = longitude, y = latitude, color = jan_max), size = 2) +
  labs(title = "January (Summer) Maximum Temperature", x = "Longitude", y = "Latitude") + 
  theme_minimal() + 
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.ticks = element_line(size = 1),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.key.size = unit(20, "in")  # Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31)) +
  theme_bw() +
  scale_color_viridis_c(
    guide = guide_colorbar(title = 'Jan Max'),
    limits = c(min(Complete_data$wind),quantile(Complete_data$wind,0.95)))

c <-  ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = Complete_data ,aes(x = longitude  , y = latitude, color = prec ),size = 2) + 
  labs(title = "Precipitation", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in"),
        panel.border = element_rect(size = 6, color = "black") # Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31)) +
  scale_fill_viridis_c(guide = guide_colorbar(title = NULL)) +
  theme_bw()+
  scale_color_viridis_c(
    guide = guide_colorbar(title = 'Prec'),
    limits = c(min(Complete_data$wind),quantile(Complete_data$wind,0.95)))

d <-  ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = Complete_data ,aes(x = longitude  , y = latitude, color = evap ),size = 2) + 
  labs(title = "Potential Evapotranspiration", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in"),
        panel.border = element_rect(size = 6, color = "black") # Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31)) +
  scale_fill_viridis_c(guide = guide_colorbar(title = NULL)) +
  theme_bw()+
  scale_color_viridis_c(
    guide = guide_colorbar(title = 'Evap'),
    limits = c(min(Complete_data$wind),quantile(Complete_data$wind,0.95)))

e <-  ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = Complete_data,aes(x = longitude  , y = latitude, color = soil ),size = 2) + 
  labs(title = "SOIL ", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in"),
        panel.border = element_rect(size = 6, color = "black") # Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31)) +
  scale_fill_viridis_c(guide = guide_colorbar(title = NULL)) +
  theme_bw()+
  scale_color_viridis_c(
    guide = guide_colorbar(title = 'Soil'),
    limits = c(min(Complete_data$wind),quantile(Complete_data$wind,0.95)))

f <- ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = Complete_data,aes(x = longitude  , y = latitude, color = wind ),size = 2) + 
  labs(title = "WIND ", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in")# Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31))+
  theme_bw()+
  scale_color_viridis_c(
    guide = guide_colorbar(title = 'Wind'),
    limits = c(min(Complete_data$wind),quantile(Complete_data$wind,0.95))
  )
## Maps for all six covariate surfaces over the CFR
library(gridExtra)
maps = grid.arrange(a,b,c,d,e,f,nrow = 3)


###### Plot for Presence of Different Species 
## Species Protea_Aurea
a = ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = data1 ,aes(x = longitude  , y = latitude),color = "red",size = 3) + 
  labs(title = "Species Protea_Aurea ", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in")# Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31))+
  theme_bw()
## Species Protea_Cyna
b = ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = data2 ,aes(x = longitude  , y = latitude),color = "orange",size = 3) + 
  labs(title = "Species Protea_Cyna ", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in")# Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31))+
  theme_bw()
## Species PRMUND
c = ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = data3 ,aes(x = longitude  , y = latitude),color = "darkblue",size = 3) + 
  labs(title = "Species PRMUND ", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in")# Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31))+
  theme_bw()

## Species PRRPED
d = ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = data4 ,aes(x = longitude  , y = latitude),color = "darkgreen",size = 3) + 
  labs(title = "Species PRRPED ", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in")# Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31))+
  theme_bw()

## Species LDSG

e = ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = data5 ,aes(x = longitude  , y = latitude),color = "brown",size = 3) + 
  labs(title = "Species LDSG", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in")# Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31))+
  theme_bw()

## Species PRPUNC

f = ggplot()+
  geom_polygon(data = south_africa, aes(x = long, y = lat, group = group),fill = "white",color = 'black', size = 1) + 
  geom_point(data = data6 ,aes(x = longitude  , y = latitude),color = "darkviolet",size = 3) + 
  labs(title = "Species PRPUNC", x = "Longitude", y = "Latitude") + 
  theme_minimal()+ 
  theme(axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        axis.ticks = element_line(size = 1) ,
        plot.title = element_text(size=20, hjust = 0.5),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(20, "in")# Adjust the size and color of the panel border
  ) +
  coord_cartesian(xlim = c(18, 27), ylim = c(-35, -31))+
  theme_bw()

## Maps for presence of species on location
library(gridExtra)
maps = grid.arrange(a,b,c,d,e,f,nrow = 3)

##############################################################
#Trace Plot
##############################################################

data_1 = read.csv('beta_samples_aurea_final.csv')
data_1 = data_1[0:750,]
data_2 = read.csv('beta_samples_aurea_1.csv')
data_2 = data_2[0:750,]
data_3 = read.csv('beta_samples_aurea_2.csv')
data_3 = data_3[0:750,]
par(mfrow=c(3,1))
tp1 = traceplot(mcmc(data_1[,3]),main='Trace Plot for first initial choice',ylab = 'Posterior beta samples')
tp2 = traceplot(mcmc(data_2[,3]),main='Trace Plot for second initial choice',ylab = 'Posterior beta samples')
tp3 = traceplot(mcmc(data_3[,3]),main='Trace Plot for third initial choice',ylab = 'Posterior beta samples')

tp = grid.arrange(tp1,tp2,tp3,nrow = 2)

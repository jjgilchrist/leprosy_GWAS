rm(list = ls(all=TRUE))
library(Rcpp)

library("ggplot2")
theme_set(theme_bw())
library("sf")

library("rnaturalearth")
library("rnaturalearthdata")
library("RColorBrewer")
library("ggrepel")

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(11,"RdBu")
africa <- ne_countries(scale = "medium", returnclass = "sf", continent = "africa")

cities <- data.frame(c("Bamako", "Karonga"),
                     c(-8.0029,33.9248),
                     c(12.6392,-9.9525))
colnames(cities) <- c("city", "lng", "lat")
cities <- st_as_sf(cities, coords = c("lng", "lat"), remove = FALSE, crs = 4326, agr = "constant")

map <- ggplot(data = africa) +
  geom_sf(aes(fill = geounit)) +
  scale_fill_manual(values = c(cols2[c(rep(6,28),3,9,rep(6,24))])) +
  theme(legend.position = "none") +
  geom_sf(data = cities) +
  geom_text_repel(data = cities, aes(x = lng, y = lat, label = city), 
                  fontface = "bold", nudge_x = c(20, -15), nudge_y = c(0.5, -5)) +
  ylab("Latitude") +
  xlab("Longitude")

map

getwd()

ggsave(
  "Fig1A.jpg",
  width = 4.5,
  height = 4.5,
  dpi = 300
) 


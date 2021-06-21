#######################################################
# COMPUTE AND MAP SHORTEST DRIVING DISTANCE           #
# Milos Popovic                                       #
# 2021-06-20                                          #
#######################################################

     # 	|￣￣￣￣￣￣|
    # 	| BEGIN      | 
   # 	| REPLICATION| 
  # 	|            |
 # 	|            | 
  # 	| ＿＿＿＿＿__| 
   # 	(\__/) || 
    # 	(•ㅅ•) || 
     # 	/ 　 づ                                                                                    

# install libraries if needed
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("osmdata")) install.packages("osmdata")
if(!require("sf")) install.packages("sf")
if(!require("dodgr")) remotes::install_git("https://git.sr.ht/~mpadge/dodgr")
if(!require("geosphere")) install.packages("geosphere")
if(!require("classInt")) install.packages("classInt")
if(!require("extrafont")) install.packages("extrafont")
if(!require("ggmap")) install.packages("ggmap")

# load libraries
library(tidyverse, quietly=T) # data processing
library(osmdata, quietly=T) # load osm data
library(sf, quietly=T) # use spatial vector data 
library(dodgr, quietly=T) # driving distance
library(geosphere, quietly=T) # aerial distance
library(classInt, quietly=T) # legend
library(extrafont, quietly=T) # font

set.seed(20210618)

#download official 2021 Serbian census circles
u <- "https://github.com/justinelliotmeyers/Official_Serbia_2021_Administrative_Boundaries/raw/main/popisni_krug-gpkg.zip"
download.file(u, basename(u), mode="wb")
unzip("popisni_krug-gpkg.zip")

#load census circles
pk <- st_read(paste0(getwd(), "/tmp/data/ready/pk/", "popisni_krug.gpkg"), stringsAsFactors = FALSE) %>% 
        st_transform(4326) %>% 
        st_as_sf()

# define Belgrade's bounding box based on ggmap's bounding box
bg_map <- ggmap::get_map(getbb("Belgrade"), 
	maptype = "toner-lite", 
	source = "stamen", 
	color="bw", 
	force=T)

bg_bbox <- attr(bg_map, 'bb')

bbox <-c(xmin=bg_bbox[,2], 
	ymin= bg_bbox[,1], 
	xmax= bg_bbox[,4],
	ymax=bg_bbox[,3])

#filter Belgrade
pkb <- st_crop(pk, bbox)
plot(pkb["objectid"])

#create centroids and their unique identifier
cen <- st_centroid(pkb)
cen$id <- 1:max(nrow(cen))

#fetch Belgrade's petrol stations
bg_amen <- opq(bbox = bbox, timeout = 180, memsize = 104857600) %>%
      add_osm_feature(
      key = 'amenity', 
      value = "fuel"  
    ) %>% 
  osmdata_sf(quiet = FALSE)
bg_pts <- bg_amen[c("osm_points")] #filter only points
bg_p <-  do.call(rbind, bg_pts) %>% #turn from list to sf object
		  select("osm_id", "geometry")

#get Belgrade's paved roads
bg_way <- opq(bbox = bbox, timeout = 120, memsize = 104857600) %>%
      add_osm_feature(
      key = 'highway') %>% 
  osmdata_sf(quiet = T)

bg_r <- bg_way$osm_lines

# let's peek into the locations of ATMs/banks and roads
ggmap::ggmap(bg_map)+
  geom_sf(data = bg_r,
          inherit.aes = F,
          col = "#3875D9",
          size = .15)+
  geom_sf(data = bg_p,
          inherit.aes = F,
          col = "#d94496",
          alpha = .6,
          size = 1.5)+
  theme_void()

# decompose roads into distinct edges
g <- weight_streetnet(bg_r, wt_profile = "motorcar", type_col = "highway")
dim(g)
head(g)

# define origin, destination and compute distance
from <- st_coordinates(cen)
to <- st_coordinates(bg_p)
d <- dodgr_dists(graph = g, from = from, to = to)
df <- apply(d, 1, FUN=min, na.rm=TRUE) %>% as.data.frame()
b <- st_sf(data.frame(pkb, df))
names(b)[14] <- "dist"

summary(b$dist)
nrow(subset(b, dist=="Inf"))

# calculate the aerial distance from every centroid to the closest station
b$id <- 1:max(nrow(b)) # create id
binf <- subset(b, dist=="Inf") #subset rows with infinite values
cinf <- subset(cen, id%in%binf$id) #filter coords by rows with infinite values
c <- st_coordinates(cinf) # get coordinates
c <- as.data.frame(cbind(c, cinf$id)) #convert into df
names(c) <- c("long", "lat", "id") #choose intuitive names

# function to find the shortest aerial distance
min_dist <- function(loc){
 sd <- c[c$id==loc,]
 sd1 <- distGeo(sd[,1:2], to[,1:2])
 sd2 <- data.frame(id = loc, 
 				   dist=min(sd1))
 return(sd2)
}
dist_mat <- dplyr::bind_rows(lapply(c$id, min_dist))

# plug new distances back into "b" and create "bb"
# use if_else to replace infinite values with aerial distance
bb <- b %>% 
  left_join(dist_mat, by = "id") %>% 
  mutate(dist = if_else(is.infinite(dist.x), dist.y, dist.x))

# let's find a natural interval with quantile breaks
ni = classIntervals(bb$dist, 
           n = 8, 
           style = 'quantile')$brks

# this function uses above intervals to create categories
labels <- c()
for(i in 1:length(ni)){
    labels <- c(labels, paste0(round(ni[i], 0), 
                             "–", 
                             round(ni[i + 1], 0)))
}
labels <- labels[1:length(labels)-1]

# finally, carve out the categorical variable based on the breaks and labels
bb$cat <- cut(bb$dist, 
              breaks = ni, 
              labels = labels, 
              include.lowest = T)
levels(bb$cat) # let's check how many levels it has (8)

# plot
p <- ggplot() +
geom_sf(data=bb, aes(fill = cat), color=NA, size=0) +
    coord_sf(crs = 4326, datum = NA) +
  scale_fill_manual(name= "meters", 
    		values = c('#ffffca', '#b9e0ad', '#72c099', '#109c99', 
    		   '#117581', '#154f68', '#122c4e', '#0c0636'),
    		labels = c("0–427",      "427–601",    "601–755",    
    				   "755–903",    "903–1066", "1066–1288",
    				   "1288–1758",  ">1758"),
    		drop = F)+
  guides(fill = guide_legend(
         direction = "horizontal",
         keyheight = unit(1.15, units = "mm"),
         keywidth = unit(20, units = "mm"),
         title.position = 'top',
         title.hjust = 0.5,
         label.hjust = .5,
         nrow = 1,
         byrow = T,
         reverse = F,
         label.position = "bottom"
          )) +
    theme_minimal() +
  theme(text=element_text(family="Georgia"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size=9, color="white", hjust=.7, vjust=200),
    axis.title.y = element_blank(),
    legend.position = c(.5, -.015),
    legend.text = element_text(size=10, color="grey20"),
    legend.title = element_text(size=11, color="grey20"),
    panel.grid.major = element_line(color = "white", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin     =   unit(c(t=0, r=0, b=0, l=0),"lines"), #added these narrower margins to enlarge map
    plot.title = element_text(face="bold", size=17, color="#095169", hjust=.5, vjust=-2),
    plot.subtitle = element_text(size=16, color="#53ba83", hjust=.5, vjust=-2),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.background = element_rect(fill = "white", color = NA), 
    legend.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank()) +
  labs(x = "©2021 Milos Popovic https://milospopovic.net\n Data: OSM Geofabrik",
    title = "Shortest driving distance to a petrol station in Belgrade", 
    subtitle = "at census circle level", 
    caption = "")

#	|￣￣￣￣￣￣ |
 #    	| END        | 
  #    	| REPLICATION| 
   #   	|            |
    #   |            | 
   #    | ＿＿＿＿＿__| 
  #    	(\__/) || 
#    	(•ㅅ•) || 
#    	/ 　 づ                                                                                    

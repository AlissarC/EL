
################## Whitakker diagram to make Figure S1  ############################ 
library(ggplot2)
library(plotbiomes)

whittaker_plot <- whittaker_base_plot() 

whittaker_plot <- whittaker_base_plot() +

    theme(legend.position = c(0.21, 0.70),
          legend.text = element_text(size = 13),
         panel.background = element_blank(),
         #panel.grid.major = element_line(gray(0.7)),
         panel.grid.major = element_blank(),
         panel.border = element_rect(fill = NA),
         axis.title.y = element_text(size = 30, colour = 'black'),
         axis.title.x = element_text(size = 30, colour = 'black'),
         axis.text.x = element_text(size = 20, colour = 'black'),
         axis.text.y = element_text(size = 20, colour = 'black'))
whittaker_plot


################ mean tmp and precip per site ####################
calculate_mean <- function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else first(x)
diagram_points <- leaf_analysis %>%
  group_by(site_code) %>%
  summarise_all(calculate_mean)
nrow(diagram_points) #26
names(diagram_points)
diagram_points$precip.cm = diagram_points$precip.mean/10
#################################################################################

####### table for longitude and latitude transformation into degree and minute ##########

# Function to convert decimal degrees to degrees and minutes
convert_to_deg_min <- function(decimal_deg, is_latitude = TRUE) {
  deg <- floor(abs(decimal_deg))
  minutes <- round((abs(decimal_deg) - deg) * 60, 2)
  
  direction <- ifelse(is_latitude, 
                      ifelse(decimal_deg < 0, "S", "N"),
                      ifelse(decimal_deg < 0, "W", "E"))
  
  return(paste0(deg, "° ", minutes, "' ", direction))
}

# Apply the function to the lat and long columns
diagram_points$lat_deg_min <- sapply(diagram_points$Lat, convert_to_deg_min, is_latitude = TRUE)
diagram_points$long_deg_min <- sapply(diagram_points$Long, convert_to_deg_min, is_latitude = FALSE)


#### Superimpose points on the diagram ############################
whittaker_ggplot_data <- whittaker_plot +
  geom_point(data = diagram_points, aes(x = tmp, y = precip.cm), color = "black", fill="white", pch=21, size = 3)
whittaker_ggplot_data


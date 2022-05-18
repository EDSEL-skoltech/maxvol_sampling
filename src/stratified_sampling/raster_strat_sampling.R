#install.packages("raster")
library(raster)


#setwd("../../source_data_2022/dem.tif")
img.dem <- raster("../../source_data_2022/dem.tif")

plot(img.dem)



smp.test <- sampleStratified(x = img.dem,
                             size = 3,
                             na.rm = TRUE,
                             sp = TRUE)

#plot(smp.test)

plot(img.dem, 
     axes = FALSE, 
     box = FALSE,
     #col = c("#fbf793", "#006601", "#bfe578", "#d00000", "#fa6700", "#6569ff")
     col = c("#fbf793", "#006601", "#bfe578")
)

points(smp.test)

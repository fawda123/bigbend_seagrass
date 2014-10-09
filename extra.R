source('.Rprofile')

##
# get data for doc

# segment polygon for old tampa bay
seg_shp <- readShapeSpatial('seagrass_gis/seg_820.shp')

# seagrass bathmetry intersected points for old tampa bay
sgpts_shp <- readShapeSpatial('seagrass_gis/sgpts_820_2006_buff.shp')

# set ggplot theme
theme_set(theme_bw())

# for debugging
grid_spc <- 0.02
grid_seed <- 1234
test_pt <- 2
radius <- 0.04
thresh <- 0.1   	
show_all <- F
    
# random points  
set.seed(grid_seed)
pts <- grid_est(seg_shp, spacing = grid_spc) 

# point from random points for buffer
test_pt <- pts[test_pt, ]

# get bathym points around test_pt
buff_pts <- buff_ext(sgpts_shp, test_pt, buff = radius)
  
##
# example

maxd <- list()
for(i in 1:length(pts)){
  
  eval_pt <- pts[i, ]
  buff_pts <- buff_ext(sgpts_shp, eval_pt, buff = radius)
	est_pts <- data.frame(buff_pts)
	est_pts$Depth <- -1 * est_pts$GRID_CODE
  ests <- doc_est(est_pts, thresh = thresh, 
  								depth_var = 'Depth', sg_var = 'SEAGRASS' 
  								)
  maxd[[i]] <- ests$ests
  
}

# combine in data frame for plotting
maxd <- data.frame(pts, zmax_all = do.call('c', maxd))

# get values for combined legend
rngs <- range(maxd$zmax_all, na.rm = T)
brks <- seq(rngs[1], rngs[2], length = 5)
labs <- format(round(brks, 1), nsmall = 1, digits =1)

# unestimable points to plot
unest <- maxd[is.na(maxd[, 'zmax_all']), ]

##
# plot

p1 <- ggplot(seg_shp, aes(long, lat)) + 
  geom_polygon(fill = 'white') +
  geom_path(color = 'black') +
  theme_classic() +
  coord_equal() +
	ylab('Latitude') + 
	xlab('Longitude') +
  geom_point(
    data = maxd, 
    aes(Var1, Var2, size = zmax_all, colour = zmax_all)
  ) +
# 				geom_point(data = unest,
# 					aes(Var1, Var2), 
# 					size = 3, colour = 'grey',
# 					pch = 1
# 					) +
  ggtitle('Depth of col (m)') +
	theme(legend.position = c(0,0), legend.justification = c(0, 0)) + 
	scale_size_continuous(name = "Depth estimate", 
		breaks = brks, 
		labels = labs,
		range = c(1, 12)) + 
	scale_colour_gradient(name = "Depth estimate", 
		breaks = brks, 
		labels = labs) + 
 	guides(colour = guide_legend())

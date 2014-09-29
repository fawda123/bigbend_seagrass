
# functions to use
source('funcs.R')

# segment polygon for old tampa bay
seg_shp <- readShapeSpatial('seagrass_gis/seg_820.shp')

# seagrass bathmetry intersected points for old tampa bay
sgpts_shp <- readShapeSpatial('seagrass_gis/sgpts_820_2006_buff.shp')

# set ggplot theme
theme_set(theme_bw())

# plotting code

# for debugging
grid_spc <- 0.02
grid_seed <- 1234
test_pt <- 19
radius <- 0.04
thresh <- 0.1
show_all <- F

# random points  
set.seed(grid_seed)
pts <- grid_est(seg_shp, spacing = grid_spc) 
 
######

# get estimates for each point
maxd <- list()
for(i in 1:length(pts)){
  
  eval_pt <- pts[i, ]
  buff_pts <- buff_ext(sgpts_shp, eval_pt, buff = radius)
  ests <- max_est(data.frame(buff_pts), thresh = thresh,  
  								depth_var = 'GRID_CODE', sg_var = 'SEAGRASS' 
  								)
  maxd[[i]] <- ests
  
}

maxd <- data.frame(do.call('rbind', maxd))
maxd <- data.frame(
	pts, 
	-1 * maxd
	)

rngs <- range(maxd$zmax_all, na.rm = T)
brks <- seq(rngs[1], rngs[2], length = 5)
labs <- format(round(brks, 1), nsmall = 1, digits =1)

unest <- maxd[is.na(maxd[, 'zmax_all']), ]

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
	geom_point(data = unest,
		aes(Var1, Var2), 
		size = 3, 
		pch = 1
		) +
	theme(legend.position = c(0,0), legend.justification = c(0, 0)) + 
  ggtitle('Depth of col (m)') +
	scale_size_continuous(name = "Depth estimate", 
		breaks = brks, 
		labels = labs,
		range = c(1, 12)) + 
	scale_colour_gradient(name = "Depth estimate", 
		breaks = brks, 
		labels = labs) + 
 	guides(colour = guide_legend())

print(p1)

######
    	
# point from random points for buffer
test_pt <- pts[test_pt, ]

# get bathym points around test_pt
buff_pts <- buff_ext(sgpts_shp, test_pt, buff = radius)

p1 <- ggplot(seg_shp, aes(long, lat)) + 
  geom_polygon(fill = 'white') +
  geom_path(color = 'black') +
  theme_classic() +
  coord_equal() +
	xlab('Longitude') + 
	ylab('Latitude') +
  geom_point(
    data = data.frame(pts), 
    aes(Var1, Var2), 
    size = 3,
    pch = 1
  ) +
  geom_point(
    data = data.frame(buff_pts),
    aes(coords.x1, coords.x2), 
    colour = 'red', 
    size = 0.3, 
    alpha = 0.7
  ) +
	geom_point(
		data = data.frame(test_pt), 
		aes(Var1, Var2), 
		size = 3, 
		pch = 1
	)

# get data used to estimate depth of col
ests <- max_est(data.frame(buff_pts), thresh = thresh, 
  								depth_var = 'GRID_CODE', sg_var = 'SEAGRASS', 
									dat_out = T)

# actual ests
act_ests <- max_est(data.frame(buff_pts), thresh = thresh, 
  								depth_var = 'GRID_CODE', sg_var = 'SEAGRASS'
  								)
if(is.na(act_ests)){ act_ests <- 'Depth of col: Not estimable'
} else { 
	act_ests <- paste('Depth of col:', -1 * round(act_ests, 1), 'm')
	}

##
# simple plot of points by depth, all pts and those with seagrass
to_plo <- ests
to_plo <- melt(to_plo, id.var = 'Depth', 
	measure.var = c('dep_cum', 'sg_cum'))
to_plo$variable <- factor(to_plo$variable, levels = c('dep_cum', 'sg_cum'), 
                            labels = c('All', 'Seagrass'))

cols  <- c('lightgreen', 'lightblue')
linesz <- 1

p2 <- ggplot(to_plo, aes(x = -1 * Depth, y = value, group = variable,
                         colour = variable)) +
  geom_line(size = linesz) +
 	ylab('No. of points') +
  xlab('Depth (m)') +
  scale_colour_manual('Point category', values = cols) +
  theme(legend.position = c(0, 1), legend.justification = c(0,1))

##
# plot slope of cumulative point curves

# treshold label for legend
thresh_lab <- paste0(round(100 * thresh), '% of all')

to_plo <- ests
to_plo <- melt(to_plo, id.var = 'Depth', 
	measure.var = c('dep_slo', 'sg_slo', 'slo_thr'))
to_plo$variable <- factor(to_plo$variable, 
	levels = c('dep_slo', 'sg_slo', 'slo_thr'), 
  labels = c('All', 'Seagrass', thresh_lab))

p3 <- ggplot(to_plo, aes(x = -1 * Depth, y = value)) +
			  geom_line(size = linesz, 
			  	aes(group = variable, colour = variable, linetype = variable)
			  	) +
			 	ylab('CDF Slope') +
			  xlab('Depth (m)') +
			  scale_colour_manual('Slope category', values = c(cols, cols[1])) +
				scale_linetype_manual('Slope category', values = c('solid', 'solid', 'dashed')) + 
			  theme(legend.position = c(1, 1), legend.justification = c(1, 1))

grid.arrange(p1,
	arrangeGrob(p2, p3, ncol = 2), 
	ncol = 1, heights = c(1.25, 1),
	main = act_ests)


   
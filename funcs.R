######
# functions for seagrass depth of col estimates

######
# function for estimating max depth of col
# 'dat_in' is data frame of binned proportions, returned from 'sg_bins' func
# 'depth_var' is name of depth column in input data
# 'sg_var' is name of seagrass column in input data
# 'cont_name' variable in 'sg_var' that defines continuous growth
# 'dat_out' is logical indicating if cumulative dist ests are returned
max_est <- function(dat_in, depth_var = 'depth', sg_var = 'seagrass', 
	dat_out = F){
  
	# order by depth, assumes column is negative
  dat_in <- dat_in[order(dat_in[, depth_var], decreasing = T), ]
	dat_in$depth <- dat_in[, depth_var]
	
	# cumulative sum of pts with all seagrass and all points
	# assumes NA is empty
	sg_pts <- rev(table(dat_in[!is.na(dat_in[, sg_var]), depth_var]))
	sg_pts <- data.frame(Depth = names(sg_pts), sg_pts = sg_pts,
		sg_cum = cumsum(sg_pts), row.names = 1:length(sg_pts))
	dep_pts <- rev(table(dat_in[, depth_var]))
	dep_pts <- data.frame(Depth = names(dep_pts), dep_pts = dep_pts, 
		dep_cum = cumsum(dep_pts), row.names = 1:length(dep_pts))
	
	# combine all pts and seagrass pts, depth as numeric
	pts <- merge(dep_pts, sg_pts, by = 'Depth', all.x = T)
	pts$Depth <- as.numeric(as.character(pts$Depth))
	pts$sg_prp <- with(pts, sg_pts/dep_pts)
	
	# add slope ests to pts
	pts$dep_slo <- with(pts, c(NA, -1 * diff(dep_cum)/diff(Depth)))
	pts$sg_slo <- with(pts, c(NA, -1 * diff(sg_cum)/diff(Depth)))
	
	# return cumulative data if T
	if(dat_out) return(pts)
	
	# get max (95th) and median (50th) depth for all seagrass
	max_depth_all <- dat_in[which.min(abs(dat_in[, 'all'] - 0.95)), 'depth']
	med_depth_all <- dat_in[which.min(abs(dat_in[, 'all'] - 0.5)), 'depth']
	
	# get max (95th) and median (50th) depth for continuous seagrass
	max_depth_cont <- dat_in[which.min(abs(dat_in[, 'cont'] - 0.95)), 'depth']
	med_depth_cont <- dat_in[which.min(abs(dat_in[, 'cont'] - 0.5)), 'depth']
		
  # return output
  out <- c(zmax_all = max_depth_all, z50_all = med_depth_all, 
  	zmax_cont = max_depth_cont, z50_cont = med_depth_cont)
  return(out)
  
}

#######
# function for creating random grid of points, bounded by polygon extent
# taken from ibi sampling manuscript functions
# 'clip_poly' is shapefile input object
# 'spacing' is spacing between points, as degrees
grid_est <- function(clip_poly, spacing = 0.03){
  
  if(!'SpatialPolygonsDataFrame' %in% class(clip_poly))
    stop('clip_poly must be of class SpatialPolygonsDataFrame')
  
  library(sp) 
  
  # extent of shapefile
  extent <- summary(clip_poly)$bbox
  
  # buffer of shapefile and random starting x/y locs
  add.on <- apply(extent, 1, diff) * 0.3
  rand <- runif(2, 0, spacing)
  
  # random points within rectangular extent
  pts<-{
    x.vals<-seq(extent[1, 1] - add.on['x'], extent[1, 2] + add.on['x'], by = spacing) + rand[1]
    y.vals<-seq(extent[2, 1] - add.on['y'], extent[2, 2] + add.on['y'], by = spacing) + rand[2]
    expand.grid(x.vals, y.vals)
  }
  
  # clip by clip_poly and return
  sel <- !is.na(SpatialPoints(pts) %over% clip_poly)[, 1]
  
  return(SpatialPoints(pts)[sel, ])
  
}

######
# function extracts bathymetric seagrass pts withing a distance from a pt
# 'pts' is spatial points to extract
# 'center' is pt from which buffer extends
# 'buff' is radius of buffer in dec degrees
buff_ext <- function(pts, center, buff = 0.03){
  
  # sanity checks
  if(!any(c('SpatialPointsDataFrame', 'SpatialPoints') %in% class(pts)))
    stop('pts must be of class SpatialPointsDataFrame or SpatialPoints')
  
  if(!any(c('SpatialPointsDataFrame', 'SpatialPoints') %in% class(center)))
    stop('center must be of class SpatialPointsDataFrame or SpatialPoints')
  
  library(rgeos)
  library(sp)
  
  # create buffer
  buffer <- gBuffer(center, width = buff)
  
  # index of pts in buffer
  sel <- !is.na(pts %over% buffer)
  
  if(sum(sel) == 0) stop('No points in buffer')
  
  return(pts[sel, ])
  
}
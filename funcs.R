######
# functions for seagrass depth of col estimates

######
# function for estimating depth of colonization
# also used for plots
# 'dat_in' is data from 'buff_ext'
# 'depth_var' is name of depth column in input data
# 'sg_var' is name of seagrass column in input data
# 'thresh' is numeric threshold value for estimating depth of col
doc_est <- function(dat_in, depth_var = 'depth', sg_var = 'seagrass',
	thresh = 0.1){
  
	# order by depth, assumes column is negative
  dat_in <- dat_in[order(dat_in[, depth_var], decreasing = T), ]
	dat_in$depth <- dat_in[, depth_var]
	
	# cumulative sum of pts with all seagrass and all points
	# assumes NA is empty
	sg_pts <- table(dat_in[!is.na(dat_in[, sg_var]), depth_var])
	sg_pts <- data.frame(Depth = names(sg_pts), sg_pts = as.numeric(sg_pts),
		sg_cum = cumsum(sg_pts), row.names = 1:length(sg_pts))
	dep_pts <- table(dat_in[, depth_var])
	dep_pts <- data.frame(Depth = names(dep_pts), dep_pts = as.numeric(dep_pts), 
		dep_cum = cumsum(dep_pts), row.names = 1:length(dep_pts))
	
	# combine all pts and seagrass pts, depth as numeric
	pts <- merge(dep_pts, sg_pts, by = 'Depth', all.x = T)
	pts$Depth <- as.numeric(as.character(pts$Depth))
	pts$sg_prp <- with(pts, sg_pts/dep_pts)
	
	# add slope ests to pts, use differences
	pts$dep_slo <- with(pts, c(NA, diff(dep_cum)/diff(Depth)))
	pts$sg_slo <- with(pts, c(NA, diff(sg_cum)/diff(Depth)))
	
	pred_ls <- vector('list', length = 2)
	names(pred_ls) <- c('sg_slo', 'dep_slo')
	
	sg_slo_max <-pts[which.max(pts[, 'sg_slo']), 'Depth']
	dep_slo_max <- pts[which.max(pts[,'dep_slo']), 'Depth']
	min_dep <- min(c(sg_slo_max, dep_slo_max))	
	
	for(var in c('sg_slo', 'dep_slo')){
		
		# scale y values for easier parameter estimates
		scl_val <- max(pts[, var], na.rm = T)

		# subset data by maximum and minimum slope value
		max_ind <- which.max(pts[, var])
		pts_sub <- pts[which.max(pts[, var]):nrow(pts), ]
		pts_sub[, var] <- pts_sub[, var]/scl_val
		min_ind <- which(pts_sub[, var] < 0.05)[1]
		if(!is.na(min_ind))
			pts_sub <- pts_sub[1:min_ind,]
	
		# function for calculating negative log likelihood
		err_est <- function(parms = list(b1, b2, b3)){
			
			b1 <- parms[[1]]
			b2 <- parms[[2]]
			b3 <- parms[[3]]
			
			to_comp <- b1 + b2 * pts_sub$Depth
			act <-  pts_sub[, var]
			
			resid <- na.omit(to_comp - act)
			resid <- suppressWarnings(dnorm(resid, 0, b3))
			
			-sum(log(resid))
			
			}

		# find estimates using optim and get predictions
		bs <- list(b1 = 1, b2 = -1, b3 = 1)
		res <- optim(par = bs, fn = err_est)$par
		names(res) <- c('b1', 'b2', 'b3')
		
		# get predictions
		max_dep <- -res[['b1']]/res[['b2']]
		new.x <- seq(min_dep, max_dep, by = 0.01)
		pred <- res['b1'] + res[['b2']] * new.x
		pred <- data.frame(Depth = new.x, scl_val * pred)
		names(pred) <- c('Depth', var)

		pred_ls[[var]] <- pred
		
	}
	
	# return output
	out <- merge(pred_ls[[1]], pred_ls[[2]], by = 'Depth', all = T)
	
	# add threshold data based on proportion of dep_slo 
	threshs <- sapply(1:length(thresh), 
		FUN = function(x) thresh[x] * out[, 'dep_slo']
		)
	threshs <- data.frame(threshs)
	names(threshs) <- paste('Threshold', thresh)
	
	out <- data.frame(out, threshs)
	
	# calculate depth of col
	doc <- sapply(thresh, 
		FUN = function(x){
			
			col <- out[, grep(x, names(out))]
			ind <- which(with(out,  sg_slo <= col))[1]
			out[ind, 'Depth']
			
			}
		)
	names(doc) <- thresh

	return(list(data = pts, thresh = out, ests = doc))
	  
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

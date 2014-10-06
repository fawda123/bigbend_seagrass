
# segment polygon for old tampa bay
seg_shp <- readShapeSpatial('seagrass_gis/seg_820.shp')

# seagrass bathmetry intersected points for old tampa bay
sgpts_shp <- readShapeSpatial('seagrass_gis/sgpts_820_2006_buff.shp')

# set ggplot theme
theme_set(theme_bw())

# for debugging
grid_spc <- 0.02
grid_seed <- 12
test_pt <- 25
radius <- 0.02
thresh <- 0.1   	
show_all <- F
    

# random points  
set.seed(grid_seed)
pts <- grid_est(seg_shp, spacing = grid_spc) 

# point from random points for buffer
test_pt <- pts[test_pt, ]

# get bathym points around test_pt
buff_pts <- buff_ext(sgpts_shp, test_pt, buff = radius)
  
# get data used to estimate depth of col
ests <- max_est(data.frame(buff_pts), thresh = thresh, 
  								depth_var = 'GRID_CODE', sg_var = 'SEAGRASS', 
									dat_out = T)
    	
ests$Depth <- -1 * ests$Depth

par(mfrow = c(1,2), family = 'serif')
plot(dep_cum ~ Depth, ests, col = 'lightblue', type = 'l')
lines(sg_cum ~ Depth,  ests, col = 'lightgreen')
plot(dep_slo ~ Depth, ests, col = 'lightblue', pch = 16)
points(sg_slo ~ Depth, ests, col = 'lightgreen', pch = 16)

##
# parameter optimization for a linear regression

set.seed(123)
vals <- 1000
x <- runif(vals, 0, 15)
y <- 20 + 5 * x + rnorm(vals, 0,10)
plot(x, y)

dat_in <- data.frame(x, y)

err_est <- function(parms = list(b1, b2), dat = dat_in, 
	var_names = c('x', 'y')){
	
	b1 <- parms[[1]]
	b2 <- parms[[2]]
	
	pred <- b1 + b2 * dat[, var_names[1]]
	act <- dat[, var_names[2]]
	
	sum((pred - act)^2)
	
	}

optim(c(0, 0), err_est)

##
# find parameters using least squares

set.seed(123)
vals <- 100
x <- runif(vals, 0, 15)
y <- dgamma(x, shape = 2, scale = 2) + rnorm(vals, 0, 0.1) # errors are norm and iid
plot(x, y)
dat_in <- data.frame(x, y)

err_est <- function(parms = list(b1, b2), dat = dat_in, 
	var_names = c('x', 'y')){
	
	b1 <- parms[[1]]
	b2 <- parms[[2]]
	
	to_comp <- dgamma(dat[, var_names[1]], shape = b1, scale = b2)
	act <- dat[, var_names[2]]
	
	# sums of squared errors
	sum((to_comp - act)^2, na.rm = T)
	
	}


res <- optim(c(2,2), err_est)$par
names(res) <- c('b1', 'b2')
new.x <- seq(0, 15, length = 500)
pred <- dgamma(new.x, shape = res['b1'], scale = res['b2'])

plot(x, y)
lines(new.x, pred, col = 'blue')

##
# find parameters using negative log-likelihood
set.seed(123)
vals <- 100
x <- runif(vals, 0, 15)
y <- dgamma(x, shape = 2, scale = 2) + rnorm(vals, 0, 0.1) # errors are norm and iid
plot(x, y)
dat_in <- data.frame(x, y)

err_est <- function(parms = list(b1, b2, b3), dat = dat_in, 
	var_names = c('x', 'y')){
	
	b1 <- parms[[1]]
	b2 <- parms[[2]]
	b3 <- parms[[3]]
	
	# estimated and actual data
	to_comp <- dgamma(dat[, var_names[1]], shape = b1, scale = b2)
	act <- dat[, var_names[2]]
	
	# residuals, then random ests from dnorm
	resids <- act - to_comp
	resids <- suppressWarnings(dnorm(resids, 0, b3))
	
	# neg log-likelihood
	-sum(log(resids))
	
	}

res <- optim(c(2, 2, 0.1), err_est)$par
names(res) <- c('b1', 'b2')
new.x <- seq(0, 15, length = 500)
pred <- dgamma(new.x, shape = res['b1'], scale = res['b2'])

plot(x, y)
lines(new.x, pred, col = 'blue')

##
# example w/ actual data

# get data used to estimate depth of col
thresh <- c(0.1, 0.5)

est_pts <- data.frame(buff_pts)
est_pts$Depth <- -1 * est_pts$GRID_CODE

# data
dat <- doc_est(est_pts, thresh = thresh, 
	depth_var = 'Depth', sg_var = 'SEAGRASS', 
	dat_out = T
	)

# actual ests
act_ests <- doc_est(est_pts, thresh = thresh,
	depth_var = 'Depth', sg_var = 'SEAGRASS'
	)


# format estimate for plot title
if(is.na(act_ests)){ act_ests <- 'Depth of col: Not estimable'
} else { 
	act_ests <- paste('Depth of col:', round(act_ests, 1), 'm')
	}

##
# simple plot of points by depth, all pts and those with seagrass
to_plo <- dat
to_plo <- melt(to_plo, id.var = 'Depth', 
	measure.var = c('dep_cum', 'sg_cum'))
to_plo$variable <- factor(to_plo$variable, levels = c('dep_cum', 'sg_cum'), 
                            labels = c('All', 'Seagrass'))

cols  <- c('lightgreen', 'lightblue')
linesz <- 1

p2 <- ggplot(to_plo, aes(x = Depth, y = value, group = variable,
                         colour = variable)) +
  geom_line(size = linesz) +
 	ylab('Cumulative points') +
  xlab('Depth (m)') +
  scale_colour_manual('Point category', values = cols) +
  theme(legend.position = c(0, 1), legend.justification = c(0,1))

##
# plot slope of cumulative point curves

# treshold label for legend
thresh_lab <- paste0(round(100 * thresh), '% of all')

to_plo <- dat
to_plo <- melt(to_plo, id.var = 'Depth', 
	measure.var = c('dep_slo', 'sg_slo'))
to_plo$variable <- factor(to_plo$variable, 
	levels = c('dep_slo', 'sg_slo'), 
  labels = c('All', 'Seagrass'))

to_plo2 <- dat
to_plo2 <- melt(to_plo2, id.var = 'Depth', 
	measure.var = grep('dep_est|sg_est|Threshold', names(dat), value = T)
	)
to_plo2$variable <- factor(to_plo2$variable, 
	labels = c('All', 'Seagrass', thresh_lab)
)

col_pts <- rep(cols[1], nrow(to_plo))
col_pts[to_plo$variable == 'All'] <- cols[2]

p3 <- ggplot(to_plo, aes(x = Depth, y = value)) +
  geom_point(size = 3, shape = 16, colour = col_pts
  	) +
	geom_line(data = to_plo2, aes(x = Depth, y = value,
		colour = variable, linetype = variable), size = linesz) +
 	ylab('CDF Slope') +
  xlab('Depth (m)') +
  scale_colour_manual('Slope category', 
  	values = c(cols[2], cols[1], cols[2], cols[2])
  	) +
	scale_linetype_manual('Slope category', 
		values = c('solid', 'solid', 'dashed', 'dashed')
		) + 
  theme(legend.position = c(1, 1), legend.justification = c(1, 1))

##
# combine all plots

grid.arrange(p1,
	arrangeGrob(p2, p3, ncol = 2), 
	ncol = 1, heights = c(1.25, 1),
	main = act_ests)

# segment polygon for old tampa bay
seg_shp <- readShapeSpatial('seagrass_gis/seg_820.shp')

# seagrass bathmetry intersected points for old tampa bay
sgpts_shp <- readShapeSpatial('seagrass_gis/sgpts_820_2006_buff.shp')

# set ggplot theme
theme_set(theme_bw())

# for debugging
grid_spc <- 0.02
grid_seed <- 12
test_pt <- 1
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
plot(sg_slo ~ Depth, ests)

mod_dep <- glm(dep_slo ~ Depth, family = Gamma(link = 'log'), ests,
	na.action = na.exclude)

gshape_dep <- gamma.shape(mod_dep)
pred_dep <- predict(mod_dep, type = "response", se = T, 
	dispersion = 1/gshape_dep$alpha)

mod_sg <- glm(sg_slo ~ Depth, family = Gamma(link = 'log'), ests,
	na.action = na.exclude)

gshape_sg <- gamma.shape(mod_sg)
pred_sg <- predict(mod_sg, type = "response", se = T, 
	dispersion = 1/gshape_sg$alpha)

ylims <- c(0, max(ests$dep_slo, na.rm = T) * 1.1)
plot(dep_slo ~ Depth, ests, ylim = ylims, col = 'blue')
points(sg_slo ~ Depth, ests, col = 'darkgreen')
lines(ests$Depth, pred_dep$fit, col = 'blue')
lines(ests$Depth, pred_sg$fit, col = 'darkgreen')




# for a linear regression

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
	
	to_comp <- b1 + b2 * dat[, var_names[1]]
	act <- dat[, var_names[2]]
	
	sum((to_comp - act)^2)
	
	}

optim(c(0, 0), err_est)

##



set.seed(123)
vals <- 100
x <- runif(vals, 0, 15)
y <- 20 + 5 * dgamma(x, shape = 2, scale = 2) + rnorm(vals, 0, 0.1)
plot(x, y)
dat_in <- data.frame(x, y)

err_est <- function(parms = list(b1, b2, b3, b4), dat = dat_in, 
	var_names = c('x', 'y')){
	
	b1 <- parms[[1]]
	b2 <- parms[[2]]
	b3 <- parms[[3]]
	b4 <- parms[[4]]
	
	to_comp <- b1 + b2 * dgamma(dat[, var_names[1]], shape = b3, scale = b4)
	act <- dat[, var_names[2]]
	
	sum((to_comp - act)^2, na.rm = T)
	
	}


res <- optim(c(0,2,2,2), err_est)$par
names(res) <- c('b1', 'b2', 'b3', 'b4')
new.x <- seq(0, 15, length = 500)
pred <- res['b1'] + res['b2'] * dgamma(new.x, shape = res['b3'], scale = res['b4'])

plot(x, y)
lines(new.x, pred, col = 'blue')

##

# get data used to estimate depth of col
ests <- max_est(data.frame(buff_pts), thresh = thresh, 
  								depth_var = 'GRID_CODE', sg_var = 'SEAGRASS', 
									dat_out = T)
    	
ests$Depth <- -1 * ests$Depth

err_est <- function(parms = list(b1, b2, b3, b4), dat = ests, 
	var_names = c('Depth', 'sg_slo')){
	
	b1 <- parms[[1]]
	b2 <- parms[[2]]
	b3 <- parms[[3]]
	b4 <- parms[[4]]
	
	to_comp <- b1 + b2 * dgamma(dat[, var_names[1]], shape = b3, scale = b4)
	act <- dat[, var_names[2]]
	
	sum((to_comp - act)^2, na.rm = T)
	
	}


res <- nlminb(c(5,5,2,0.2), err_est)$par
names(res) <- c('b1', 'b2', 'b3', 'b4')
new.x <- seq(0, 15, length = 500)
pred <- res['b1'] + res['b2'] * dgamma(new.x, shape = res['b3'], scale = res['b4'])

plot(sg_slo ~ Depth, ests)
lines(new.x, pred, col = 'blue')

##
# gamma fun, manual

gamma_fun <- function(x, alpha, beta){
	
	numer <- (beta^alpha) * (x^(alpha - 1)) * exp(-x * beta)
	denom <- gamma(alpha)
	
	out <- numer/denom
	
	return(out)
	
	}

x <- sort(runif(1000, 0, 15))
y <- gamma_fun(x, 2, 0.2)

plot(x, y,  type = 'l')

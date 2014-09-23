# packages to use
library(maptools)
library(reshape2) 
library(plyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(rgdal)
library(rgeos)
library(gridExtra)

# functions to use
source('funcs.R')

# segment polygon for old tampa bay
seg_shp <- readShapeSpatial('seagrass_gis/seg_820.shp')

# seagrass bathmetry intersected points for old tampa bay
sgpts_shp <- readShapeSpatial('seagrass_gis/sgpts_820_2006.shp')

# set ggplot theme
theme_set(theme_bw())

# Define server logic required to generate and plot data
shinyServer(function(input, output) {
  
  # dynamic controls
  # pick test pt once pts are selected
  output$reserveControls <- renderUI({

    grid_spc <- input$grid_spc
    grid_seed <- input$grid_seed
    set.seed(grid_seed)
    pts <- grid_est(seg_shp, spacing = grid_spc) 
    
    selectInput(inputId = 'test_pt',
                label = h3('Test point'),
                choices = 1:length(pts)
      )
    
    })
  
  output$simplot <- renderPlot({
    
    # plotting code

#     # for debugging
#     grid_spc <- 0.02
#     grid_seed <- 12
#     test_pt <- 1
#     radius <- 0.02
#     show_all <- F
#     est <- 'all - zmax'
    
    # input from ui
    grid_spc <- input$grid_spc
    grid_seed <- input$grid_seed
    test_pt <- input$test_pt
    radius <- input$radius
    show_all <- input$show_all
    est <- input$est
    
    # random points  
    set.seed(grid_seed)
    pts <- grid_est(seg_shp, spacing = grid_spc) 
    
    if(show_all){
      
      est <- strsplit(est, ' - ')[[1]]
      
    	# get estimates for each point
      maxd <- list()
      for(i in 1:length(pts)){
        
        eval_pt <- pts[i, ]
        buff_pts <- buff_ext(sgpts_shp, eval_pt, buff = radius)
        ests <- max_est(data.frame(buff_pts), 
        								depth_var = 'GRID_CODE', sg_var = 'SEAGRASS', 
        								cont_name = 'Continuous')
      	if(length(ests) == 0) ests <- rep(NA_real_, 4)
        maxd[[i]] <- ests
        
      }
      
    	maxd <- data.frame(do.call('rbind', maxd))
      maxd <- data.frame(
      	pts, 
      	ptsz = maxd[, grepl(paste0(est[2], '_', est[1]), names(maxd))]
      	)
      
      p1 <- ggplot(seg_shp, aes(long, lat)) + 
        geom_polygon(fill = 'white') +
        geom_path(color = 'black') +
        theme_classic() +
        coord_equal() +
        geom_point(
          data = maxd, 
          aes(Var1, Var2, size = ptsz, colour = ptsz)
        ) +
        scale_size(range = c(1, 12)) +
        theme(legend.title = element_blank()) + 
        ggtitle(paste(est[1], '-', est[2], '(m)'))
      
      print(p1)
      
    } else {
    
      # point from random points for buffer
      test_pt <- pts[test_pt, ]
  
      # get bathym points around test_pt
      buff_pts <- buff_ext(sgpts_shp, test_pt, buff = radius)
      
      p1 <- ggplot(seg_shp, aes(long, lat)) + 
        geom_polygon(fill = 'white') +
        geom_path(color = 'black') +
        theme_classic() +
        coord_equal() +
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
        )
      
      # get data used to estimate depth of col
			ests <- max_est(data.frame(buff_pts), 
        								depth_var = 'GRID_CODE', sg_var = 'SEAGRASS', 
        								cont_name = 'Continuous',
												dat_out = T)
    	
      # simple plot of proportion of 'occupied' points at depth bins
      to_plo <- ests
      to_plo <- melt(to_plo, id.var = 'depth')
      to_plo$growth_cat <- factor(to_plo$variable, levels = c('all', 'cont'), 
                                  labels = c('All', 'Continuous'))
      
      # actual ests
    	act_ests <- max_est(data.frame(buff_pts), 
        								depth_var = 'GRID_CODE', sg_var = 'SEAGRASS', 
        								cont_name = 'Continuous')
    	act_ests <- round(act_ests, 2)
      
      ##
      # plot cumulative distribution curves from data
      
      cols  <- c('lightgreen', 'lightblue')
      linesz <- 1
      
      p2 <- ggplot(to_plo, aes(x = value, y = depth, group = variable,
                               colour = growth_cat)) +
        geom_line(size = linesz) +
       	xlab('Cumulative frequency') +
        ylab('Depth (m)') +
        scale_colour_manual('Growth\ncategory', values = cols) +
        theme(legend.position = c(0.8, 0.8))
      
      ests <- paste0('All: Zmax ', act_ests[1], ' Z50 ', act_ests[2], '\n', 
                     'Continuous: Zmax ', act_ests[3], ' Z50 ', act_ests[4])
      
      grid.arrange(p1, p2, 
                   ncol = 1, 
                   main = ests
      )
      
    }
    
    },height = 700, width = 600)

    })
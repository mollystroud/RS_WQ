################################################################################
# Code written by Molly Stroud on 11/3/25
################################################################################

################################################################################
# the below code creates the bounding boxes for the Roanoke reservoirs
################################################################################
# beaverdam reservoir bbox
bvr_box <- c(xmin = -79.827088, 
             ymin = 37.311798, 
             xmax = -79.811865, 
             ymax = 37.321694)
# falling creek reservoir bbox
fcr_box <- c(xmin = -79.840037, 
             ymin = 37.301435, 
             xmax = -79.833651, 
             ymax = 37.311487)
# carvins cove reservoir bbox
ccr_box <- c(xmin = -79.981728, 
             ymin = 37.367522, 
             xmax = -79.942552, 
             ymax = 37.407255)

# convert the bounding boxes to the correct UTM projection
bvr_box_utm <- sf::st_bbox(
  sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(xmin = -79.827088, 
                                               ymin = 37.311798, 
                                               xmax = -79.811865, 
                                               ymax = 37.321694), 
                                             crs = "EPSG:4326")), "EPSG:32617")
)
fcr_box_utm <- sf::st_bbox(
  sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(xmin = -79.840037, 
                                               ymin = 37.301435, 
                                               xmax = -79.833651, 
                                               ymax = 37.311487), 
                                             crs = "EPSG:4326")), "EPSG:32617")
)
ccr_box_utm <- sf::st_bbox(
  sf::st_transform(sf::st_as_sfc(sf::st_bbox(c(xmin = -79.981728, 
                                               ymin = 37.367522, 
                                               xmax = -79.942552, 
                                               ymax = 37.407255), 
                                             crs = "EPSG:4326")), "EPSG:32617")
)



#
# CVS - $Id: lib_data_processing.R,v 1.1 2010/01/19 15:41:31 lmilazzo Exp $
#



#
# Library for Data Processing [vers. 0.7.0]
#


#
# Let us consider a set of ordered data and a grid.
#
# The data set contains N data structures. Each data structure
# contains N' elements - the first two elements are the X and Y coordinates,
# the last element is the data structure ID.
#
# The grid can be square or rectangular. The grid point at the bottom left
# corner is the first grid point (ID = 1) and it is associated with the
# spatial location (1,1) - in general, a `coordinate transformation' is needed.
#
# Let us associate each data structure with a grid point. The (ordered) data
# structures are associated with the grid points row by row, starting from the
# bottom left corner and ending to the top right corner of the grid.
#

#
# data structure format:
#        [x, y, p1, p2, p3, id]
# where
#        x  = X coord.
#        y  = Y coord.
#        p1 = parameter 1
#        p2 = parameter 2
#        p3 = parameter 3
#        id = data structure ID
#



#--

#
# Given a data structure associated with a spatial location, the function
# `getGridPointID()' evaluates its grid location (ID of the grid point).
# Note that the first two elements of the data structure are the X and Y coord.
#
# input:
#   x            = X coord. (in original units) 
#   y            = Y coord. (in original units)
#   grid_specs   = grid specs: grid resolution (nc x nr) and grid size step
#
#
getGridPointID = function(x, y, grid_specs){

  offset = grid_specs[1]

  grid_point_id = round((x-1)/grid_specs[3] + 1) + round((y-1)/grid_specs[3]) * offset

  return(grid_point_id)

}



#
# Given a data structure associated with a grid location, the function
# `getCoords()' evaluates its spatial location (X and Y coord.).
# Note that the last element of the data structure is the ID of the grid point.
#
# input:
#   grid_point_id  = grid point ID
#   grid_specs     = grid specs: grid resolution (nc x nr) and grid size step
#
#
getCoords = function(grid_point_id, grid_specs){

  offset = grid_specs[1]

  # distance of the grid point from the left edge of the grid
  dist_from_le = grid_point_id - 1
  while(dist_from_le>=offset){
    dist_from_le = dist_from_le - offset 
  }
  x = (dist_from_le * grid_specs[3]) + 1

  # distance of the grid point from the bottom edge of the grid
  dist_from_be = 0
  scan_dist_be = grid_point_id
  while(scan_dist_be>offset){
    scan_dist_be = scan_dist_be - offset
    dist_from_be = dist_from_be + 1
  }
  y = (dist_from_be * grid_specs[3]) + 1

  return(c(x, y))

}



#--


#
# The function `getROIPointsID()' evaluates the grid locations (IDs
# of the grid points) within a given ROI.
# 
# input:
#   no_roi_points  = no. of grid points within the ROI
#   offset         = no. of grid points for each row (grid width)
#   roix           = no. of grid points within the ROI in the X direction (ROI width)
#   lbound         = first grid location at the bottom left corner of the ROI
#
#
getROIPointsID = function(no_roi_points, offset, roix, lbound){

  # define data structure for the IDs of the
  # points within the ROI
  roi_points_id = rep(0,no_roi_points)
  roi_i = lbound
  j = 0

  for(i in 1:no_roi_points){
    roi_points_id[i] = roi_i
    roi_i = roi_i + 1
    j = j + 1
    if(j == roix){
      roi_i = roi_i + offset - roix
      j = 0
    }
  }

  return(roi_points_id)

}



#
# ROI_SFD = the Region of Interest is a square and it is localized within
#           a finite domain
#
# Let us consider the neighbourhood (Region of Interest, ROI) of a grid point
# (i.e. data structure). The ROI is a square centred at the grid point and
# containing the neighbour points.
# The radius of the ROI is `roi_r'. The area of the ROI depends on the radius
# and the location within the domain.
#

#
# Given a set of the data structures associated with a ROI, the function
# `evaluateROI_SFD()' evaluates the grid locations (IDs of the grid points).
# 
# input:
#   grid_res      = grid resolution (nc x nr)
#   ref_point_id  = ID of the reference point (center of the ROI)
#   roi_r         = radius of the ROI (in grid units)
#
#
evaluateROI_SFD = function(grid_res, ref_point_id, roi_r){

  # grid specs
  no_grid_points = grid_res[1] * grid_res[2]
  offset = grid_res[1]
  
  # distance of the reference point from the left edge of the grid
  dist_from_le = ref_point_id - 1
  while(dist_from_le>=offset){
    dist_from_le = dist_from_le - offset 
  }

  # distance of the reference point from the right edge of the grid
  dist_from_re = ref_point_id
  while(dist_from_re>offset){
    dist_from_re = dist_from_re - offset 
  }
  dist_from_re = offset - dist_from_re

  # distance of the reference point from the bottom edge of the grid
  dist_from_be = 0
  scan_dist_be = ref_point_id
  while(scan_dist_be>offset){
    scan_dist_be = scan_dist_be - offset
    dist_from_be = dist_from_be + 1
  }

  # distance of the reference point from the top edge of the grid
  dist_from_te = 0
  scan_dist_te = ref_point_id
  while(scan_dist_te<=no_grid_points){
    scan_dist_te = scan_dist_te + offset
    dist_from_te = dist_from_te + 1
  }
  dist_from_te = dist_from_te - 1


  # input validation
  if(roi_r>dist_from_le & roi_r>dist_from_re &
     roi_r>dist_from_be & roi_r>dist_from_te){ 
    cat("== error == unacceptable value for roi_r (the ROI cannot be bigger than the grid)\n")
    stop("script stopped due to an error")
  }


  #--
  # bottom left corner

  if(roi_r>dist_from_le & roi_r>dist_from_be){ 
    roix = dist_from_le + roi_r + 1
    roiy = dist_from_be + roi_r + 1
    no_roi_points = roix*roiy
    lbound = 1
  }

  #--
  # bottom right corner

  if(roi_r>dist_from_re & roi_r>dist_from_be){
    roix = dist_from_re + roi_r + 1
    roiy = dist_from_be + roi_r + 1
    no_roi_points = roix*roiy
    lbound = scan_dist_be - roi_r
  }

  #--
  # top left corner

  if(roi_r>dist_from_le & roi_r>dist_from_te){
    roix = dist_from_le + roi_r + 1
    roiy = dist_from_te + roi_r + 1
    no_roi_points = roix*roiy
    lbound = ref_point_id - roi_r*offset - dist_from_le
  }

  #--
  # top right corner

  if(roi_r>dist_from_re & roi_r>dist_from_te){
    roix = dist_from_re + roi_r + 1
    roiy = dist_from_te + roi_r + 1
    no_roi_points = roix*roiy
    lbound = ref_point_id - roi_r*offset - roi_r
  }

  #--
  # bottom edge

  if(roi_r>dist_from_be & roi_r<=dist_from_le & roi_r<=dist_from_re){ 
    roix = roi_r + roi_r + 1
    roiy = dist_from_be + roi_r + 1
    no_roi_points = roix*roiy
    lbound = scan_dist_be - roi_r
  }

  #--
  # top edge

  if(roi_r>dist_from_te & roi_r<=dist_from_le & roi_r<=dist_from_re){
    roix = roi_r + roi_r + 1
    roiy = dist_from_te + roi_r + 1
    no_roi_points = roix*roiy
    lbound = ref_point_id - roi_r*offset - roi_r
  }

  #--
  # left edge

  if(roi_r>dist_from_le & roi_r<=dist_from_be & roi_r<=dist_from_te){
    roix = dist_from_le + roi_r + 1
    roiy = roi_r + roi_r + 1
    no_roi_points = roix*roiy
    lbound = ref_point_id - roi_r*offset - dist_from_le
  }

  #--
  # right edge

  if(roi_r>dist_from_re & roi_r<=dist_from_be & roi_r<=dist_from_te){
    roix = dist_from_re + roi_r + 1
    roiy = roi_r + roi_r + 1
    no_roi_points = roix*roiy
    lbound = ref_point_id - roi_r*offset - roi_r
  }

  #--
  # central area

  if(roi_r<=dist_from_be & roi_r<=dist_from_te &
     roi_r<=dist_from_le & roi_r<=dist_from_re){
    roix = roi_r + roi_r + 1
    roiy = roi_r + roi_r + 1
    no_roi_points = roix*roiy
    lbound = ref_point_id - roi_r*offset - roi_r
  }



  #--


  roi_points_id = getROIPointsID(no_roi_points, offset, roix, lbound)

  return(roi_points_id) 

}



#
# ROI_RID = the Region of Interest is a rectangle (or a square) and it is
#           localized within a infinite domain
#
# Let us consider the neighbourhood (Region of Interest, ROI) of a grid point
# (i.e. data structure). The ROI is a rectangle centred at the grid point and
# containing the neighbour points.
# The width and height of the ROI are (roi_le+roi_re) and (roi_be+roi_te).
# The area of the ROI depends only on the width and height.
#

#
# Given a set of the data structures associated with a ROI, the function
# `evaluateROI_RID()' evaluates the grid locations (IDs of the grid points).
# 
# input:
#   grid_res     = grid resolution (nc x nr)
#   ref_point_id = ID of the reference point (center of the ROI)
#   roi_le       = distance of the reference point from the left edge
#                  of the ROI (in grid units)
#   roi_re       = distance of the reference point from the right edge
#                  of the ROI (in grid units)
#   roi_be       = distance of the reference point from the bottom edge
#                  of the ROI (in grid units)
#   roi_te       = distance of the reference point from the top edge
#                  of the ROI (in grid units)
#
#
evaluateROI_RID = function(grid_res, ref_point_id, roi_le, roi_re, roi_be, roi_te){

  # grid specs
  no_grid_points = grid_res[1] * grid_res[2]
  offset = grid_res[1]
  
  roix = roi_le + roi_re + 1
  roiy = roi_be + roi_te + 1
  no_roi_points = roix*roiy
  lbound = ref_point_id - roi_be*offset - roi_le

  roi_points_id = getROIPointsID(no_roi_points, offset, roix, lbound)

  return(roi_points_id) 

}



#
# ROI_BID = the Region of Interest is a buffer around a line and it is
#           localized within a infinite domain
#
# Let us consider a `buffer' around a line (series of discrete points associated
# with spatial locations);
# The dimensions of the ROI are the distances `roi_le', `roi_re', `roi_be',
# and `roi_te'.
# The area of the ROI depends only on the dimensions of the ROI.
#

#
# Given a set of the data structures associated with a ROI, the function
# `evaluateROI_BID()' evaluates the grid locations (IDs of the grid points).
# 
# input:
#   grid_res    = grid resolution (nc x nr)
#   l_points_id = IDs of the grid points associated with the line
#   roi_le       = distance of the line from the left edge
#                  of the ROI (in grid units)
#   roi_re       = distance of the line from the right edge
#                  of the ROI (in grid units)
#   roi_be       = distance of the line from the bottom edge
#                  of the ROI (in grid units)
#   roi_te       = distance of the line from the top edge
#                  of the ROI (in grid units)
#
#
evaluateROI_BID = function(grid_res, l_points_id, roi_le, roi_re, roi_be, roi_te){

  roi_points_id = evaluateROI_RID(grid_res, l_points_id[1], roi_le, roi_re, roi_be, roi_te)
  for (i in 2:length(l_points_id)) {
    id_temp = evaluateROI_RID(grid_res, l_points_id[i], roi_le, roi_re, roi_be, roi_te)
    roi_points_id = append(roi_points_id, id_temp)
  }

  roi_points_id = unique(roi_points_id)
  roi_points_id = sort(roi_points_id)

  return(roi_points_id) 

}



#--

#
# data structure format:
#        [x, y, p1, p2, p3, id]
# where
#        x  = X coord.
#        y  = Y coord.
#        p1 = parameter 1
#        p2 = parameter 2
#        p3 = parameter 3
#        id = data structure ID
# 
# the function `updateDataStructures()' updates the subset of
# data structures associated with a given ROI by assigning a value
# to the parameter.
#
# input:
#   data_structs   = set of data structures
#   no_ds_el       = no. of elements for each data structure (4<=no_ds_el<=6)
#   roi_points_id  = grid locations (IDs of the grid points) of the
#                    data structures within the ROI
#   value1         = value to assign to the param. 1 (set always)
#   value2         = value to assign to the param. 2 (set to NULL, if no_ds_el<5)
#   value3         = value to assign to the param. 3 (set to NULL, if no_ds_el<6)
#
#
updateDataStructures = function(data_structs, no_ds_el, roi_points_id,
                                value1, value2=NULL, value3=NULL){

  # each data structure contains no_ds_el elements,
  # the last element is the data structure ID

  # no. of data structures
  no_dstrs = dim(data_structs)[1]
  # no. of grid points within the ROI
  no_roi_points = length(roi_points_id)

  ubound = roi_points_id[no_roi_points]
  if (ubound > no_dstrs) {
    ubound = no_dstrs  
  }

  for(i in 1:ubound){
    for(j in 1:no_roi_points){
      if (data_structs[i,no_ds_el] == roi_points_id[j]){
        data_structs[i,3] = value1
        if (length(value2)) data_structs[i,4] = value2
        if (length(value3)) data_structs[i,5] = value3
        break
      }
    }
  }

  return(data_structs)

}



#
# CVS - $Id: update_procedure.R,v 1.2 2010/02/24 13:39:51 lmilazzo Exp $
#



#--
#
# -- conventions --
#
# `s'   = survey
# `p'   = prior
# `h'   = holes
# `bw'  = bandwidth
# 'st'  = stratified
# 'mst' = multiple stratified
# `sl'  = selected
# `m'   = modified
# `u'   = updated
# `sm'  = smoothed
#
#--

#--
#
# -- functions --
#
# evaluatePosterior()
# updateDataSet()
# getExtraPriorPoints()
# searchHoles()
#
# runUpdateProcedure()
# runUpdateProcedure_st()
#
#--





#--

#' Function to evaluate the `posterior' density surfacer; it updates the prior surface with the information associated with the new data
#'
#' @param data_new_y vector of new density data
#' @param data_new_cv vector of coefficients of variation values for the new data (between 0 and 25)
#' @param data_prior_y vector of prior density data
#' @param data_prior_cv vector of coefficients of variation values for the prior data (between 0 and 25)
#'
#' @return
#' Returns a list object: posterior densities, posterior CVs and posterior variances
#' @export
#'
#'@examples
#' x<-seq(0, 10, length=15)
#' y1<-rep(5, 15)
#' cv1<-rep(0.6, 15)
#'
#' y2<- (evaluateGaussianKernel(x, 5, 1)*6) +2
#' cv2<-rep(0.1, 15)
#'
#' posterior<-evaluatePosterior(y2, cv2, y1, cv1)
#'
#' plot(x,y1, ylim=c(0,6), ylab='y', xlab='x')
#' points(x, y2, pch=20)
#' points(x, posterior$post_y, pch=20, col='red')
#'
evaluatePosterior<- function(data_new_y, data_new_cv, data_prior_y, data_prior_cv){

  no_el<- length(data_prior_y)
  post_y<- vector(length=no_el)
  post_sigma<- vector(length=no_el)
  post_cv<- vector(length=no_el)

  for(i in 1:no_el){

    y = data_new_y[i]

    if(is.na(y)==T){
      post_y[i]<- NA
      post_sigma[i]<- NA
      post_cv[i]<- NA
    }
    else{
      if(is.na(data_prior_y[i])==T){
        post_y[i]<- data_new_y[i]
        post_sigma[i]<- NA
        post_cv[i]<- data_new_cv[i]
      }else{
        if(data_prior_cv[i]==0 & data_prior_y[i]==0){
        post_y[i]<- data_new_y[i]
        post_sigma[i]<- NA
        post_cv[i]<- data_new_cv[i]
      }else{
        if(y==0){
          y<- 1e-20
        }
      sigma2logy = log(data_new_cv[i]^2+1)                # variance of log(y) (taken as known)
      # prior mean and cv of lognormal
      mutheta = log(data_prior_y[i]) - 0.5*log(data_prior_cv[i]^2+1)     # prior mean for theta=mean[log(y)]
      sigma2theta = log(data_prior_cv[i]^2+1)               # prior variance for theta=mean[log(y)]

      # calculate posterior on log scale
      taulogy = 1/sigma2logy                     # precision of logy
      tau0 = 1/sigma2theta                       # precision of mean[log(y)]
      taupost = taulogy + tau0                   # posterior precision of mean[log(y)]
      # posterior mean of mean[log(y)]
      mupost = (taulogy*log(y) + tau0*mutheta)/taupost
      # posterior variance of mean[log(y)]
      sigma2post = 1/taupost
      # convert posterior of  mean[log(y)] to posterior of exp{mean[log(y)]} (i.e. posterior on scale of y)
      estpost = fromNorm2logNorm(list(mun=mupost, sigma2n=sigma2post))

      post_y[i]<- estpost$muln
      post_sigma[i]<- estpost$sigma2ln
      post_cv[i]<- sqrt(estpost$sigma2ln)/estpost$muln   # calculate posterior cv of y
    }}}

  }


  return(list(post_y=post_y, post_cv=post_cv, post_sigma=post_sigma))

}



#
# this function updates the data associated with a fine resolution grid
# with the information associated with a lower resolution grid;
#
# inputs:
#   data_full       = data structure containing all information;
#                     data associated with the low resolution grid;
#   data_incomplete = data structure to be updated;
#                     data associated with the fine resolution grid;
#   grid_specs      = nc, nr, and resolution (grid size step)
#   returnID        = if TRUE, grid locations (IDs of the grid points)
#                     are returned
#
updateDataSet<- function(data_full, data_incomplete, grid_specs, returnID=FALSE,DEBUG_MODE=FALSE){

  # grid resolution
  grid_res = c(grid_specs[1], grid_specs[2])

  # no. of grid points
  no_grid_points = grid_specs[1] * grid_specs[2]

  # offset (grid width)
  offset = grid_specs[1]

  #--

  # the main data are grouped into a data structure (ds0) of 5 elements:
  # Longitude, Latitude, PredictedDensity, QualityDensity, ID
  #
  # no. of elements in data structure `ds0'
  ds0_el = 5

  # define a data structure containing a set of data structures `ds0'
  # containing the info associated with the low resolution grid
  no_dstrs_1 = dim(data_full)[1]
  ds1 = array(dim=c(no_dstrs_1, ds0_el))
  ds1[,1] = data_full$Longitude
  ds1[,2] = data_full$Latitude
  ds1[,3] = data_full$PredictedDensity
  ds1[,4] = data_full$QualityDensity
  ds1[,5] = NA

  # evaluating the grid point IDs for the data structures
  for(i in 1:no_dstrs_1){
    ds1[i,ds0_el] = getGridPointID(ds1[i,1], ds1[i,2], grid_specs)
  }

  # define a data structure containing a set of data structures `ds0'
  # to be updated
  no_dstrs_2 = dim(data_incomplete)[1]
  ds2 = array(dim=c(no_dstrs_2, ds0_el))
  ds2[,1] = data_incomplete$Longitude
  ds2[,2] = data_incomplete$Latitude
  ds2[,3] = NA
  ds2[,4] = NA
  ds2[,5] = NA

  # evaluating the grid point IDs for the data structures
  for(i in 1:no_dstrs_2){
    ds2[i,5] = getGridPointID(ds2[i,1], ds2[i,2], grid_specs)
  }

  # let us consider the neighbourhood (Region of Interest, ROI) of a grid point
  # (i.e. data structure). The ROI is a square centred at the grid point and
  # containing the neighbour points. The radius of the ROI is `roi_r' and the
  # ROI area is (2*roi_r+1)*(2*roi_r+1), in grid units.
  #
  # radius of the ROI (in grid units)
  roi_r = round((ds1[2,1]/grid_specs[3] - ds1[1,1]/grid_specs[3])/2)

  #--
  #
  if(DEBUG_MODE==TRUE) {
    cat(paste("\n> DEBUG -\nROI radius:", roi_r, "\n"))
  }
  #--

  for(i in 1:no_dstrs_1){
    # evaluating the grid locations (IDs of the grid points) of the
    # data structures within the ROI
    roi_points_id = evaluateROI_SFD (grid_res, ds1[i,5], roi_r)
    # updating the set of data structures
    ds2 = updateDataStructures(ds2, ds0_el, roi_points_id, ds1[i,3], ds1[i,4])
  }

  # if the prior resolution is different from the survey resolution, it is
  # necessary to apply a correction factor to the Quality Densities
  if(res_p!=res_s){
    ds2[,4]<- ds2[,4] * sqrt((roi_r*2)^2)
  }

  if(returnID==F) return(list(ds2=ds2))
  else return(list(ds2=ds2, ds1_id=ds1[,ds0_el], ds2_id=ds2[,ds0_el]))

}



#
# this function is used to get extra prior points when the survey runs over
# the edge of the prior surface. Extra points are identified, added to the
# prior data and assigned values of 10^-15 for density (cannot be zero as
# issue with  evaluating the posterior) and 80 for CV (very uncertain as we do
# not know that the surface should be zero here).
#
# inputs:
#  data_u     = original prior data
#  data_s_m   = survey data with range outside prior data
#  grid_specs = (nc, nr, res) of virtual grid
#
getExtraPriorPoints<- function(data_u, data_s_m, grid_specs, DEBUG_MODE=FALSE){

  data_u_temp = na.omit(data_u)

  # is the range of the survey greater than the range of the prior

  range_data_s = rbind(range(na.omit(data_s_m$Longitude)), range(na.omit(data_s_m$Latitude)))
  range_data_p = rbind(range(na.omit(data_u_temp$Longitude)), range(na.omit(data_u_temp$Latitude)))

  if(range_data_p[1,1] > range_data_s[1,1]) {le = 1}else {le = 0}
  if(range_data_p[1,2] < range_data_s[1,2]) {re = 1}else {re = 0}
  if(range_data_p[2,1] > range_data_s[2,1]) {be = 1}else {be = 0}
  if(range_data_p[2,2] < range_data_s[2,2]) {te = 1}else {te = 0}

  # if so then take survey points that are outside prior
  if(any(c(le,re,be,te)==1)){
    overlap <<- T
    extra_edge = which(c(le,re,be,te)==1)
    if(extra_edge[1]==1){
      data_extra = data_s_m[which(data_s_m$Longitude < range_data_p[1,1]),]
    }
    if(extra_edge[1]==2){
      data_extra = data_s_m[which(data_s_m$Longitude > range_data_p[1,2]),]
    }
    if(extra_edge[1]==3){
      data_extra = data_s_m[which(data_s_m$Latitude < range_data_p[2,1]),]
    }
    if(extra_edge[1]==4){
      data_extra = data_s_m[which(data_s_m$Latitude > range_data_p[2,2]),]
    }
    if(length(extra_edge)>1){
      for(i in 2:length(extra_edge)){
        if(extra_edge[i]==2){
          data_extra = rbind(data_extra,data_s_m[which(data_s_m$Longitude > range_data_p[1,2]),])
        }
        if(extra_edge[i]==3){
          data_extra = rbind(data_extra,data_s_m[which(data_s_m$Latitude < range_data_p[2,1]),])
        }
        if(extra_edge[i]==4){
          data_extra = rbind(data_extra,data_s_m[which(data_s_m$Latitude > range_data_p[2,2]),])
        }
      }
    }


    #--
    # -- plots --
    #
    if(DEBUG_MODE==TRUE) {
      plot(data_u$Lon,data_u$Lat,
           xlim = c(min(range(na.omit(data_s_m$Longitude))[1], range(na.omit(data_u$Longitude))[1]),
           max(range(na.omit(data_s_m$Longitude))[2], range(na.omit(data_u$Longitude))[2])),
           ylim= c(min(range(na.omit(data_s_m$Latitude))[1], range(na.omit(data_u$Latitude))[1]),
           max(range(na.omit(data_s_m$Latitude))[2], range(na.omit(data_u$Latitude))[2])))
      points(data_s_m$Lon,data_s_m$Lat,pch=20)
      points(data_extra$Lon,data_extra$Lat,pch=20,col='red')
      cat("\n> DEBUG - plotting data ... done\n")
    }

    #--


    # find virtual id's of points that are extra
    VirtID = vector(length=dim(data_extra)[1])
    for(i in 1:dim(data_extra)[1]){
      VirtID[i] = getGridPointID(data_extra$Longitude[i],data_extra$Latitude[i],grid_specs)
    }

    # if id's equal in virtID and data u then replace density and QD in data_u with those from data_extra
    for(i in 1:dim(data_u)[1]){
      for(j in 1:length(VirtID)){
        if(data_u$VirtID[i] == VirtID[j]){
          data_u$PredictedDensity[i] = 10^(-15)
          data_u$QualityDensity[i] = 80
        }
      }
    }

  }else{
    overlap <<- F
  }


  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    plot(data_u$Lon,data_u$Lat)
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--


  return(data_u)

}



#
# this function search for holes in the `prior' data;
# if there are holes, a data structure associated with the holes is generated
# (it contains `NAs') and the holes data are appended to the prior data;
#
# inputs:
#   data_p_m  = prior data (modified, after coord. transformation)
#   range_lon = max range of longitude between prior and survey
#   range_lat = max range of latitude between prior and survey
#
searchHoles<- function(data_p_m, range_lon, range_lat){

  #
  # project the prior data to a `virtual' grid ...
  #
  # (note that the grid used here is smaller than the grid used for
  # the main calculations)
  #

  # grid specifications: nc, nr, and resolution (grid size step)
  grid_specs<- c((round((range_lon-1)/res_p))+1,
                 (round((range_lat-1)/res_p))+1,
                 res_p)

  no_el = dim(data_p_m)[1]
  # define a data structure for the grid locations (IDs of the grid points)
  # associated with the prior data
  gpoint_id_p_m<- vector(length=no_el)

  # evaluate the the grid locations associated with the prior data
  for(i in 1:no_el){
      gpoint_id_p_m[i]<- getGridPointID(data_p_m$Longitude[i], data_p_m$Latitude[i], grid_specs)
  }

  # define a data structure for the grid locations (IDs of the grid points)
  # associated with the holes
  gpoint_id_h<- NULL

  # evaluate the the grid locations associated with the holes
  for (r in 1:(no_el-1)){
    if(r==1) {
      gpoint_id_h<- NULL
    }
    if((gpoint_id_p_m[r+1]-gpoint_id_p_m[r])>1 & (gpoint_id_p_m[r+1]-gpoint_id_p_m[r])<grid_specs[1]){
      gpoint_id_h<- c(gpoint_id_h, seq((gpoint_id_p_m[r]+1), (gpoint_id_p_m[r+1]-1), by=1))
    }
  }


  # if there are not holes in the prior data ...
  if(is.null(gpoint_id_h)){

    return(list(data_p_m=data_p_m, gpoint_id_h=NULL))

  }
  # if there are holes in the prior data ...
  else{

    if(gpoint_id_p_m[length(gpoint_id_p_m)]!=(grid_specs[1]*grid_specs[2])){
      gpoint_id_h<- c(gpoint_id_h, (gpoint_id_p_m[length(gpoint_id_p_m)]+1):(grid_specs[1]*grid_specs[2]))
    }

    no_el_h = length(gpoint_id_h)
    # define a data structure for the holes;
    # the d.s. has the same format of the d.s. associated with the prior data;
    data_h<- matrix(NA, no_el_h, 9)
    data_h[,1]<- data_p_m[1,1]
    data_h[,2]<- NA
    for(i in 1:no_el_h){
      data_h[i,4:3]<- getCoords(gpoint_id_h[i], grid_specs)
    }
    data_h[,5]<- data_p_m[1,5]
    data_h[,6:9]<- NA
    data_h<- data.frame(data_h)
    names(data_h)<- c('Species_ID', 'AreaCode', 'Latitude', 'Longitude', 'TimePeriod',
                      'PredictedDensity', 'QualityDensity', 'RES','QualityRES')

    # append the holes data to the prior data
    # (note that the holes data contain `NAs')
    data_p_m_u<- rbind(data_p_m, data_h)
    data_p_m_u<- data_p_m_u[order(data_p_m_u$Longitude),]
    data_p_m_u<- data_p_m_u[order(data_p_m_u$Latitude),]

    return(list(data_p_m=data_p_m_u, gpoint_id_h=gpoint_id_h))

  }

}




#' Main function to run the update and smoothing of two surfaces.
#'
#' @param data_s The survey data (new data)
#' @param data_p The prior data (existing surface, e.g. RES surface)
#' @param filepath Pathway for the csv files to be stored. defaults to the working directory.
#' @param DEBUG_MODE Logical stating whether additional plots are returned.
#'
#' @return
#' various plots, two csv files into the working directory. \code{densities_sm} gives the updated survey information in the region of the survey (plus some buffer). \code{densities_sm} gives the updated and smoothed survey embedded in the prior surface.
#'
#' @export
#'
#' @author Lindesay Scott-Hayward
#' @author Lorenzo Milazzo
#'
#' @examples
#' data(data_s) # survey data
#' data(data_p) # prior data
#'
#' runUpdateProcedure(data_s, data_p)
#'
runUpdateProcedure<- function(data_s, data_p, filepath='./', DEBUG_MODE=FALSE, filetag=''){

  library(fields)
  library(splancs)

  # checking the order of the data
  if((data_s$Longitude[2]-data_s$Longitude[1])==0){
    cat("\n> warning - survey is ordered wrongly\n")
  }
  if((data_p$Longitude[2]-data_p$Longitude[1])==0){
    cat("\n> warning - prior is ordered wrongly\n")
    data_p<- data_p[order(data_p$Latitude),]
  }
  if((data_p$Longitude[2]-data_p$Lon[1])<0){
    cat("\n> warning - prior is ordered wrongly\n")
    data_p<- data_p[order(data_p$Lon),]
    data_p<- data_p[order(data_p$Latitude),]
  }


  #--
  # -- plots --
  #
  quilt.plot(data_p$Longitude, data_p$Latitude, data_p$PredictedDensity,
             nrow=length(unique(data_p$Longitude)), ncol=length(unique(data_p$Latitude)),
             zlim=c(min(range(na.omit(data_s$PredictedDensity))[1], range(na.omit(data_p$PredictedDensity))[1]),
                    max(range(na.omit(data_s$PredictedDensity))[2], range(na.omit(data_p$PredictedDensity))[2])),
             main="prior and survey",
             xlim= c(min(range(na.omit(data_s$Longitude))[1], range(na.omit(data_p$Longitude))[1]),
                    max(range(na.omit(data_s$Longitude))[2], range(na.omit(data_p$Longitude))[2])),
             ylim= c(min(range(na.omit(data_s$Latitude))[1], range(na.omit(data_p$Latitude))[1]),
                    max(range(na.omit(data_s$Latitude))[2], range(na.omit(data_p$Latitude))[2])))
  quilt.plot(data_s$Longitude,data_s$Latitude, data_s$PredictedDensity,
             nrow=length(unique(data_s$Longitude)), ncol=length(unique(data_s$Latitude)),
             zlim=c(min(range(na.omit(data_s$PredictedDensity))[1], range(na.omit(data_p$PredictedDensity))[1]),
                    max(range(na.omit(data_s$PredictedDensity))[2], range(na.omit(data_p$PredictedDensity))[2])),
             add=T)
  points(data_s$Lon, data_s$Lat, pch=20, cex=0.5)
  cat("\n++ plotting data ... done\n")

  #--


  # check for 0 or NA values
  if(length(which(data_s$PredictedDensity==0))>=1){
    cat("\n> warning - there are 0's in survey\n")
  }
  if(length(which(is.na(data_s)))){
    cat("\n> warning - there are NA's in survey data\n")
  }

  # check to see if the area under analysis crosses the 180/-180 line
  # and change accordingly
  if((length(which(data_s$Longitude >=-180 & data_s$Longitude <=-179)) +
      length(which(data_s$Longitude<=180 & data_s$Longitude >=179))) >
      length(which(data_s$Longitude >=-180 & data_s$Longitude <=-179))) {

    cat("\n> warning - the area crosses the 180/-180 line\n")

    # make all values positive to stitch world back together for survey and prior
    for (i in 1:length(data_s$Longitude)){
      if(data_s$Longitude[i]<0) {
        data_s$Longitude[i]<- data_s$Longitude[i] + 360
      }
    }
    for (i in 1:length(data_p$Longitude)){
      if(data_p$Longitude[i]<0) {
        data_p$Longitude[i]<-data_p$Longitude[i] + 360
      }
    }
    # sort data again so in order (bottom left as start)
    data_p<- data_p[order(data_p$Longitude),]
    data_p<- data_p[order(data_p$Latitude),]
    data_s<- data_s[order(data_s$Longitude),]
    data_s<- data_s[order(data_s$Latitude),]


    #--
    # -- plots --
    #
    quilt.plot(data_p$Longitude, data_p$Latitude, data_p$PredictedDensity,
               nrow=length(unique(data_p$Longitude)), ncol=length(unique(data_p$Latitude)),
               zlim=c(min(range(data_s$PredictedDensity)[1], range(data_p$PredictedDensity)[1]),
                      max(range(data_s$PredictedDensity)[2], range(data_p$PredictedDensity)[2])),
               main='180 adjustment')
    quilt.plot(data_s$Longitude, data_s$Latitude,data_s$PredictedDensity,
               nrow=length(unique(data_s$Longitude)), ncol=length(unique(data_s$Latitude)),
               zlim=c(min(range(data_s$PredictedDensity)[1], range(data_p$PredictedDensity)[1]),
                      max(range(data_s$PredictedDensity)[2], range(data_p$PredictedDensity)[2])),
               add=T)
    cat("\n++ plotting 180/-180 data ... done\n")

    #--

  }


  # select a subset of prior data that corresponds to a region
  # including the survey area;
  # more precisely ...
  # . extend the survey area by adding 4 degrees in all directions;
  # . select the prior data that corresponds to this region;
  range_s<- rbind(range(data_s$Longitude), range(data_s$Latitude))
  range_p<- range_s
  range_p[,1]<- range_s[,1] - 4
  range_p[,2]<- range_s[,2] + 4
  data_p_sl<- data_p[data_p$Longitude>=range_p[1,1] &
                     data_p$Longitude<=range_p[1,2] &
                     data_p$Latitude>=range_p[2,1] &
                     data_p$Latitude<=range_p[2,2],]
  range_p<- rbind(range(data_p_sl$Longitude), range(data_p_sl$Latitude))


  #
  # projecting the data (selected prior data and survey data) to a `virtual' grid ...
  #
  # Let us associate each location with a grid point.
  # The locations are associated with the grid points row by row, starting from the
  # bottom left corner and ending to the top right corner of the grid.
  #
  # The grid point at the bottom left corner is the first grid point (ID = 1)
  # and it is associated with the spatial location (1,1).
  # Note that in general a `coordinate transformation' is needed.
  #

  # coordinate transformation; translation factors
  k_lon<- (1 - min(data_p_sl$Longitude[1], data_s$Longitude[1]))
  k_lat<- (1 - min(data_p_sl$Latitude[1], data_s$Latitude[1]))

  data_p_m<- data_p_sl
  data_p_m$Longitude<- data_p_sl$Longitude + k_lon
  data_p_m$Latitude<- data_p_sl$Latitude + k_lat
  # resolution (grid size step) associated with the prior
  gridsteps<-dist(cbind(data_p_m$Longitude, data_p_m$Latitude))
  res_p_temp<- min(gridsteps[gridsteps>0])
  # store the prior resolution into a global variable
  res_p<<- res_p_temp

  data_s_m<- data_s
  data_s_m$Longitude<- data_s$Longitude + k_lon
  data_s_m$Latitude<- data_s$Latitude + k_lat
  data_s_m$Latitude<- round(data_s_m$Latitude, 4)
  data_s_m$Longitude<- round(data_s_m$Longitude, 4)
  # resolution (grid size step) associated with the survey
  res_s_temp<- min(dist(cbind(data_s_m$Longitude, data_s_m$Latitude)))
  # store the survey resolution into a global variable
  res_s<<- res_s_temp

  # minimum value of the resolution (grid size step)
  res_min_temp<- min(res_s, res_p)
  # store the minimum resolution into a global variable
  res_min<<- res_min_temp

  #--
  #
  if(DEBUG_MODE==TRUE) {
    cat(paste("\n> DEBUG -\nresolutions (grid size steps):", res_p, "(prior),", res_s, "(survey),", res_min, "(minimum value)\n"))
  }

  #--


  #
  # note that here the term `resolution' refers to the `grid size step'
  # (or `grid size');
  # a low value of the `resolution' corresponds to a high resolution grid;
  # a high value of the `resolution' corresponds to a low resolution grid;
  #

  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    quilt.plot(data_s_m$Longitude, data_s_m$Latitude, data_s_m$PredictedDensity,
               nrow=length(unique(data_s_m$Longitude)), ncol=length(unique(data_s_m$Latitude)),
               main='Survey on virtual grid')
    points(data_s_m$Longitude, data_s_m$Latitude)
    quilt.plot(data_p_m$Longitude, data_p_m$Latitude, data_p_m$PredictedDensity,
               nrow=length(unique(data_p_m$Longitude)), ncol=length(unique(data_p_m$Latitude)),
               main='Prior on virtual grid')
    points(data_p_m$Longitude,data_p_m$Latitude)
    points(data_s_m$Longitude, data_s_m$Latitude, pch=20)
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--


  range_lon<- max(range(data_p_m$Longitude)[2], range(data_s_m$Longitude)[2])
  range_lat<- max(range(data_p_m$Latitude)[2], range(data_s_m$Latitude)[2])

  # search for holes in the prior data (the holes indicate `islands')
  ds_h<- searchHoles(data_p_m, range_lon, range_lat)
  data_p_m<- ds_h$data_p_m

  # grid specifications: nc, nr, and resolution (grid size step)
  grid_specs<- c((round((range_lon)/res_min)*2),
                 (round((range_lat)/res_min)*2),
                 res_min/2)

  #--
  #
  if(DEBUG_MODE==TRUE) {
    cat(paste("\n> DEBUG -\ngrid specifications:",
               grid_specs[1], grid_specs[2], grid_specs[3], "\n"))
  }
  #--


  # redefine the data structure by making space for the IDs
  data_p_m<- cbind(data_p_m, VirtID=vector(length=dim(data_p_m)[1]))
  # evaluate the grid locations (IDs of the grid points)
  # associated with the prior data
  for(i in 1:dim(data_p_m)[1]){
    data_p_m$VirtID[i] = getGridPointID(data_p_m$Longitude[i], data_p_m$Latitude[i], grid_specs)
  }

  # grid point IDs associated with the holes within
  # the prior area
  gpoint_id_h<- data_p_m$VirtID[which(is.na(data_p_m$PredictedDensity))]

  # evaluate the grid locations (IDs of the grid points)
  # associated with the holes within the prior area
  if(is.null(gpoint_id_h)==F){
    gpoint_id_h<- data_p_m$VirtID[which(is.na(data_p_m$PredictedDensity))]
  }

  # redefine the data structure by making space for the IDs
  data_s_m<- cbind(data_s_m, VirtID=vector(length=dim(data_s_m)[1]))
  # evaluate the grid locations (IDs of the grid points)
  # associated with the survey data
  for(i in 1:dim(data_s_m)[1]){
    data_s_m$VirtID[i] = getGridPointID(data_s_m$Longitude[i], data_s_m$Latitude[i], grid_specs)
  }


  if(res_s==res_p){
    cat("\n> the grid size step of the survey is equal to the grid size step\n  of the prior (the survey has an equal resolution)\n")
    data_full<- data_p_m
    data_incomplete<- data_s_m
  }

  if(res_s<res_p){
    cat("\n> the grid size step of the survey is less than the grid size step\n  of the prior (the survey has a higher resolution)\n")
    data_full<- data_p_m
    data_incomplete<- data_s_m
  }

  if(res_s>res_p){
    cat("\n> the grid size step of the survey is greater than the grid size step\n  of the prior (the survey has a lower resolution)\n")
    data_full<- data_s_m
    data_incomplete<- data_p_m
  }

  #if(res_s<res_p | res_s>res_p){
    # update the data associated with a fine resolution grid with the information
    # associated with a lower resolution grid
    update<- updateDataSet(data_full, data_incomplete, grid_specs, returnID=T, DEBUG_MODE=DEBUG_MODE)
    data_u<- data.frame(update$ds2)
    names(data_u)<- c('Longitude', 'Latitude', 'PredictedDensity', 'QualityDensity','VirtID')
 # }


  data_u<- getExtraPriorPoints(data_u,data_s_m, grid_specs, DEBUG_MODE=DEBUG_MODE)

  # if there are holes in the survey, translate this to holes in the updated prior
  if(length(which(is.na(data_s$PredictedDensity)))>0){
    gpoint_id_h<- unique(c(gpoint_id_h, update$ds2_id[which(is.na(data_s$PredictedDensity))]))
    for(i in 1:length(gpoint_id_h)){
      for(j in 1:dim(data_u)[1]){
        if(gpoint_id_h[i]==update$ds2_id[j]){
          data_u$PredictedDensity[j]<- NA
          data_u$QualityDensity[j]<- NA
        }
      }
    }
    gpoint_id_h<- sort(gpoint_id_h)
  }else{
    gpoint_id_h<- gpoint_id_h
  }


  # define `new' and `prior' data
  if(res_s<res_p){
    data_new<- data_s_m
    data_prior<- data_u
    gpoint_id<- update$ds2_id
  }
  if(res_s>res_p){
    data_new<- data_u
    data_prior<- data_s_m
    gpoint_id<- update$ds1_id
  }

  if(res_s==res_p){
    data_new<- data_s_m
    data_prior<- data_u
    gpoint_id<- update$ds2_id
  }

  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    quilt.plot(data_new$Longitude, data_new$Latitude, data_new$PredictedDensity,
               nrow=length(unique(data_new$Longitude)), ncol=length(unique(data_new$Latitude)),
               zlim=range(na.omit(data_new$PredictedDensity)),
               main="updated survey resolution surface")
    quilt.plot(data_prior$Longitude, data_prior$Latitude, data_prior$PredictedDensity,
               nrow=length(unique(data_prior$Longitude)), ncol=length(unique(data_prior$Latitude)),
               zlim=range(na.omit(data_u$PredictedDensity)),
               main="updated prior resolution surface")
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--


  # evaluate `posterior' by combining `prior' and `new' data;
  # update the `prior' surface with the `new' data;
  data_post<- evaluatePosterior(data_new$PredictedDensity, data_new$QualityDensity,
                                data_prior$PredictedDensity, data_prior$QualityDensity)
  cat("\n++ evaluating `posterior' by combining `prior' and `new' data ... done\n")

  if(any(na.omit(data_post$post_y)<0)) {
    cat("\n> warning - zeros in posterior\n")
  }

  # check for updated densities greater or less than max/min of inputs
  if (max(na.omit(data_post$post_y)) > max(na.omit(data_p$PredictedDensity), na.omit(data_s$PredictedDensity)) |
      min(na.omit(data_post$post_y)) < min(na.omit(data_p$PredictedDensity), na.omit(data_s$PredictedDensity))) {
    cat("\n> warning - updated data outside range of inputs\n")
  }


  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {

    if(is.null(gpoint_id_h)){
      quilt.plot(data_p_m$Longitude, data_p_m$Latitude, data_p_m$PredictedDensity,
                 nrow=length(unique(data_p_m$Longitude)), ncol=length(unique(data_p_m$Latitude)),
                 zlim=c(min(range(na.omit(data_prior$PredictedDensity))[1], range(data_p_sl$PredictedDensity)[1]),
                        max(range(na.omit(data_prior$PredictedDensity))[2], range(data_p_sl$PredictedDensity)[2])),
                 main="updated surface on prior")
      quilt.plot(data_prior$Longitude, data_prior$Latitude, data_post$post_y,
                 nrow=length(unique(data_s_m$Longitude)), ncol=length(unique(data_s_m$Latitude)),
                 zlim=c(min(range(na.omit(data_prior$PredictedDensity))[1], range(data_p_sl$PredictedDensity)[1]),
                        max(range(na.omit(data_prior$PredictedDensity))[2], range(data_p_sl$PredictedDensity)[2])),
                 add=T)
    }else{

      if(length(which(is.na(data_p_m$PredictedDensity)))>0){
        n_row<- length(unique(data_p_m$Longitude[-which(is.na(data_p_m$PredictedDensity))]))
        n_col<- length(unique(data_p_m$Latitude[-which(is.na(data_p_m$PredictedDensity))]))
      }else{
        n_row=length(unique(data_p_m$Longitude))
        n_col=length(unique(data_p_m$Latitude))
      }
      quilt.plot(data_p_m$Longitude, data_p_m$Latitude, data_p_m$PredictedDensity,
                 nrow=n_row, ncol=n_col,
                 zlim=c(min(range(na.omit(data_prior$PredictedDensity))[1], range(na.omit(data_p_m$PredictedDensity))[1]),
                        max(range(na.omit(data_prior$PredictedDensity))[2], range(na.omit(data_p_m$PredictedDensity))[2])),
                 xlim= c(min(range(na.omit(data_new$Longitude))[1], range(na.omit(data_p_m$Longitude))[1]),
                         max(range(na.omit(data_new$Longitude))[2], range(na.omit(data_p_m$Longitude))[2])),
                 ylim= c(min(range(na.omit(data_new$Latitude))[1], range(na.omit(data_p_m$Latitude))[1]),
                         max(range(na.omit(data_new$Latitude))[2], range(na.omit(data_p_m$Latitude))[2])),
                 main="updated surface on prior")

     if(length(which(is.na(data_prior$PredictedDensity)))>0){
        n_row<- length(unique(data_prior$Longitude[-which(is.na(data_prior$PredictedDensity))]))
        n_col<- length(unique(data_prior$Latitude[-which(is.na(data_prior$PredictedDensity))]))
     }else{
        n_row=length(unique(data_prior$Longitude))
        n_col=length(unique(data_prior$Latitude))
     }
     quilt.plot(data_prior$Longitude, data_prior$Latitude, data_post$post_y,
                nrow=n_row, ncol=n_col,
                zlim=c(min(range(na.omit(data_prior$PredictedDensity))[1], range(na.omit(data_p_m$PredictedDensity))[1]),
                       max(range(na.omit(data_prior$PredictedDensity))[2], range(na.omit(data_p_m$PredictedDensity))[2])),
                add=T)
    }

    points(data_new$Longitude, data_new$Latitude, pch='x', cex=0.5)
    points(data_new$Longitude[which(is.na(data_post$post_y))],
           data_new$Latitude[which(is.na(data_post$post_y))], pch=20)
    cat("\n> DEBUG - plotting data ... done\n")

  }

  #--


  # updated densities (posterior)
  densities_u<- data.frame(data_prior$Latitude, data_prior$Longitude,
                           data_post$post_y, data_post$post_cv)
  names(densities_u)<- c('Latitude', 'Longitude',

                         'PredictedDensity', 'QualityDensity')

  # coordinate transformation; convert back from virtual grid
  densities_u$Longitude<- densities_u$Longitude - k_lon
  densities_u$Latitude<- densities_u$Latitude - k_lat
  # convert back to degrees and if at 180, make appropriate -/+ split
  if(range(densities_u$Longitude)[2]>180){
    for(i in 1:length(densities_u$Longitude)){
      if(densities_u$Longitude[i]>180) densities_u$Longitude[i]<- (densities_u$Longitude[i] - 360)
    }
  }

  #
  # output processing
  #
  # output:
  #  . updated densities
  write.csv(densities_u, file=paste(filepath, 'densities_u_', filetag ,'.csv', sep=''), row.names=FALSE)
  cat("\n++ writing output (updated densities) ... done\n")


  #--

  # selecting the bandwidth
  data_bw<- selectBandwidth(grid_specs, k_lon, k_lat,
                            range_lon, range_lat,
                            data_full, data_prior, data_new, data_post,
                            data_s, data_p_m,
                            gpoint_id, gpoint_id_h, DEBUG_MODE=DEBUG_MODE)
  cat("\n++ selecting the bandwidth ... done\n")

  #--

  gpoints<- matrix(as.matrix(data_bw[,1:2]), ncol=2)
  # data to be smoothed
  data_y<- data_bw$PredictedDensity
  data_h<- cbind(data_bw$h1, data_bw$h2)
  data_id<- data_bw$id

  # grid point IDs associated with the holes within
  # the area to be smoothed
  # (note that `gpoint_id_h_sm' is a subset of `gpoint_id_h')
  gpoint_id_h_sm<- data_id[which(is.na(data_y))]
  if(length(gpoint_id_h_sm)<=1){
    gpoint_id_h_sm<- NULL
  }

  # doing the kernel smoothing
  data_sm_y<- doKernelSmoothing(gpoints, data_y, data_h, data_id, gpoint_id_h_sm)
  cat("\n++ doing the kernel smoothing ... done\n")

  # check for densities produced less than zero
  if(any(na.omit(data_sm_y) < 0)){
    cat("\n> warning - smoothed densities less than zero\n")
    for (i in 1:length(data_sm_y)){
      if(is.na(data_sm_y[i]) == F){
        if(data_sm_y[i] < 0){
          data_sm_y[i] <- 10^(-15)
        }
      }
    }
  }

  #--
  # -- plots --
  #
  if(is.null(gpoint_id_h)==T){
    quilt.plot(data_p_m$Longitude, data_p_m$Latitude, data_p_m$PredictedDensity,
               nrow=length(unique(data_p_m$Lon)), ncol=length(unique(data_p_m$Lat)),
               zlim=c(min(range(na.omit(data_p_m$PredictedDensity))[1], range(na.omit(data_sm_y))[1]),
                      max(range(na.omit(data_p_m$PredictedDensity))[2], range(na.omit(data_sm_y))[2])),
               xlim=range(data_p_m$Longitude), ylim=range(data_p_m$Latitude),
               main="smoothed surface on prior")
    quilt.plot(data_bw$Longitude, data_bw$Latitude, data_sm_y,
               nrow=length(unique(data_bw$Lon)), ncol=length(unique(data_bw$Lat)),
               zlim=c(min(range(na.omit(data_p_m$PredictedDensity))[1], range(na.omit(data_sm_y))[1]),
                      max(range(na.omit(data_p_m$PredictedDensity))[2], range(na.omit(data_sm_y))[2])),
               add=T)
  }else{
    if(length(which(is.na(data_bw$PredictedDensity)))>0){
      n_row_b<- length(unique(data_bw$Longitude[-which(is.na(data_bw$PredictedDensity))]))
      n_col_b<- length(unique(data_bw$Latitude[-which(is.na(data_bw$PredictedDensity))]))
    }else{
      n_row_b=length(unique(data_bw$Longitude))
      n_col_b=length(unique(data_bw$Latitude))
    }
    quilt.plot(data_bw$Longitude, data_bw$Latitude, data_sm_y,
               nrow=n_row_b, ncol=n_col_b,
               main="smoothed surface (survey area and buffer edge)")
    if(length(which(is.na(data_p_m$PredictedDensity)))>0){
      n_row<- length(unique(data_p_m$Longitude[-which(is.na(data_p_m$PredictedDensity))]))
      n_col<- length(unique(data_p_m$Latitude[-which(is.na(data_p_m$PredictedDensity))]))
    }else{
      n_row=length(unique(data_p_m$Longitude))
      n_col=length(unique(data_p_m$Latitude))
    }
    quilt.plot(data_p_m$Longitude, data_p_m$Latitude, data_p_m$PredictedDensity,
               nrow=n_row, ncol=n_col,
               zlim=c(min(range(na.omit(data_p_m$PredictedDensity))[1], range(na.omit(data_sm_y))[1]),
                      max(range(na.omit(data_p_m$PredictedDensity))[2], range(na.omit(data_sm_y))[2])),
               xlim=c(1,range_lon), ylim=c(1,range_lat),
               main="smoothed surface on prior")
    quilt.plot(data_bw$Lon, data_bw$Lat, data_sm_y,
               nrow=n_row_b, ncol=n_col_b,
               zlim=c(min(range(na.omit(data_p_m$PredictedDensity))[1], range(na.omit(data_sm_y))[1]),
                      max(range(na.omit(data_p_m$PredictedDensity))[2], range(na.omit(data_sm_y))[2])),
               add=T)
  }
  polymap(cbind(data_b$Longitude+k_lon, data_b$Latitude+k_lat), add=T)
  cat("\n++ plotting data ... done\n")

  #--


  data_sm<- data_bw[,1:4]
  # smoothed data
  data_sm$PredictedDensity<- data_sm_y
  # coordinate transformation; convert back from virtual grid
  data_sm$Longitude<- data_sm$Longitude - k_lon
  data_sm$Latitude<- data_sm$Latitude - k_lat
  # convert back to degrees and if at 180, make appropriate -/+ split
  if(range(data_sm$Longitude)[2]>180){
    for(i in 1:length(data_sm$Longitude)){
      if(data_sm$Longitude[i]>180) data_sm$Longitude[i]<- (data_sm$Longitude[i] - 360)
    }
  }

  densities_sm<- data.frame(data_sm$Latitude, data_sm$Longitude,
                            data_sm$PredictedDensity, data_sm$QualityDensity)
  names(densities_sm)<- c('Latitude', 'Longitude', 'PredictedDensity', 'QualityDensity')

  #
  # output processing
  #
  # output:
  #  . smoothed densities
  write.csv(densities_sm, file=paste(filepath, 'densities_sm_',filetag,'.csv', sep=''), row.names=FALSE)
  cat("\n++ writing output (smoothed densities) ... done\n")


}


#
#        --   main function (stratified data)   --
#
runUpdateProcedure_st<- function(data_s_st, data_b, data_p){

  #
  # note that, in the case of `stratified' data, a data set for
  # the survey has to be generated;
  #


  # check that prior is ordered correctly so that survey is made correctly too
  if((data_p$Longitude[2]-data_p$Lon[1])<0){
    cat("\n> warning - prior is ordered wrongly\n")
    data_p<- data_p[order(data_p$Lon),]
    data_p<- data_p[order(data_p$Latitude),]
  }


  # select a subset of prior data that corresponds to a region
  # including the survey area;
  # more precisely ...
  # . use the data from the boundaries file and define the survey area;
  # . extend the survey area by adding 3 degrees in all directions;
  # . select the prior data that corresponds to this region;
  range_s<- rbind(range(data_b$Longitude), range(data_b$Latitude))
  range_p<- range_s
  range_p[,1]<- range_s[,1] - 3
  range_p[,2]<- range_s[,2] + 3
  data_p_sl<- data_p[data_p$Longitude>=range_p[1,1] &
                     data_p$Longitude<=range_p[1,2] &
                     data_p$Latitude>=range_p[2,1] &
                     data_p$Latitude<=range_p[2,2],]
  range_p<- rbind(range(data_p_sl$Longitude), range(data_p_sl$Latitude))


  # check for 0 or NA values
  if(length(which(data_p_sl$PredictedDensity==0))>=1){
    cat("\n> warning - there are 0's in prior\n")
  }
  if(length(which(is.na(data_p_sl)))){
    cat("\n> warning - there are NA's in prior data\n")
  }

  # check to see if the area under analysis crosses the 180/-180 line
  # and change accordingly
  if((length(which(data_p_sl$Longitude >=-180 & data_p_sl$Longitude <=-179)) +
      length(which(data_p_sl$Longitude<=180 & data_p_sl$Longitude >=179))) >
      length(which(data_p_sl$Longitude >=-180 & data_p_sl$Longitude <=-179))){

    cat("\n> warning - the area crosses the 180/-180 line\n")

    # make all values positive to stitch world back together for prior
    for (i in 1:length(data_p_sl$Longitude)){
      if(data_p_sl$Longitude[i]<0) {
        data_p_sl$Longitude[i]<- data_p_sl$Longitude[i] + 360
      }
    }

    # sort data again so in order (bottom left as start)
    data_p_sl<- data_p[order(data_p_sl$Longitude),]
    data_p_sl<- data_p[order(data_p_sl$Latitude),]


    #--
    # -- plots --
    #
    quilt.plot(data_p_sl$Longitude, data_p_sl$Latitude, data_p_sl$PredictedDensity,
               nrow=length(unique(data_p_sl$Longitude)), ncol=length(unique(data_p_sl$Latitude)),
               zlim=range(data_p_sl$PredictedDensity),
               main='180 adjustment')
    cat("\n++ plotting 180/-180 data ... done\n")

    #--

  }


  #--
  # -- plots --
  #
  quilt.plot(data_p$Longitude, data_p$Latitude, data_p$PredictedDensity,
             nrow=length(unique(data_p$Longitude)), ncol=length(unique(data_p$Latitude)),
             zlim= range(na.omit(data_p$PredictedDensity)),
             main="prior and survey area")
  polymap(data_b[,3:2], add=T)
  cat("\n++ plotting data ... done\n")

  #--


  #
  # projecting the data (only the selected prior data) to a `virtual' grid ...
  #
  # Let us associate each location with a grid point.
  # The locations are associated with the grid points row by row, starting from the
  # bottom left corner and ending to the top right corner of the grid.
  #
  # The grid point at the bottom left corner is the first grid point (ID = 1)
  # and it is associated with the spatial location (1,1).
  # Note that in general a `coordinate transformation' is needed.
  #

  # coordinate transformation; translation factors
  k_lon<- (1 - data_p_sl$Longitude[1])
  k_lat<- (1 - data_p_sl$Latitude[1])

  data_p_m<- data_p_sl
  data_p_m$Longitude<- data_p_sl$Longitude + k_lon
  data_p_m$Latitude<- data_p_sl$Latitude + k_lat
  # resolution (grid size step) associated with the prior
  res_p_temp<- min(dist(cbind(data_p_m$Longitude, data_p_m$Latitude)))
  # store the prior resolution into a global variable
  res_p<<- res_p_temp

  # resolution (grid size step) associated with the survey
  # (in the case of stratified data, prior and survey resolution are the same)
  # store the prior resolution into a global variable
  res_s<<- res_p

  # minimum value of the resolution (grid size step)
  res_min_temp<- res_p
  # store the minimum resolution into a global variable
  res_min<<- res_min_temp


  #--
  #
  if(DEBUG_MODE==TRUE) {
    cat(paste("\n> DEBUG -\nresolutions (grid size steps):", res_p, "(prior),", res_s, "(survey),", res_min, "(minimum value)\n"))
  }
  #--

  #
  # note that here the term `resolution' refers to the `grid size step'
  # (or `grid size');
  # a low value of the `resolution' corresponds to a high resolution grid;
  # a high value of the `resolution' corresponds to a low resolution grid;
  #


  range_lon<- range(data_p_m$Longitude)[2]
  range_lat<- range(data_p_m$Latitude)[2]

  # search for holes in the prior data (the holes indicate `islands')
  ds_h<- searchHoles(data_p_m, range_lon, range_lat)
  data_p_m<- ds_h$data_p_m
  # grid point IDs associated with the holes within
  # the prior area
  gpoint_id_h<- ds_h$gpoint_id_h

  grid_specs<- c((round((range_lon)/res_min)*2),
                 (round((range_lat)/res_min)*2),
                 res_min/2)

  #--
  #
  if(DEBUG_MODE==TRUE) {
    cat(paste("\n> DEBUG -\ngrid specifications:",
               grid_specs[1], grid_specs[2], grid_specs[3], "\n"))
  }
  #--


  # redefine the data structure by making space for the IDs
  data_p_m<- cbind(data_p_m, VirtID=vector(length=dim(data_p_m)[1]))
  # evaluate the grid locations (IDs of the grid points)
  # associated with the prior data
  for(i in 1:dim(data_p_m)[1]){
    data_p_m$VirtID[i] = getGridPointID(data_p_m$Longitude[i], data_p_m$Latitude[i], grid_specs)
  }

  # evaluate the grid locations (IDs of the grid points)
  # associated with the holes within the prior area
  if(is.null(gpoint_id_h)==F){
    gpoint_id_h<- data_p_m$VirtID[which(is.na(data_p_m$PredictedDensity))]
  }


  # note that, in the case of `stratified' data, a data set for
  # the survey has to be generated;
  #
  # these are the steps ...
  # . consider the survey area
  #   (locations of the boundaries after coord. transformation);
  # . select the grid points that correspond to the survey area;
  # . use these grid points to generate the survey data from the prior data;
  #
  data_b_m<- data_b
  data_b_m$Longitude<- data_b$Longitude + k_lon
  data_b_m$Latitude<- data_b$Latitude + k_lat
  data_p_m_sl<- data_p_m[inout(cbind(data_p_m$Latitude,data_p_m$Longitude),
                               cbind(data_b_m$Latitude[data_b_m$Block==1],
                                     data_b_m$Longitude[data_b_m$Block==1]),
                               bound=T),]
  data_s<- data_p_m_sl[,1:7]
  data_s$PredictedDensity<- data_s_st$PredictedDensity[1]
  # note that it is necessary to apply a correction factor to the Quality Densities
  data_s$QualityDensity<- data_s_st$QualityDensity[1] * sqrt(dim(data_s)[1])
  data_s<- cbind(data_s, VirtID=data_p_m_sl$VirtID)


  # in the case of two or more stratified areas ...
  if(nrow(data_s_st)>1){
    for (i in 2:nrow(data_s_st)){
      data_p_m_sl_temp<- data_p_m[ inout(cbind(data_p_m$Latitude,data_p_m$Longitude),
                                         cbind(data_b_m$Latitude[data_b_m$Block==i],
                                               data_b_m$Longitude[data_b_m$Block==i]),
                                         bound=T),]
      data_s_temp<- data_p_m_sl_temp[,1:7]
      data_s_temp$PredictedDensity<- data_s_st$PredictedDensity[i]
      # if the prior resolution is different from the survey resolution, it is
      # necessary to apply a correction factor to the Quality Densities
      data_s_temp$QualityDensity<- data_s_st$QualityDensity[i] * sqrt(dim(data_s_temp)[1])
      data_s_temp = cbind(data_s_temp, VirtID = data_p_m_sl_temp$VirtID)
      data_s<- rbind(data_s, data_s_temp)
      data_p_m_sl = rbind(data_p_m_sl, data_p_m_sl_temp)
    }
    # re-order the data
    data_s<- data_s[order(data_s$Longitude),]
    data_s<- data_s[order(data_s$Latitude),]
    data_p_m_sl<- data_p_m_sl[order(data_p_m_sl$Longitude),]
    data_p_m_sl<- data_p_m_sl[order(data_p_m_sl$Latitude),]
  }


  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    quilt.plot(data_p_m$Longitude, data_p_m$Latitude, data_p_m$PredictedDensity,
               nrow=length(unique(data_p_m$Longitude)), ncol=length(unique(data_p_m$Latitude)),
               zlim=c(min(range(na.omit(data_s$PredictedDensity))[1], range(na.omit(data_p_m$PredictedDensity))[1]),
                      max(range(na.omit(data_s$PredictedDensity))[2], range(na.omit(data_p_m$PredictedDensity))[2])),
               main="transformed prior and constructed survey")
    quilt.plot(data_s$Longitude,data_s$Latitude, data_s$PredictedDensity,
               nrow=length(unique(data_s$Longitude)), ncol=length(unique(data_s$Latitude)),
               zlim=c(min(range(na.omit(data_s$PredictedDensity))[1], range(na.omit(data_p_m$PredictedDensity))[1]),
                      max(range(na.omit(data_s$PredictedDensity))[2], range(na.omit(data_p_m$PredictedDensity))[2])),
               add=T)
    points(data_s$Lon, data_s$Lat, pch=20, cex=0.5)
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--


  # define `new' and `prior' data
  data_new<- data_s
  data_prior<- data_p_m_sl
  gpoint_id<- data_s$VirtID


  # evaluate `posterior' by combining `prior' and `new' data;
  # update the `prior' surface with the `new' data;
  data_post<- evaluatePosterior(data_new$PredictedDensity, data_new$QualityDensity,
                                data_prior$PredictedDensity, data_prior$QualityDensity)
  cat("\n++ evaluating `posterior' by combining `prior' and `new' data ... done\n")

  if(any(na.omit(data_post$post_y)<0)) {
    cat("\n> warning - zeros in posterior\n")
  }

  # check for updated densities greater or less than max/min of inputs
  if (max(data_post$post_y) > max(data_p$PredictedDensity, data_s$PredictedDensity) |
      min(data_post$post_y) < min(data_p$PredictedDensity, data_s$PredictedDensity) ) {
    cat("\n> warning - updated data outside range of inputs\n")
  }


  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    if(is.null(gpoint_id_h)==T){
      quilt.plot(data_p_m$Longitude, data_p_m$Latitude, data_p_m$PredictedDensity,
                 nrow=length(unique(data_p_m$Longitude)), ncol=length(unique(data_p_m$Latitude)),
                 zlim=c(min(range(na.omit(data_new$PredictedDensity))[1], range(data_p_sl$PredictedDensity)[1]),
                        max(range(na.omit(data_new$PredictedDensity))[2], range(data_p_sl$PredictedDensity)[2])),
                 main='updated surface on prior')
      quilt.plot(data_prior$Longitude, data_prior$Latitude, data_post$post_y,
                 nrow=length(unique(data_new$Longitude)), ncol=length(unique(data_new$Latitude)),
                 zlim=c(min(range(na.omit(data_new$PredictedDensity))[1], range(data_p_sl$PredictedDensity)[1]),
                        max(range(na.omit(data_new$PredictedDensity))[2], range(data_p_sl$PredictedDensity)[2])),
                 add=T)
    }else{
      if(length(which(is.na(data_p_m$PredictedDensity)))>0){
        n_row<- length(unique(data_p_m$Longitude[-which(is.na(data_p_m$PredictedDensity))]))
        n_col<- length(unique(data_p_m$Latitude[-which(is.na(data_p_m$PredictedDensity))]))
      }else{
        n_row=length(unique(data_p_m$Longitude))
        n_col=length(unique(data_p_m$Latitude))
      }
      quilt.plot(data_p_m$Longitude, data_p_m$Latitude, data_p_m$PredictedDensity,
                 nrow=n_row, ncol=n_col,
                 zlim=c(min(range(na.omit(data_new$PredictedDensity))[1], range(data_p_sl$PredictedDensity)[1]),
                        max(range(na.omit(data_new$PredictedDensity))[2], range(data_p_sl$PredictedDensity)[2])),
                 main='updated surface on prior')
      if(length(which(is.na(data_prior$PredictedDensity)))>0){
        n_row<- length(unique(data_prior$Longitude[-which(is.na(data_prior$PredictedDensity))]))
        n_col<- length(unique(data_prior$Latitude[-which(is.na(data_prior$PredictedDensity))]))
      }else{
        n_row=length(unique(data_prior$Longitude))
        n_col=length(unique(data_prior$Latitude))
      }
      quilt.plot(data_prior$Longitude, data_prior$Latitude, data_post$post_y,
                 nrow= n_row, ncol= n_col,
                 zlim=c(min(range(na.omit(data_new$PredictedDensity))[1], range(data_p_sl$PredictedDensity)[1]),
                        max(range(na.omit(data_new$PredictedDensity))[2], range(data_p_sl$PredictedDensity)[2])),
                 add=T)
    }

    points(data_new$Longitude, data_new$Latitude)
    points(data_new$Longitude[which(is.na(data_post$post_y))],
           data_new$Latitude[which(is.na(data_post$post_y))], pch=20)
    cat("\n> DEBUG - plotting data ... done\n")

  }

  #--


  # updated densities (posterior)
  densities_u<- data.frame(data_prior$Latitude, data_prior$Longitude,
                           data_post$post_y, data_post$post_cv)
  names(densities_u)<- c('Latitude', 'Longitude', 'PredictedDensity', 'QualityDensity')

  # coordinate transformation; convert back from virtual grid
  densities_u$Longitude<- densities_u$Longitude - k_lon
  densities_u$Latitude<- densities_u$Latitude - k_lat
  # convert back to degrees and if at 180, make appropriate -/+ split
  if(range(densities_u$Longitude)[2]>180){
    for(i in 1:length(densities_u$Longitude)){
      if(densities_u$Longitude[i]>180) densities_u$Longitude[i]<- (densities_u$Longitude[i] - 360)
    }
  }

  #
  # output processing
  #
  # output:
  #  . updated densities
  write.csv(densities_u, file="densities_u.csv", row.names=FALSE)
  cat("\n++ writing output (updated densities) ... done\n")



  #--


  data_full<- data_p_m

  # selecting the bandwidth
  data_bw<- selectBandwidth(grid_specs, k_lon, k_lat,
                            range_lon, range_lat,
                            data_full, data_prior, data_new, data_post,
                            data_s, data_p_m,
                            gpoint_id, gpoint_id_h, DEBUG_MODE=DEBUG_MODE)
  cat("\n++ selecting the bandwidth ... done\n")



  #--

  gpoints<- matrix(as.matrix(data_bw[,1:2]), ncol=2)
  # data to be smoothed
  data_y<- data_bw$PredictedDensity
  data_h<- cbind(data_bw$h1, data_bw$h2)
  data_id<- data_bw$id

  # grid point IDs associated with the holes within
  # the area to be smoothed
  # (note that `gpoint_id_h_sm' is a subset of `gpoint_id_h')
  gpoint_id_h_sm<- data_id[which(is.na(data_y))]
  if(length(gpoint_id_h_sm)<=1){
    gpoint_id_h_sm<- NULL
  }

  # doing the kernel smoothing
  data_sm_y<- doKernelSmoothing(gpoints, data_y, data_h, data_id, gpoint_id_h_sm)
  cat("\n++ doing the kernel smoothing ... done\n")

  # check for densities produced less than zero
  if(any(na.omit(data_sm_y) < 0)){
    cat("\n> warning - smoothed densities less than zero")
    for (i in 1:length(data_sm_y)){
      if(is.na(data_sm_y[i]) == F){
        if(data_sm_y[i] < 0){
          data_sm_y[i] <- 10^(-15)
        }
      }
    }
  }


  #--
  # -- plots --
  #
  if(is.null(gpoint_id_h)==T){
    quilt.plot(data_p_m$Longitude, data_p_m$Latitude, data_p_m$PredictedDensity,
               nrow=length(unique(data_p_m$Lon)), ncol=length(unique(data_p_m$Lat)),
               zlim=c(min(range(na.omit(data_p_m$PredictedDensity))[1], range(na.omit(data_sm_y))[1]),
                      max(range(na.omit(data_p_m$PredictedDensity))[2], range(na.omit(data_sm_y))[2])),
               xlim=range(data_p_m$Longitude), ylim=range(data_p_m$Latitude),
               main="smoothed surface on prior")
    quilt.plot(data_bw$Longitude, data_bw$Latitude, data_sm_y,
               nrow=length(unique(data_bw$Lon)), ncol=length(unique(data_bw$Lat)),
               zlim=c(min(range(na.omit(data_p_m$PredictedDensity))[1], range(na.omit(data_sm_y))[1]),
                      max(range(na.omit(data_p_m$PredictedDensity))[2], range(na.omit(data_sm_y))[2])),
               add=T)
  }else{
    if(length(which(is.na(data_bw$PredictedDensity)))>0){
      n_row_b<- length(unique(data_bw$Longitude[-which(is.na(data_bw$PredictedDensity))]))
      n_col_b<- length(unique(data_bw$Latitude[-which(is.na(data_bw$PredictedDensity))]))
    }else{
      n_row_b=length(unique(data_bw$Longitude))
      n_col_b=length(unique(data_bw$Latitude))
    }
    quilt.plot(data_bw$Longitude, data_bw$Latitude, data_sm_y,
               nrow=n_row_b, ncol=n_col_b,
               main="smoothed surface (survey area and buffer edge)")
    if(length(which(is.na(data_p_m$PredictedDensity)))>0){
      n_row<- length(unique(data_p_m$Longitude[-which(is.na(data_p_m$PredictedDensity))]))
      n_col<- length(unique(data_p_m$Latitude[-which(is.na(data_p_m$PredictedDensity))]))
    }else{
      n_row=length(unique(data_p_m$Longitude))
      n_col=length(unique(data_p_m$Latitude))
    }
    quilt.plot(data_p_m$Longitude, data_p_m$Latitude, data_p_m$PredictedDensity,
               nrow=n_row, ncol=n_col,
               zlim=c(min(range(na.omit(data_p_m$PredictedDensity))[1], range(na.omit(data_sm_y))[1]),
                      max(range(na.omit(data_p_m$PredictedDensity))[2], range(na.omit(data_sm_y))[2])),
               xlim=c(1,range_lon), ylim=c(1,range_lat),
               main="smoothed surface on prior")
    quilt.plot(data_bw$Lon, data_bw$Lat, data_sm_y,
               nrow=n_row_b, ncol=n_col_b,
               zlim=c(min(range(na.omit(data_p_m$PredictedDensity))[1], range(na.omit(data_sm_y))[1]),
                      max(range(na.omit(data_p_m$PredictedDensity))[2], range(na.omit(data_sm_y))[2])),
               add=T)
  }
  points(cbind(data_b$Longitude+k_lon, data_b$Latitude+k_lat))
  cat("\n++ plotting data ... done\n")

  #--


  data_sm<- data_bw[,1:4]
  # smoothed data
  data_sm$PredictedDensity<- data_sm_y
  # coordinate transformation; convert back from virtual grid
  data_sm$Longitude<- data_sm$Longitude - k_lon
  data_sm$Latitude<- data_sm$Latitude - k_lat
  # convert back to degrees and if at 180, make appropriate -/+ split
  if(range(data_sm$Longitude)[2]>180){
    for(i in 1:length(data_sm$Longitude)){
      if(data_sm$Longitude[i]>180) data_sm$Longitude[i]<- (data_sm$Longitude[i] - 360)
    }
  }


  densities_sm<- data.frame(data_sm$Latitude, data_sm$Longitude,
                            data_sm$PredictedDensity, data_sm$QualityDensity)
  names(densities_sm)<- c('Latitude', 'Longitude', 'PredictedDensity', 'QualityDensity')

  #
  # output processing
  #
  # output:
  #  . smoothed densities
  write.csv(densities_sm, file="densities_sm.csv", row.names=FALSE)
  cat("\n++ writing output (smoothed densities) ... done\n")


}

#' Functions for bandwidth selection
#'
#'
#'

#
# CVS - $Id: bandwidth_selection.R,v 1.2 2010/02/24 13:41:09 lmilazzo Exp $
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
# extrapolateBandwidth()
# bandwidths()
# getRowAuto()
# getColAuto()
# selectDataLowRes()
# getDataRegularDomain()
# getBufferData()
# getColumnIndex()
# selectBandwidth()
#
#--





#--


#
# this function extrapolates the bandwidth
#
# inputs:
#   data_bw  = data associated with the bandwith
#   d1       =
#   d2       =
#   buf_dat  =
#   rowauto  =
#   colauto  =
#   edge_res =
#   t        =
#   b        =
#   r        =
#   l        =
#
extrapolateBandwidth<- function(data_bw, d1, d2, buf_dat, rowauto, colauto,
                                edge_res, t, b, r, l,DEBUG_MODE=FALSE){

  #
  # changing h1 in buffer zone
  #
  # top and bottom continuation, declining to resolution
  if(b>0){
    data_bw$h1[1:(d1*b)]<- rep(seq(edge_res, rowauto$row_band[1], length=b+1)[1:b],
                               each=d1)
  }
  if(t>0){
    data_bw$h1[(dim(data_bw)[1]-(d1*t)+1):(dim(data_bw)[1])]<- rep(seq(rowauto$row_band[length(rowauto$row_band)],
                                                                       edge_res,
                                                                       length=t+1)[1:t],
                                                                   each=d1)
  }
  # left hand side
  for(i in 1:(d2-b-t)){
    if(edge_res==data_bw$h1[buf_dat$zone==2][(colauto$index[i,1])]){
      data_bw$h1[((d1*(b+(i-1)))+1):((d1*(b+(i-1)))+l)]<- rep(edge_res, length=l+1)[1:l]
    }
    else{
      data_bw$h1[((d1*(b+(i-1)))+1):((d1*(b+(i-1)))+l)]<- seq(edge_res,
                                                              data_bw$h1[buf_dat$zone==2][(colauto$index[i,1])],
                                                              length=l+1)[1:l]
    }
  }
  # right hand side
  for(i in 1:(d2-b-t)){
    if(edge_res==data_bw$h1[buf_dat$zone==2][(colauto$index[i,1])]){
      data_bw$h1[((d1*(b+(i-1)))+d1-r+1):((d1*(b+(i-1)))+d1)]<- rep(edge_res,length=r+1)[2:(r+1)]
    }
    else{
      data_bw$h1[((d1*(b+(i-1)))+d1-r+1):((d1*(b+(i-1)))+d1)]<- seq(data_bw$h1[buf_dat$zone==2][(colauto$index[i,1])],
                                                                    edge_res,
                                                                    length=r+1)[2:(r+1)]
    }
  }

  #
  # changing h2 in buffer zone
  #
  # left
  for(i in 1:(d2)){
    data_bw$h2[((d1*(i-1))+1):((d1*(i-1))+l)]<- seq(edge_res,colauto$col_band[1],length=l+1)[1:l]}
  # right
  for(i in 1:(d2)){
    data_bw$h2[(d1*(i)-r+1):(d1*(i))]<- seq(colauto$col_band[length(colauto$col_band)],
                                            edge_res,
                                            length=r+1)[2:(r+1)]
  }
  # bottom
  if(b>0){
    if(b>1){
      bottom<- (l+1):(d1-r)
      for(i in 2:b){
        temp<- (d1*(i-1)+l+1):(d1*i-(r))
        bottom<- rbind(bottom, temp)
      }
      for(i in 1:(d1-r-l)){
        data_bw$h2[bottom[,i]]<- seq(edge_res, data_bw$h2[buf_dat$zone==2][i], length=b+1)[1:b]
      }
    }
    if(b==1){
      bottom<- (l+1):(d1-r)
      for(i in 1:d1){
        data_bw$h2[bottom[i]]<- seq(edge_res, data_bw$h2[buf_dat$zone==2][i], length=b+1)[1:b]
      }
    }
  }
  # top
  if(t>0){
    if(t>1){
      top<- (((d2-t)*d1)+l+1):(((d2-t)*d1)+(d1-r))
      for(i in 2:t){
        temp<- (((d2-t+(i-1))*d1)+l+1):(((d2-t+(i-1))*d1)+(d1-r))
        top<- rbind(top, temp)
      }
      for(i in 1:(d1-r-l)){
        data_bw$h2[top[,i]]<- seq(data_bw$h2[buf_dat$zone==2][i], edge_res, length=t+1)[2:(t+1)]
      }
    }
    if(t==1){
      top<- (((d2-t)*d1)+l+1):(((d2-t)*d1)+(d1-r))
        for(i in 1:d1){
          data_bw$h2[top[i]]<- seq(data_bw$h2[buf_dat$zone==2][i], edge_res, length=t+1)[2:(t+1)]
        }
    }
  }


  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    quilt.plot(data_bw$Longitude, data_bw$Latitude,data_bw$h1, ncol=d2, nrow=d1, main='h1 in whole area')
    quilt.plot(data_bw$Longitude, data_bw$Latitude,data_bw$h2, ncol=d2, nrow=d1, main='h2 in whole area')
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--


  return(data_bw)

}



#
# this function evaluates the bandwidths;
# calculating bandwidths from autocorrelation measures calc in rowauto and colauto
#
# inputs:
#   row.band  = vector containing bandwidths for a seletion of rows (up to 10)
#   row.ind   = vector containing index of rows used to find row.band
#   col.band  = vector containing bandwidths for a seletion of columns (up to 10)
#   index     = (total number of rows x number of columns selected) index of points in each column
#   ind       = vector of id of columns used to get col.band
#   lengths   = number of data points in survey
#   number.rowcol = (1 x 2) number rows, number columns  in survey
#
bandwidths<-function(row.band, row.ind, col.band, index, ind, lengths, number.rowcol) {

  number.row<- number.rowcol[1]
  number.col<- number.rowcol[2]

  # fill in rows
  m<- row.ind[2]-row.ind[1]     #gap between rows used
  temp.h1<- rep(0,number.row)
  temp.h1[1:row.ind[1]]<- row.band[1]

  # fill in rows between row indexes in row.ind
  for (i in 1:(length(row.band)-1)){
      m<- round((row.ind[i+1]-row.ind[i]))
      temp.h1[(row.ind[i]):(row.ind[i+1])]<- seq(row.band[i], row.band[i+1], length=m+1)
  }
  temp.h1[row.ind[length(row.band)]:number.row]<- row.band[length(row.band)]
  h1<- rep(temp.h1, each=(number.col))

  # fill in columns
  h2<- rep(0,lengths)
  mc<- index[1,2]-index[1,1]

  temp.h2<- rep(0,number.col)
  temp.h2[1:index[1,1]]<- col.band[1]
  for (i in 1:(dim(index)[2]-1)){
    mc<- index[1,i+1] - index[1,i]
    if(col.band[i]==col.band[i+1]){
      temp.h2[ind[i]:ind[i+1]]<- rep(col.band[i], length=mc+1)
    }
    else{
      temp.h2[ind[i]:ind[i+1]]<- seq(col.band[i], col.band[i+1], length=mc+1)
    }
  }

  temp.h2[ind[length(ind)]:number.col]<- col.band[length(ind)]
  h2<- rep(temp.h2, number.row)

  return(list(h1=h1, h2=h2))

}



#
# this function evaluates the autocorelation of n rows of data;
# function for finding the autocorelation of n rows of data
#
# inputs:
#   grid         = N x 2 matrix of coordinates
#   y            = vector of length N of animal densities
#   res          = resolution of the grid - ie gap between points
#   rowAvailable = binary vector of rows in the dataset that have enough points to find autocorrelation (1 indicates availability)
#
getRowAuto<- function(grid, y, res, rowAvailable){

  nr<- length(unique(grid[,2]))
  nc<- length(unique(grid[,1]))
  avail<- (1:nr)*rowAvailable
  if(nr<10){
    n<-nr         # if 10 rows not available use as many as there are
  }
  else{
    n<-10
  }            # else use 10 rows
  gap<- round(length(avail[avail>0])/(n+1))
  if(gap==1){
    row_ind<- sample(sample(avail[avail>0],n)) # if gap of 1 randomly sample available indexes
  }else{
    row_ind<- rep(0,n)
    row_ind[1]<- sample(avail[avail>0][1:gap],1)
    for (i in 2:n){
      row_ind[i]<- row_ind[1]+gap*(i-1)
      if(row_ind[i]>nr){
        # if end up with column greater than allowed sample from avail columns not already used
        row_ind[i]<- sample(avail[-row_ind][avail[-row_ind]>0],1)
      }
    }
  }

  row_ind<- row_ind[order(row_ind)]     # order the rows chosen
  row_start<- (row_ind*nc)-(nc-1)       # calculate id's of points in each row from row_ind
  row_end<- (row_ind*nc)

  # autocorrelation for N lines
  row_band<- vector(length=n)

  for(r in 1:n){
    row_band[r]<- evaluateAutocorrelation(na.omit(y[row_start[r]:row_end[r]]))[1]
  }
  row_band<- (row_band*res)/3
  for(r in 1:n){
    if(row_band[r]<=res){
      row_band[r]<- res
    }
  }
  b<- round((row_band[1]*3)/res)
  t<- round((row_band[n]*3)/res)

  return(list(row_band=row_band, start=row_start, end=row_end,
              ind=row_ind, b=b, t=t))

}


#
# this function evaluates the autocorelation of 10 columns of data;
# function for finding the autocorelation of 10 columns of data
#
# inputs:
#   grid         = N x 2 matrix of coordinates
#   y            = vector of length N of animal densities
#   res          = resolution of the grid - ie gap between points
#   rowAvailable = binary vector of rows in the dataset that have enough points to find autocorrelation (1 indicates availability)
#
getColAuto<- function(grid, y, res, columnAvailable){

  nr<- length(unique(grid[,2]))
  nc<- length(unique(grid[,1]))
  avail<- (1:nc)*columnAvailable
  if(nc<10){
    n<- nc         # if 10 rows not available use as many as there are
  }
  else{
    n<- 10
  }            # else use 10 rows
  gap<- round(length(avail[avail>0])/(n+1))
  if(gap==1){
    col_ind<- sample(sample(avail[avail>0],length(avail[avail>0]))) # if gap of 1 randomly sample available indexes
    n<- length(avail[avail>0])
  }
  else{
    col_ind<- rep(0,n)
    col_ind[1]<- sample(avail[avail>0][1:gap],1)      # find start point
    for(i in 2:n){
      col_ind[i]<- col_ind[1]+(gap*(i-1))             # repeatedly add gap
      if(col_ind[i]>nc){
      # if end up with column greater than allowed sample from avail columns not already used
        col_ind[i]<- sample(avail[-col_ind][avail[-col_ind]>0],1)
      }
    }
  }

  col_ind<- col_ind[order(col_ind)]
  col_index<- getColumnIndex(nr=nr,nc=nc,id=col_ind)

  # autocorelation calculation for 10 lines
  col_band<- rep(0,n)
  for(r in 1:n){
    col_band[r]<- evaluateAutocorrelation(na.omit(y[col_index[,r]]))[1]
    if(col_band[r]<=(res/3)){
      col_band[r]<- res/3
    }
  }

  col_band<- (col_band*res)/3
  for(r in 1:n){
    if(col_band[r]<=res){
      col_band[r]<- res
    }
  }
  l<- round((col_band[1]*3)/res)
  r<- round((col_band[n]*3)/res)

  return(list(col_band=col_band, index=col_index, ind=col_ind,
              number_rowcol=c(nr,nc), l=l, r=r))

}


#
# this function selects the correct points for the buffer area;
#
selectDataLowRes = function(data, no_el, iw,
                                  edge1, edge2, no_gpoints_dx, no_gpoints_hx,
                                  no_gp_brx, no_gp_blx){

  # input validation
  if (!no_gpoints_dx %% 2) {
    cat("\n> error - it must be an odd number\n")
  }
  if (!no_gpoints_hx %% 2) {
    cat("\n> error - it must be an odd number\n")
  }
  if (no_gp_brx %% 2) {
    cat("\n> error - it must be an even number\n")
  }
  if (no_gp_blx %% 2) {
    cat("\n> error - it must be an even number\n")
  }

  # define data structure for selected data
  data_sl = c(0)

  # offset associated with the edge 1
  offset_e1 = no_gpoints_dx
  # offset associated with the edge 2
  offset_e2 = no_gpoints_hx

  # no. of selected grid points within the buffer
  # on the right side, in the X direction
  no_gp_brx_sl = no_gp_brx/2

  isl = 1
  ie1 = 1
  ie2 = 1

  for (i in 1:no_el){

    data_sl[isl] = data[iw]
    isl = isl + 1

    if(data[iw] == edge1[ie1]){
      if((no_gp_blx == 0) & (ie2 > 1) & (ie2 < length(edge2))){
        iw = iw + offset_e1 + offset_e2
        ie1 = ie1 + 2
        ie2 = ie2 + 2
      }
      else{
        iw = iw + offset_e1 - 1
        ie1 = ie1 + 2
      }
    }
    iw = iw + 2


    if(iw>no_el) break

    if (ie2 <= length(edge2)) {
      if(data[iw] == edge2[ie2] & no_gp_blx != 0){
        iw = iw + offset_e2 - 1
        ie2 = ie2 + 2
        for (j in 1:no_gp_brx_sl){
          iw = iw + 2
          data_sl[isl] = data[iw]
          isl = isl + 1
          if (j == no_gp_brx_sl) {
            iw = iw + offset_e1 + 1
            ie1 = ie1 + 2
          }
        }
      }
      if(data[iw] == edge2[ie2] & no_gp_blx == 0){
        iw = iw + offset_e2 - 1
        ie2 = ie2 + 2
        for (j in 1:no_gp_brx_sl){
          iw = iw + 2
          data_sl[isl] = data[iw]
          isl = isl + 1
          if (j == no_gp_brx_sl) {
            iw = iw + offset_e1 + offset_e2 + 2
            ie1 = ie1 + 2
          }
        }
      }

    }

    if(iw>no_el) break

  }

  return(data_sl)

}



#
# in the case in which the data are associated with a irregular domain (an
# irregular shaped area), this function evaluates the data associated with an
# area with a regular shape (square or rectangular) that includes the original data;
#
# inputs:
#   grid_specs  = nc, nr, and resolution (grid size step)
#   data_prior  = data set containing points around which the regular domain is to be created
#   data_post   = output from density updater function containing updated densities and cv scores
#   data_full   = data set containing information outside of survey area for getting buffer densities and cv's
#
getDataRegularDomain<- function(a1, p1, grid_specs,
                                data_prior, data_post, data_full, DEBUG_MODE=FALSE){

  nc = length(a1)  - 1
  nr = length(p1)  - 1
  data = 1: (nc * nr)
  no_el = length(data)

  iw = 1
  no_gpoints_hx = 0
  no_gpoints_dx = nc

  no_gp_brx = floor(nc/2)
  no_gp_blx = floor(nc/2)

  # check left and right buffer are even numbers
  if(!no_gp_brx %% 2 == F){
    no_gp_brx = no_gp_brx - 1
  }
  if(!no_gp_blx %% 2 == F){
    no_gp_blx = no_gp_blx - 1
  }

  edge1 = nc

  for(i in 2:nr){
   edge1 = c(edge1, nc*i)
  }

  # internal edge does not exist so is zero
  edge2 = 0

  # find id's of points within whole square
  data_sl = selectDataLowRes(data, no_el, iw,
                             edge1, edge2, no_gpoints_dx, no_gpoints_hx,
                             no_gp_brx, no_gp_blx)

  # find coordinates
  no_dstrs_2 = length(data_sl)
  ds2 = array(dim=c(no_dstrs_2, 5))
  for(i in 1:no_dstrs_2){
    loc = getCoords(data_sl[i], c(nc,nr,grid_specs[3]))
    ds2[i,1] = loc[1]
    ds2[i,2] = loc[2]
    ds2[i,3] = 0
    ds2[i,4] = NA
    ds2[i,5] = data_sl[i]
  }

  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    plot(ds2[,1], ds2[,2], pch=20, col="blue")
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--

  # transform back to same virtual grid as data
  startcoord = a1[1,]

  ds2[,1] = ds2[,1] + startcoord[1] - 1
  ds2[,2] = ds2[,2] + startcoord[2] - 1

  for(i in 1:no_dstrs_2){
    ds2[i,5] = getGridPointID(ds2[i,1], ds2[i,2], grid_specs)
  }

  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    plot(ds2[,1], ds2[,2], pch=20, col="blue")
    points(data_prior$Longitude, data_prior$Latitude, pch=20)
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--


  # find points that overlap between survey and square area
  overlap = no_dstrs_2 - dim(data_prior)[1]
  counter=1
  for (i in 1:no_dstrs_2){
    for(j in 1:dim(data_prior)[1]){
      if(ds2[i,5]==data_prior$VirtID[j]){
          overlap[counter]<- i
          counter<-counter + 1
      }
    }
  }

  points(ds2[-overlap,1],ds2[-overlap,2],pch=20, col='red')
  sq_id = ds2[-overlap,5]

  # make new dataset of survey and buffer.  update buffer to get response values
  sq_up<- updateDataSet(data_full, data.frame(Longitude=ds2[-overlap,1], Latitude=ds2[-overlap,2]),
                         grid_specs, returnID=T)

  sq_dat<- data.frame(rbind(cbind(sq_up$ds2, rep(1,length=length(sq_id))),
                             cbind(data_prior$Longitude, data_prior$Latitude,
                                   data_post$post_y, data_post$post_cv,
                                   data_prior$VirtID, rep(2,dim(data_prior)[1]))))

  names(sq_dat)<- c('Longitude','Latitude','PredictedDensity','QualityDensity','VirtID','zone')
  sq_dat<- sq_dat[order(sq_dat$VirtID),]

  return(sq_dat)

}


#
# this function evaluates the buffer area;
#
# inputs:
#   poly       = (n x 2) points on the boundary
#   roi_le     = distance of the line from the left edge of the ROI (in number of points)
#   roi_re     = distance of the line from the right edge of the ROI (in number of points)
#   roi_be     = distance of the line from the bottom edge of the ROI (in number of points)
#   roi_te     = distance of the line from the top edge of the ROI (in number of points)
#   grid_specs = nc, nr, resolution
#   data_prior = data set containing points around which buffer is to be created
#   data_post  = output from density updater function containing updated densities and cv scores
#   dat_full  = data set containing information outside of survey area for getting buffer densities and cv's
#   unique_lon, unique_lat = unique values of Lon (and Lat) in survey area
#
getBufferData<- function(poly,roi_le, roi_re, roi_be, roi_te, grid_specs,
                         data_prior, data_post, data_full,
                         unique_lon, unique_lat, DEBUG_MODE=FALSE){

  # find id's of polygon points
  no_dstrs_1 = dim(poly)[1]
  ds1 = array(dim=c(no_dstrs_1,4))
  ds1[,1] = poly[,1]
  ds1[,2] = poly[,2]
  ds1[,3] = 0
  ds1[,4] = NA

  for(i in 1:no_dstrs_1){
    ds1[i,4] = getGridPointID(ds1[i,1], ds1[i,2], grid_specs)
  }

  l_points_id = ds1[1,4]
  for(i in 2:no_dstrs_1){
    id_temp = ds1[i,4]
    l_points_id = append(l_points_id, id_temp)
  }

  if(!roi_be%%2 == F){
    roi_be = roi_be-1
    }

  if(!roi_te%%2 == F){
    roi_te = roi_te-1
    }

  if(!roi_le%%2 == F){
    roi_le = roi_le-1
    }

  if(!roi_re%%2 == F){
    roi_re = roi_re-1
  }

  # evaluating the grid locations (IDs of the grid points) within the buffer
  roi_points_id = evaluateROI_BID(grid_specs[1:2], l_points_id, roi_le, roi_re, roi_be, roi_te)

  for (i in 1:length(roi_points_id)){
    if(roi_points_id[i]<0){
      roi_points_id[i]<- NA
    }
  }
  roi_points_id<- na.omit(roi_points_id)

  startpoint = roi_points_id[1]
  endpoint = max(roi_points_id)

  a = getGridPointID(unique_lon[1],unique_lat[1],grid_specs)
  b = getGridPointID(max(unique_lon),unique_lat[1],grid_specs)

  # no. of grid points in the X direction (hole)
  # (this number must be a odd number)
  no_gpoints_hx = b - a + 1

  # no. of grid points in the X direction (full domain)
  no_gpoints_dx = no_gpoints_hx + roi_le + roi_re


  # no. of grid points within the buffer
  # on the right side, in the X direction
  no_gp_brx = roi_re
  #no_gp_brx = 5
  # no. of grid points within the buffer
  # on the left side, in the X direction
  no_gp_blx = roi_le

  # number of rows and columns in buffer area
  nc = no_gpoints_dx
  nr = (length(unique_lat)*2) - 1 + roi_te + roi_be

  data = c(1:(nc*nr))

  # no. of elements (grid points) to process
  no_el = length(data)

  # working index
  # (it refers to the first grid point ID to store)
  iw = 1

  # grid points IDs associated with the right edge of the full domain
  edge1 = nc
  for(i in 2:nr){
    edge1 = c(edge1, nc*i)
  }

  # grid points IDs associated with the left edge of the hole
  edge2 = (nc * roi_be) + roi_le + 1
  for(i in 2:(nr - roi_be - roi_te)){
    edge2 = c(edge2, edge2[i-1] + nc)
  }

  data_sl = selectDataLowRes(data, no_el, iw,
                             edge1, edge2, no_gpoints_dx, no_gpoints_hx,
                             no_gp_brx, no_gp_blx)

  no_dstrs_2 = length(data_sl)
  ds2 = array(dim=c(no_dstrs_2, 5))

  for(i in 1:no_dstrs_2){
    loc = getCoords(data_sl[i], c(nc,nr,grid_specs[3]))
    ds2[i,1] = loc[1]
    ds2[i,2] = loc[2]
    ds2[i,3] = 0
    ds2[i,4] = NA
    ds2[i,5] = data_sl[i]
  }


  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    plot(ds2[,1],ds2[,2],pch=20, col="blue",
         main="buffer on the buffer virtual grid")
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--


  # transform back to same virtual grid as data
  startcoord = getCoords(startpoint,grid_specs)

  ds2[,1] = ds2[,1] + startcoord[1] - 1
  ds2[,2] = ds2[,2] + startcoord[2] - 1

  for(i in 1:no_dstrs_1){
    ds2[i,5] = getGridPointID(ds2[i,1], ds2[i,2], grid_specs)
  }

  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    plot(ds2[,1], ds2[,2], pch=20, col="blue",
         main="buffer on the survey virtual grid, with survey points")
    points(data_prior$Longitude, data_prior$Latitude)
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--


  buffer_ID<- ds2[,5]

  # make new dataset of survey and buffer.  update buffer to get response values
  buf_up<- updateDataSet(data_full, data.frame(Longitude=ds2[,1], Latitude=ds2[,2]),
                         grid_specs, returnID=T)

  buf_dat<- data.frame(rbind(cbind(buf_up$ds2, rep(1,length=length(buffer_ID))),
                             cbind(data_prior$Longitude, data_prior$Latitude,
                                   data_post$post_y, data_post$post_cv,
                                   data_prior$VirtID, rep(2,dim(data_prior)[1]))))


  names(buf_dat)<- c('Longitude','Latitude','PredictedDensity','QualityDensity','id','zone')
  buf_dat<- buf_dat[order(buf_dat$id),]

  return(list(buf_dat = buf_dat, roi_re = roi_re, roi_le = roi_le,
              roi_be = roi_be, roi_te = roi_te))

}


#
# this  function finds the index of points in a particular column
#
# inputs:
#   nr = number of rows in dataset
#   nc = number of columns in dataset
#   id = id of columns from 1 to nc for which points need identifying
#
getColumnIndex<- function(nr, nc, id){

  col_index<- matrix(NA, nr, length(id))

  for(i in 1:length(id)){
    for(r in 1:nr){
      col_index[r,i]<- id[i] + (nc*(r-1))
    }
  }

  return(col_index)

}



#
# this function selects the bandwidth
#
# inputs:
#   grid_specs   = nc, nr, and resolution (grid size step)
#   k_lon, k_lat = translations factors
#   range_lon, range_lat = max ranges of long. (and lat.) between prior and survey
#   data_full    = data associated with the low resolution grid
#   data_prior   = prior data
#   data_new     = new data
#   data_post    = posterior data
#   data_s       = survey data (unmodified, before coord. transformation)
#   data_p_m     = prior data  (modified, after coord. transformation)
#   gpoint_id    = grid point IDs (high resolution grid)
#   gpoint_id_h  = grid point IDs (holes)
#
selectBandwidth<-function(grid_specs, k_lon, k_lat,
                          range_lon, range_lat,
                          data_full, data_prior, data_new, data_post,
                          data_s, data_p_m,
                          gpoint_id, gpoint_id_h,DEBUG_MODE=FALSE){

  # if irregularly shaped, make area square/rectangular to work out bandwidths
  a1<- cbind(seq(min(data_new$Longitude), max(data_new$Longitude), by=res_s),
             rep(min(data_new$Latitude), 1))
  b1<- cbind(seq(min(data_new$Longitude), max(data_new$Longitude), by=res_s),
             rep(max(data_new$Latitude), 1))
  p1<- cbind(rep(min(data_new$Longitude), 1),
             seq(min(data_new$Latitude), max(data_new$Latitude), by= res_s))
  q1<- cbind(rep(max(data_new$Longitude), 1),
             seq(min(data_new$Latitude), max(data_new$Latitude), by= res_s))

  poly<- rbind(a1, q1[2:(dim(q1)[1]-1),], b1, p1[2:(dim(p1)[1]-1),])

  s_col<- length(unique(data_s$Lon))
  s_row<- length(unique(data_s$Lat))

  rowGreater5<- rep(0, s_row)
  unique_lat<- sort(round(unique(data_new$Latitude),4))
  for(i in 1:s_row){      # find number of rows where there are points greater than 5
    if(length(na.omit(data_new$PredictedDensity[data_new$Latitude == unique_lat[i]]))>5){
      rowGreater5[i]<- 1
    }
  }
  colGreater5<- rep(0, s_col)
  unique_lon<- sort(round(unique(data_new$Longitude),4))
  # find number of rows where there are points greater than 5
  col_index<- getColumnIndex(nr=s_row, nc=s_col, id=1:dim(data_new)[1])
  for(i in 1:s_col){
    if(length(na.omit(data_new$PredictedDensity[data_new$Longitude== unique_lon[i]]))>5){
      colGreater5[i]<- 1
    }
  }


  # if there are not enough rows/cols with enough points, extend the data set
  # for the `new' data;
  # consider the data associated with a regular domain (square or rectangular)
  # that includes the irregular shaped area associated iwth the original data

  if(sum(rowGreater5)<10 | sum(colGreater5)<10 | s_col*s_row!=dim(data_new)[1]){

    sq_dat<- getDataRegularDomain(a1, p1, grid_specs,
                                  data_prior, data_post, data_full, DEBUG_MODE = DEBUG_MODE)

    #--
    # -- plots --
    #
    if(DEBUG_MODE==TRUE) {
      plot(data_p_m[,4:3], main="buffer (red) and survey points (black) on virtual grid")
      points(sq_dat[sq_dat$zone==2,1:2],pch=20)
      points(poly,pch='x')
      points(sq_dat[sq_dat$zone==1,1:2],pch=20,col='red')
      cat("\n> DEBUG - plotting data ... done\n")
    }

    #--
    rowauto<- getRowAuto(cbind(sq_dat$Longitude, sq_dat$Latitude),
                               sq_dat$PredictedDensity,
                               min(res_s, res_p),
                               rep(1, s_row))
    roi_te<- rowauto$t
    if(roi_te>14) roi_te<- 14
    roi_be<- rowauto$b
    if(roi_be>14) roi_be<- 14

    colauto<- getColAuto(cbind(sq_dat$Longitude, sq_dat$Latitude),
                               sq_dat$PredictedDensity,
                               min(res_s, res_p),
                               rep(1, s_col))
    roi_re<-colauto$r
    if(roi_re>14) roi_re<- 14
    roi_le<- colauto$l
    if(roi_le>14) roi_le<- 14

    # if survey lies along edge of prior surface then there is no data to have buffer
    if(data_new$Longitude[1] == 1){ roi_le = 0}
    if(data_new$Latitude[1] == 1){ roi_be = 0}

    survey_bands<- bandwidths(rowauto$row_band, rowauto$ind,colauto$col_band,
                              colauto$index, colauto$ind,dim(sq_dat)[1], colauto$number_rowcol)

    if(DEBUG_MODE==TRUE) {
    cat(paste("\n> DEBUG -\n`b', `t', `l', `r' values:",
               roi_be, roi_te, roi_le, roi_re, "\n"))
    }

    buf_dat_temp<- getBufferData(poly, roi_le, roi_re, roi_be, roi_te, grid_specs, sq_dat[,1:5],
                            data.frame(post_y=sq_dat$PredictedDensity, post_cv=sq_dat$QualityDensity),
                            data_full,
                            unique_lon, unique_lat, DEBUG_MODE=DEBUG_MODE)

  }else{

    rowauto<- getRowAuto(cbind(data_new$Longitude, data_new$Latitude),
                         data_post$post_y, min(res_s, res_p), rowGreater5)
    roi_te<- rowauto$t
    if(roi_te>14) roi_te<- 14
    roi_be<- rowauto$b
    if(roi_be>14) roi_be<- 14
    colauto<- getColAuto(cbind(data_new$Longitude,data_new$Latitude),
                         data_post$post_y, min(res_s, res_p),
                         colGreater5)
    roi_re<-colauto$r
    if(roi_re>14) roi_re<- 14
    roi_le<- colauto$l
    if(roi_le>14) roi_le<- 14

    # if survey lies along edge of prior surface then there is no data to have buffer
    if(data_new$Longitude[1] == 1){ roi_le = 0}
    if(data_new$Latitude[1] == 1){ roi_be = 0}

    if(DEBUG_MODE==TRUE) {
    cat(paste("\n> DEBUG -\n`b', `t', `l', `r' values:",
               roi_be, roi_te, roi_le, roi_re, "\n"))
    }

    survey_bands<- bandwidths(rowauto$row_band, rowauto$ind,
                              colauto$col_band, colauto$index, colauto$ind,
                              dim(data_new)[1], colauto$number_rowcol)

    buf_dat_temp<- getBufferData(poly,roi_le, roi_re, roi_be, roi_te, grid_specs, data_prior,
                            data_post,
                            na.omit(data_full),
                            unique_lon, unique_lat,DEBUG_MODE=DEBUG_MODE)

  }

  buf_dat = buf_dat_temp$buf_dat
  roi_le =  buf_dat_temp$roi_le
  roi_be =  buf_dat_temp$roi_be
  roi_te =  buf_dat_temp$roi_te
  roi_re =  buf_dat_temp$roi_re

  buf_dat$Longitude = round(buf_dat$Longitude,4)
  buf_dat$Latitude = round(buf_dat$Latitude,4)

  #--
  #
  if(DEBUG_MODE==TRUE) {
    cat(paste("\n> DEBUG -\n`b', `t', `l', `r' values:",
               roi_be, roi_te, roi_le, roi_re, "\n"))
  }
  #--

  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    plot(data_p_m$Lon, data_p_m$Lat,
         main="buffer (red) and survey points (black) on virtual grid")
    points(buf_dat[buf_dat$zone==2,1:2], pch=20)
    points(buf_dat[buf_dat$zone==1,1:2], pch=20, col='red')
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--

  #--
  # -- plots --
  #
  if(length(which(is.na(buf_dat$PredictedDensity)))>0){
    quilt.plot(buf_dat$Longitude, buf_dat$Latitude, buf_dat$PredictedDensity,
               nrow=length(unique(buf_dat$Longitude[-which(is.na(buf_dat$PredictedDensity))])),
               ncol=length(unique(buf_dat$Latitude[-which(is.na(buf_dat$PredictedDensity))])),
               main="buffer surface containing posterior surface",
               ylim=range(buf_dat$Lat), xlim=range(buf_dat$Lon))
  }else{
    quilt.plot(buf_dat$Longitude, buf_dat$Latitude, buf_dat$PredictedDensity,
               nrow=length(unique(buf_dat$Longitude)), ncol=length(unique(buf_dat$Latitude)),
               main="buffer surface containing posterior surface",
               ylim=range(buf_dat$Lat), xlim=range(buf_dat$Lon))
  }
  points(data_prior$Longitude, data_prior$Latitude)
  points(data_p_m$Longitude, data_p_m$Latitude, pch='x')
  cat("\n++ plotting data ... done\n")

  #--

  # combine coordinates, densities and cv's
  data_bw<- data.frame(buf_dat, h1=rep(0), h2=rep(0))

  data_bw$h1[buf_dat$zone==2]<- survey_bands$h1
  data_bw$h2[buf_dat$zone==2]<- survey_bands$h2


  #--
  # -- plots --
  #
  if(DEBUG_MODE==TRUE) {
    if(length(which(is.na(buf_dat$PredictedDensity)))>0){
      quilt.plot(data_bw$Longitude, data_bw$Latitude, data_bw$PredictedDensity,
                 nrow=length(unique(data_bw$Longitude[-which(is.na(data_bw$PredictedDensity))])),
                 ncol=length(unique(data_bw$Latitude[-which(is.na(data_bw$PredictedDensity))])),
                 main="posterior surface on prior buffer")
    }
    else{
      quilt.plot(data_bw$Longitude,data_bw$Latitude, data_bw$PredictedDensity,
                 nrow=length(unique(data_bw$Longitude)), ncol=length(unique(data_bw$Latitude)),
                 main="posterior surface on prior buffer")
    }
    quilt.plot(data_bw$Longitude,data_bw$Latitude, data_bw$h1,
               nrow=length(unique(data_bw$Longitude)), ncol=length(unique(data_bw$Latitude)),
               main="h1")
    quilt.plot(data_bw$Longitude,data_bw$Latitude, data_bw$h2,
               nrow=length(unique(data_bw$Longitude)), ncol=length(unique(data_bw$Latitude)),
               main="h2")
    cat("\n> DEBUG - plotting data ... done\n")
  }

  #--

  d1<- length(unique(data_bw$Longitude))
  d2<- length(unique(data_bw$Latitude))
  edge_res<- res_s
  if(buf_dat$id[1+roi_le]==gpoint_id[1]) roi_be=0

  # extrapolating bandwidth
  data_bw = extrapolateBandwidth(data_bw, d1, d2, buf_dat, rowauto, colauto,
                                 edge_res, roi_te/2, roi_be/2, roi_re/2, roi_le/2, DEBUG_MODE=DEBUG_MODE)


 return(data_bw)

}

#' 2D Nonparametric Regression (Kernel Smoothing)
#'
#' A function to calculate the regression estimate for each point p and bandwiths h1 and h2. Taken from A.W. Bowman and A. Azzalini, "Applied Smoothing Techniques for Data Analysis", Oxford, Clarendon Press, 1997, pag. 53
#'
#' @param p (x, y) coordinates at which to calculate the regression estimate.
#' @param y Original value of the parameter associated with the point {\code{p}}
#' @param loc locations (matrix of coordinates)
#' @param Y data to be smoothed
#' @param h_values bandwidths
#' @param ind index of points to use to calculate the regression estimate
#'
#' @author Lindesay Scott-Hayward
#'
#' @export
#'
evaluate2DNonParamRegression<- function(p, y, loc, Y, h_values, ind){


  # number of points to use for the calculation
  no_el<- length(ind)

  # define the `design matrix' (or `X matrix')
  x_mx<- matrix(NA, no_el, 3)
  x_mx[,1]<- rep(1, no_el)
  # distances in X direction
  x_mx[,2]<- loc[ind,1] - p[1]
  # distances in Y direction
  x_mx[,3]<- loc[ind,2] - p[2]

  # check for column of zeros in `x_mx' which would lead to singularities;
  # in this case, return the original value
  if(length(which(x_mx[,2]==0))==dim(x_mx)[1] | length(which(x_mx[,3]==0))==dim(x_mx)[1]){
    y
  }
  else{
    # define the `weight matrix' (or `W matrix')
    w_mx<- matrix(0, no_el, no_el)

    # calculate the weights
    w1<- evaluateGaussianKernel(p[1], loc[ind,1], h_values[1])
    w2<- evaluateGaussianKernel(p[2], loc[ind,2], h_values[2])

    diag(w_mx)<- (w1 * w2)

    # make kernel estimate
    solve((t(x_mx) %*% w_mx %*% x_mx)) %*% t(x_mx) %*% w_mx %*% Y[ind]

  }

}


#' Kernel Smoothing
#'
#' @param gpoints grid points (locations)
#' @param data_y data Y to be smoothed
#' @param data_h bandwidths
#' @param data_id IDs associated with the data Y
#' @param hole_id IDs associated with the holes
#'
#' @author Lindesay Scott-Hayward
#'
#' @export
#'
doKernelSmoothing<- function(gpoints, data_y, data_h, data_id, hole_id){

  no_el<- length(data_y)
  no_cols<- length(unique(gpoints[,1]))
  no_rows<- length(unique(gpoints[,2]))

  # define data structure for the smoothed data
  data_sm_y<- vector(length=no_el)

  for(i in 1:no_el){

    if(is.na(data_y[i])){
      data_sm_y[i]<- NA
    }
    else{
      h_values<- data_h[i,]
      # find IDs of points in box around point of interest `i'
      roi_id <- evaluateROI_SFD(c(no_cols, no_rows), i, round(max(h_values)*3/res_min))
      if(is.null(hole_id)==F){
        # make sure any holes in the data are given NA's in corresponding vector of ID's
        for(j in 1:length(roi_id)){
          for(s in 1:length(hole_id)){
            if(data_id[roi_id[j]]==hole_id[s]){
              roi_id[j]<- NA
              break
            }
          }
        }
      }

      # if there are not enough points to smooth just use the original data
      if(length(na.omit(roi_id))<4){
        data_sm_y[i]<- data_y[i]
      }
      else{
        data_sm_y[i]<- evaluate2DNonParamRegression(gpoints[i,], data_y[i],                                               gpoints, data_y,
            h_values, na.omit(roi_id))[1]
      }
    }
  }

  return(data_sm_y)

}

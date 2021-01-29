#..............................................................................................
# R A C O R D  -  Fuctions
#..............................................................................................
# Last update: 2020/04/08

# Colours:
colSource <- rgb(0,0,0,0.3)
colTarget <- rgb(1,0,0,0.3)

# FUNCTIONS FROM OTHER PACKAGES
# kmlShape package (for 'DouglasPeuckerNbPoints' function):
DouglasPeuckerNbPoints <- function(trajx,trajy,nbPoints,spar=NA) {
  missings <- is.na(trajx)|is.na(trajy)
  if(any(missings)){
    trajx <- trajx[!missings]
    trajy <- trajy[!missings]
  }else{}
  if(!is.na(spar)){trajy <- smooth.spline(trajx,trajy,spar=spar)[["y"]]}else{}
  result <- matrix(NA,nbPoints-1,4)
  colnames(result) <- c("firstPoint","lastPoint","middle","maxDist")
  result[1,] <- c(1,length(trajx),findFarestPoint(trajx,trajy))
  for(i in 2:(nbPoints-1)){
    if(max(result[,"maxDist"],na.rm=TRUE)==0){
      warning("[DouglasPeukerNbPoints] the simplified curve perfectly fit with the original with only ",i," points")
      break;
    }else{}
    lineToSplit <- which.max(result[,"maxDist"])
    firstPoint <- result[lineToSplit,1]
    lastPoint <- result[lineToSplit,2]
    middlePoint <- result[lineToSplit,3]
    result[lineToSplit,] <- c(firstPoint,middlePoint,findFarestPoint(trajx[firstPoint:middlePoint],trajy[firstPoint:middlePoint])+c(firstPoint-1,0))
    result[i,] <- c(middlePoint,lastPoint,findFarestPoint(trajx[middlePoint:lastPoint],trajy[middlePoint:lastPoint])+c(middlePoint-1,0))
  }
  x <- c(sort(result[,"firstPoint"]),length(trajx))
  return(data.frame(x=trajx[x],y=trajy[x]))
}

findFarestPoint <- function(trajx,trajy) {
  dmax <- 0
  index <- 1
  end <- length(trajx)
  if(end==2){
    index <- 1
    dmax <- 0
  }else{
    for(i in 2:(end-1)){
      d <- shortestDistanceToLines(Mx=trajx[i],My=trajy[i], Ax=trajx[1],Ay=trajy[1], Bx=trajx[end],By=trajy[end])
      if ( d > dmax ) {
        index <- i
        dmax <- d
      }else{}
    }
  }
  return(c(index=index,dmax=dmax))
}

shortestDistanceToLines <- function(Mx,My,Ax,Ay,Bx,By) {
  aire <- abs((By-Ay)*(Mx-Ax)-(Bx-Ax)*(My-Ay))
  return(  aire / sqrt((Bx-Ax)^2 + (By-Ay)^2))
}


# SUPPORT FUNCTIONS
rot2dM <- function(M, theta) {
  # Rotation of the matrix by an angle
  # Created: 2016/12/01
  # Last Update: 2016/12/01
  
  # Arguments:
  #   M:      matrix of xy-coordinates ('matrix', num, Nx2)
  # Value:
  #   res:      matrix of xy-coordinates ('matrix', num, Nx2)
  
  x <- M[,1]
  y <- M[,2]
  xx <- x*cos(theta)-y*sin(theta)
  yy <- y*cos(theta)+x*sin(theta)
  res <- cbind(xx,yy)
  return(res)
}

deg2rad <- deg.to.rad <- function(deg) {
  radians <- deg*pi/180; return(radians)
}

rmsd <- function(x) {
  # RMSD (Root-mean-square deviation)
  # Created: 2016/12/16
  # Last Update: 2016/12/16
  
  # Arguments:
  #   x:          vector of numbers ('vector', num, N)
  # Value:
  #   rmsd:       root-mean-square-deviation ('vector', num, 1)
  
  rmsd <- sqrt(sum(x^2)/length(x))
  return(rmsd)
}

rmsd.norm <- function(x) {
  # NRMSD=NRMSE (Normalized root-mean-square deviation==Normalized root-mean-square error)
  # Created: 2016/12/16
  # Last Update: 2016/12/16
  # Notes: NRMSD may be calculated as
  #       (i)  RMSD/(max(x)-min(x))   == NRMSD
  #       (ii) RMSD/mean(x)           == NRMSD == CV(RMSD) == coefficient of variation of the RMSD
  # Function calculates (i); for (ii) see 'rmsd.cv' function
  
  # Arguments:
  #   x:          vector of numbers ('vector', num, N)
  # Value:
  #   rmsd.norm:  root-mean-square-deviation ('vector', num, 1)
  
  rmsd.norm <- (sqrt(sum(x^2)/length(x))) / (max(x)-min(x))
  return(rmsd.norm)
}

rmsd.cv <- function(x) {
  # CV(RMSD) Coefficient of variation of the RMSD
  # see 'rmsd.norm' and 'rmsd' functions for more details
  
  rmsd.cv <- (sqrt(sum(x^2)/length(x))) / mean(x)
  return(rmsd.cv)
}

dot.prod <- function(a, b) {
  # Dot/scalar product of two vectors
  # Return multiplied vector (i.e. scalar)
  
  dot.prod <- a%*%b
  return(dot.prod)
}

rechant <- function(M, pts=1000, method="po", plots=F) {
  # Outline resampling
  # Created: 2012
  # Last Update: 2017/02/25
  # Dependencies: 'sp'
  # Note: based on Claude (2008)
  
  # Arguments:
  #   M:          matrix of xy-coordinates ('matrix', num, Nx2)
  #   pts:        number of points ('vector', num, 1)
  #   method:     method used for segmentation ('vector', char, 1)
  #               'po'/'plus.one' to: (i) close outline, (ii) segment to pts+1 segments and (iii) take out the last point
  #               'n' for none
  #               'sp' Claude 2008
  #   plots:      visualisation of process/results ('vector', logical, 1')
  # Value:
  #   res:      matrix of xy-coordinates ('matrix', num, Nx2)
  
  require(sp)
  M.orig=M
  if (method=="n") { M <- M[seq(1,dim(M)[1],length=pts),] }
  else if (method=="plus.one" || method=="po") {
    if (sum(M[1,]==M[dim(M)[1],])<2) { M <- rbind(M,M[1,]) } # check if the outline is closed and if not, close it
    M <- M[seq(1,dim(M)[1],length=pts+1),]  # segment to pts+1 segments
    M <- M[-dim(M)[1],]                     # take out the last point
  }
  else if (method=="sp") {
    Ldig <- Line(M)
    M <- spsample(Ldig, pts, type="regular", offset=c(0,1))@coords
  }
  M <- as.matrix(M)
  if (plots==T) {
    plot(M.orig, asp=1, cex=0.5)
    points(M, col="red", pch=3)
    points(M[1,1], M[1,2], col="blue", cex=1.5)
    points(M[pts,1], M[pts,2], col="green", cex=1.5)
    plot(M.orig, asp=1, type="l", cex=0.5)
    lines(M, col="red", pch=3)
    points(M[1,1], M[1,2], col="blue", cex=1.5)
    points(M[pts,1], M[pts,2], col="green", cex=1.5)
    
  }
  return(M)
}

rechant.equidistant_v2 <- function(M, dista, plots=F) {
  # Create equidistant points along (Poly)line in R
  # Author: mjaskowski
  # Last Update: 2019/06/06
  # source: https://stackoverflow.com/questions/33281319/create-equidistant-points-along-polyline-in-r
  
  # Arguments:
  #   M:          matrix of xy-coordinates ('matrix', num, Nx2)
  #   dista:      distance between points ('vector', num, 1)
  #   plots:      visualisation of process/results ('vector', logical, 1')
  # Value:
  #   eqpoints:   matrix of xy-coordinates ('matrix', num, ?x2)
  
  M_orig = M
  norm_vec <- function(x) sqrt(sum(x^2))
  new_point <- function(p0, p1, di) { # Finds point in distance di from point p0 in direction of point p1
    v = p1 - p0
    u = v / norm_vec(v)
    return (p0 + u * di)
  }
  result = M[1,,drop=FALSE] 
  # for all subsequent points p1, p2 in this data.frame norm_vec(p2 - p1) = M at all times
  equidistantPoints = M[1,,drop=FALSE] 
  M = tail(M, n = -1)
  accDist = 0
  while (length(M) > 0) {
    point = M[1,]
    lastPoint = result[1,]
    dist = norm_vec(point - lastPoint)    
    if ( accDist + dist > dista ) {
      np = new_point(lastPoint, point, dista - accDist)
      equidistantPoints = rbind(np, equidistantPoints) # add np to equidistantPoints
      result = rbind(np, result) # add np to result
      accDist = 0 # reset accDist
    } else {
      #move point from river to result  
      M = tail(M, n = -1)
      result = rbind(point, result)    
      #update accDist  
      accDist = accDist + dist
    }
  }
  # allPoints = result[NROW(result):1,] # reverse result
  # return(list(newPoints = equidistantPoints, allPoints = allPoints))
  if (plots==T) {
    plot(M_orig, asp=1)
    points(M, col='red')
  }
  eqPoints <- equidistantPoints[dim(equidistantPoints)[1]:1,]
  eqPoints <- matrix(eqPoints,dim(eqPoints)[1],dim(eqPoints)[2])
  return(eqPoints)
}

plot.blank <- function() {
  # Blank plot
  
  plot(0, type="n", axes=F, xlab="", ylab="")
}

outline <- function(M) {
  # Create outline
  # Created: 2012

  # Arguments:
  #   M:          outline ('matrix' num, xxx)
  # Value:
  #   res:        matrix containing coordinates('matrix', num, xxx)
  
  Mm <- cbind(M,c(1:(dim(M)[1])))
  y.max <-  Mm[which(Mm[,2]==max(Mm[,2])),]    # matrix of maxs
  y.max <- as.matrix(y.max)
  if (length(y.max)==3) {
    y.max.ind <- y.max[3]
  }
  # if there is two or more points, find the middle one
  else if (length(y.max)>3) {
    x.mean <- mean(y.max[,1])
    y.max.ind <- y.max[which.min(abs(y.max[,1]-x.mean)),3]
  }
  y.max.ind <- as.numeric(y.max.ind)
  OUT <- M[c(y.max.ind:dim(M)[1], 1:(y.max.ind-1)),]
  return(OUT)
}

angle.2vectors_v4 <- angle.2vectors_best <- function(a,b) {
  # Vector between two angles in radians
  # accounts for orientation and can treat higher angles then 0,pi
  # 2019/06/05
  # originally: angle2 <- function(M,N){ atan2(N[2],N[1]) - atan2(M[2],M[1]) }
  
  angle <- atan2(b[2],b[1]) - atan2(a[2],a[1])
  return(angle %% (2*pi))
}

normals2D_points <- function(M, plots=F, norm_len=1) {
  # Calculate normals to an outline (2D points)
  # https://stackoverflow.com/questions/16417891/how-can-i-find-a-normal-vector-of-a-2d-line
  
  # Arguments:
  #   M:          matrix of xy-coordinates ('matrix', num, Nx2)
  #   plots:      visualisation of process/results ('vector', logical, 1)
  #   norm_len:   length of normal vectors ('vector', num, 1)
  # Value:
  #   list containing 2 elements  ('list')
  #             '$normals' Matrix of normals coordinates ('matrix', num, Nx2)
  #             '$normals_inv' Matrix of inverse normals coordinates ('matrix', num, Nx2)
  
  x = M[,1]
  y = M[,2]
  n = dim(M)[1]
  m = dim(M)[2] 
  # calculate all normals
  normals = normals_inv = matrix(0,n,m)
  for (i in 2:(n-1)){
    # MM[i]
    dx = (x[i+1])-(x[i-1])
    dy = (y[i+1])-(y[i-1])
    normals[i,] = c(-dy,dx)
    normals_inv[i,] = c(dy,-dx)
  }
  # calculate the first normal
  dx1 = (x[2])-(x[1])
  dy1 = (y[2])-(y[1])
  normals[1,] = c(-dy1,dx1)
  normals_inv[1,] = c(dy1,-dx1)
  # calculate the last normal
  dxn = (x[n])-(x[n-1])
  dyn = (y[n])-(y[n-1])
  normals[n,] = c(-dyn,dxn)
  normals_inv[n,] = c(dyn,-dxn)
  unit_vector <- function(x) {x / sqrt(sum(x^2))}
  normals = t(apply(normals,1,unit_vector))
  normals_inv = t(apply(normals_inv,1,unit_vector))
  if (plots) {
    plot(x,y,asp=1)
    lines(x,y, col='lightblue')
    for (i in 1:n) {
      lines(rbind(M[i,], (M[i,]+(normals[i,]*norm_len)) ), col='red')
      lines(rbind(M[i,], (M[i,]+(normals_inv[i,]*norm_len)) ), col='blue')
    }
    legend("bottomright", pch=19, legend=c('normals', 'normals_inv'), col=c('red','blue'), cex=0.9)
  }
  return(list(normals=normals, normals_inv=normals_inv))
}

tangents2D_points <- function(M, plots=F, tan_len=1) {
  # Calculate tangent vectors to an outline (2D points)
  # direction goes clockwise
  # Created: 2019/06/05
  
  # Arguments:
  #   M:          matrix of xy-coordinates ('matrix', num, Nx2)
  #   plots:      visualisation of process/results ('vector', logical, 1)
  #   norm_len:   length of tangent vectors ('vector', num, 1)
  # Value:
  #   tangents:   matrix of tangent coordinates ('matrix', num, Nx2)
  
  normals <- normals2D_points(M)$normals
  tangents = cbind(normals[,2],-normals[,1])
  if (plots) {
    plot(M[,1],M[,2],asp=1)
    lines(M[,1],M[,2], col='lightblue')
    for (i in 1:dim(tangents)[1]) {
      lines(rbind(M[i,], (M[i,]+(tangents[i,]*tan_len)) ), col='violet')
    }
    legend("bottomright", pch=19, legend=c('tangents'), col=c('violet'), cex=0.9)
  }
  return(tangents)
}

cumchord <- function(M) {
  # Length of the chord
  # Claude 2008
  
  cumsum(sqrt(apply((M-rbind(M[1,], M[-(dim(M)[1]),]))^2,1,sum)))
}

surf <- function(M) {
  # Calculate surface of the polygon
  # Dependencies: SpatialPolygons {sp}
  
  # Arguments:
  #   M:          matrix of xy-coordinates ('matrix', num, Nx2)
  # Value:
  #   S:          Surface of the polygon ('vector', num, 1)
  
  require(sp)
  if (sum(M[1,]==M[dim(M)[1],])==2) {       # ring must be closed!
    M <- M
  }        
  else {
    M <- rbind(M,M[1,])
  }
  SSp <- SpatialPolygons(list(Polygons(list(Polygon(M)),1)))
  S <- SSp@polygons[[1]]@area
  return(S)
}


# BEZIER
bezier <- function(M, n = dim(M)[1]) {
  # Bezier polynomials
  # Authors: J. Claude (2008)
  # Created: 2012
  # Last Update: 2017/02/25
  # Dependencies: MASS
  # Note: based on Claude (2008)
  
  # Arguments:
  #   M:          outline ('matrix' num, xxx)
  #   n:          Number of estimated Bezier vertices minus one
  # Value:
  #   profil:   list containing 4 elements  ('list')
  #             '$M' Matrix of point coordinates ('matrix', num, Nx2)
  #             '$n' Number of estimated Bezier vertices minus one
  #             '$J' Matrix of Bezier coefficients
  #             '$B' Coordinates of Bezier vertices
  
  require(MASS)
  p <- dim(M)[1]
  if (n != p) {n<-n+1}
  M1 <- M/cumchord(M)[p]
  t1 <- 1-cumchord(M1)
  J <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      J[i, j] <- (factorial(p-1)/(factorial(j-1)* factorial(p-j)))*(((1-t1[i])^(j-1))*t1[i]^(p-j))
    }
  }
  B <- ginv(t(J[,1:n])%*%J[,1:n])%*%(t(J[,1:n]))%*%M
  M<-J[,1:n]%*%B
  list(J=J, B=ginv(t(J[,1:n])%*%J[,1:n])%*%(t(J[,1:n]))%*%M)
}

compare.bez <- function(P.source, P.target, method="rmsd", plots=F) {
  # Comparaison of profiles based on Bezier
  # Function calculates different coeffs/distances between P.source and P.target objects
  # Created: 2014
  # Last Update: 2019/06/11
  # Note: does not perform well on data
  
  # Arguments:
  #   P.source:   'pottery' object of source - list containing 6 elements  ('list')
  #   P.target:   'pottery' object of target - list containing 6 elements  ('list')
  #   method:     minfun:     metric ('vector', char):
  #               'rmsd' for RMSD
  #               'rmsd.cv' for Coefficient of variation of the RMSD
  #               'rmsd.norm' for Normalized RMSD
  #               'ssd' for Summed squared distances
  #               'var' for variance
  #               'cor' for correlation
  #               'cosine_similarity' for cosine similarity
  #   plots:      visualisation of process/results ('vector', logical, 1)
  # Value:
  #   result:     value coeff/distance ('vector' num, 1)
  
  PS <- P.source
  PT <- P.target
  s_min <- which(PT$s==min(PS$s))
  s_max <- which(PT$s==max(PS$s))
  sel <- s_min:s_max
  
  # need resampling because factorial(171) gives Inf
  PS.val <- rechant(PS$extint_seg, 170, "n")
  PS.val <- bezier(PS.val)$B
  PS.val <- c(t(PS.val))
  PT.val <- rechant(PT$extint_seg[sel,], 170, "n")
  PT.val <- bezier(PT.val)$B
  PT.val <- c(t(PT.val))
  if (method=="rmsd") { result <- rmsd(PS.val-PT.val) }
  if (method=="rmsd.cv") { result <- rmsd.cv(PS.val-PT.val) }
  if (method=="rmsd.norm") { result <- rmsd.norm(PS.val-PT.val) }
  if (method=="ssd") { result <- sum(PS.val-PT.val^2) }
  if (method=="var") { result <- var(PS.val-PT.val) }
  if (method=="cor") { result <- cor(PS.val, PT.val) }
  if (method=="cosine_similarity") { result <- dot.prod(PS.val, PT.val) }
  if (plots==T) {
    plot(PS.val, type='l', col=colSource, xlab="Bezier coefficients", ylab="value", lwd=2)
    lines(PT.val, col=colTarget, lwd=2)
  }
  return(result)
}

# DCT (Discrete cosine transform)
dct <- function(M) {
  # Discrete cosine transform
  # Created: 2012/04
  # Last Update: 2012/05/10
  # Note: Originated and inspired by
  #       (i)   CDFT Matlab function (ver 2.7.05) created by Cyril Dommergues (Calcul of DFT for closed outlines)
  #       (ii)  CDFT Matlab function (ver 2.7.4) created by Cyril Dommergues (Calcul of DFT and DCT for open outlines)
  
  # Arguments:
  #   M:          matrix of xy coordinates ('matrix', num, Nx2)
  # Value:
  #   mat:        matrix of xy coefficients ('matrix', num, Nx2)
  
  N <- dim(M)[1]
  st <- complex(real=M[,1],imaginary=M[,2])
  ck <- rep (sqrt(2/N),N)
  ck[1]<-1/sqrt(N)
  Sv <- function(k,N) sum(st*cos(((2*seq(0,N-1)+1)*k*pi)/(2*N)))
  Sk <- mapply(Sv,k=seq(0,N-1),MoreArgs=list(N=N))
  Sk <- ck * Sk
  #    Coefficients are X and Y; out list contains also Amplitudes and Phases
  out <- list(X=Re(Sk), Y=Im(Sk), Amp=Mod(Sk), Phas=Arg(Sk))
  mat <- cbind(out$X, out$Y)
  colnames(mat) <- c("X", "Y")
  return(mat) 
}

idct <- function(Sk, harm, N) {
  # Inverse of the Discrete cosine transform
  # Created: 2012/04
  # Last Update: 2012/05/10
  # Note: Originated and inspired by
  #       (i)   CDFT Matlab function (ver 2.7.05) created by Cyril Dommergues (Calcul of DFT for closed outlines)
  #       (ii)  CDFT Matlab function (ver 2.7.4) created by Cyril Dommergues (Calcul of DFT and DCT for open outlines)
  
  # Arguments:
  #   Sk:         matrix of xy coefficients ('matrix', num, Nx2)
  #   harm:       number of harmonics ('vector', num, 1)
  #   N:          number of points to be reconstructed ('vector', num, 1)
  # Value:
  #   mat:        matrix of ab coefficients ('matrix', num, Nx2)
  
  ck <- rep (sqrt(2/N),harm)
  ck[1] <- 1/sqrt(N)
  Sk <- complex(real=Sk[,1],imaginary=Sk[,2])
  iSv <- function(n,N,H,ck,Sk) sum(ck*Sk*cos(((2*n+1)*seq(0,H-1)*pi)/(2*N)))
  sn <- mapply(iSv,n=seq(0,N-1),MoreArgs=list(N=N,H=harm,ck=ck,Sk=Sk))
  M <- matrix(c(Re(sn),Im(sn)),,2)
  colnames(M) <- c("x", "y")
  return(M)
}

compare.dct <- function(P.source, P.target, nharm, method="rmsd", plots=F) {
  # Comparaison of profiles based on DCT
  # Function calculates different coeffs/distances between P.source and P.target objects
  # Created: 2014
  # Last Update: 2019/06/11
  
  # Arguments:
  #   P.source:   'pottery' object of source - list containing 6 elements  ('list')
  #   P.target:   'pottery' object of target - list containing 6 elements  ('list')
  #   nharm:      number of harmonics ('vector', numeric, 1)
  #   method:     minfun:     metric ('vector', char):
  #               'rmsd' for RMSD
  #               'rmsd.cv' for Coefficient of variation of the RMSD
  #               'rmsd.norm' for Normalized RMSD
  #               'ssd' for Summed squared distances
  #               'var' for variance
  #               'cor' for correlation
  #               'cosine_similarity' for cosine similarity
  #   plots:      visualisation of process/results ('vector', logical, 1)
  # Value:
  #   result:     value coeff/distance ('vector' num, 1)
  
  PS <- P.source
  PT <- P.target
  s_min <- which(PT$s==min(PS$s))
  s_max <- which(PT$s==max(PS$s))
  sel <- s_min:s_max
  PS.val <- dct(PS$extint_seg)
  PS.val <- c(t(PS.val))
  PS.val <- PS.val[1:(nharm*2)]
  PT.val <- dct(PT$extint_seg[sel,])
  PT.val <- c(t(PT.val))
  PT.val <- PT.val[1:(nharm*2)]
  if (method=="rmsd") { result <- rmsd(PS.val-PT.val) }
  if (method=="rmsd.cv") { result <- rmsd.cv(PS.val-PT.val) }
  if (method=="rmsd.norm") { result <- rmsd.norm(PS.val-PT.val) }
  if (method=="ssd") { result <- sum(PS.val-PT.val^2) }
  if (method=="var") { result <- var(PS.val-PT.val) }
  if (method=="cor") { result <- cor(PS.val, PT.val) }
  if (method=="cosine_similarity") { result <- dot.prod(PS.val, PT.val) }
  if (plots==T) {
    plot(PS.val, type='l', col=colSource, xlab="DCT coefficients", ylab="value", lwd=2)
    lines(PT.val, col=colTarget, lwd=2)
  }
  return(result)
}



# RDP (Ramer-Douglas-Peucker algorithm)
compare.rdp <- function(P.source, P.target, method="rmsd", nbPoints=20, plots=F) {
  # Comparaison of profiles based on RDP
  # Function calculates different distances between P.source and P.target objects
  # Created: 2019/07/03
  # Last Update: 2019/12/22
  # Dependencies: 'DouglasPeuckerNbPoints' {kmlShape}
  # More details: see Lucena et al. 2016
  
  # Arguments:
  #   P.source:   'pottery' object of source - list containing 6 elements  ('list')
  #   P.target:   'pottery' object of target - list containing 6 elements  ('list')
  #   method:     minimised function ('vector', char, 1'):
  #               'rmsd' for RMSD
  #               'rmsd.cv' for Coefficient of variation of the RMSD
  #               'rmsd.norm' for Normalized RMSD
  #               'ssd' for Summed squared distances
  #               'ssd.normBySurf' for Summed squared distances normalized by surface
  #               'var' for variance
  #   plots:      visualisation ('vector', logical, 1)
  # Value:
  #   result:     value coeff ('vector' num, 1)
  
  PS <- P.source
  PT <- P.target
  s_min <- which(PT$s==min(PS$s))
  s_max <- which(PT$s==max(PS$s))
  sel <- s_min:s_max
  PS.val <- DouglasPeuckerNbPoints(PS$extint_seg[,1], PS$extint_seg[,2], nbPoints)
  PT.val <- DouglasPeuckerNbPoints(PT$extint_seg[sel,1], PT$extint_seg[sel,2], nbPoints)
  minDist <- sqrt(apply((PS.val-PT.val)^2,1,sum))
  # # !!!! what to do with NA
  # minDist[is.na(minDist)] = median(minDist[!is.na(minDist)])
  if (method=="rmsd") { result <- rmsd(minDist) }
  if (method=="rmsd.cv") { result <- rmsd.cv(minDist) }
  if (method=="rmsd.norm") { result <- rmsd.norm(minDist) }
  if (method=="ssd") { result <- sum(minDist^2) }
  if (method=="ssd.normBySurf") { result <- sum(minDist^2)/surf(PS.val) }
  if (method=="var") { result <- var(minDist) }
  if (plots==T) {
    plot(PS.val, type='l', col=colSource, xlab="r (cm)", ylab="z (cm)", asp=1, lwd=2)
    lines(PT.val, col=colTarget, lwd=2)
  }
  return(result)
}

# RTC (Radius, tangent, curvature)
radius_representation <- function(M, plots=F) {
  # Radius representation
  # Created: 2019/02/03
  # Last Update: 2019/06/12
  # More details: see Adan-Bayewitz et al. 2009; should be also in Gilboa et al. 2004 (NO), RTC et al. 2005 (NO)
  
  # Arguments:
  #   M:          matrix of xy coordinates ('matrix', num, Nx2)
  # Value:
  #   Rs:         Radius representation ('vector', num, N)
  
  Rs <- M[,1]
  if (plots==T) { plot(Rs) }
  return(Rs)
}

tangent_representation <- function(M, plots=F) {
  # The tangent representation
  # Created: 2019/02/03
  # Last Update: 2019/06/12
  # More details: see Saragusti et al. 2005
  # Note: theta(s) = arctan((dy/ds)/(dx/ds))
  
  # Arguments:
  #   M:          matrix of xy coordinates ('matrix', num, Nx2)
  # Value:
  #   Ts:         Tangent representation ('vector', num, N)
  
  p <- dim(M)[1]
  dx <- M[,1]-M[c(p, (1:p-1)),1]
  dy <- M[,2]-M[c(p, (1:p-1)),2]
  ds <- sqrt(dx^2+dy^2)
  Ts <- (atan2((dy/ds),(dx/ds)))
  Ts[1] <- Ts[2]    # "correction"
  if(plots==T){
    plot(Ts, type='l')
  }
  return (Ts)
}

curvature_CurvatureOfTheCurve <- curvature_tangent <- function(M, plots=F) {
  # The curvature (curvature of the curve)
  # Created: 2019/02/03
  # Last Update: 2019/06/12
  # More details: see Gilboa et al. 2004, Saragusti et al. 2005; https://www.math24.net/curvature-radius/
  # Note: k(s) = (dAngle/ds)
  
  # Arguments:
  #   M:          matrix of xy coordinates ('matrix', num, Nx2)
  # Value:
  #   ks:         Curvature representation ('vector', num, N)
  
  tangents <- tangents2D_points(M,plots=F)
  ks <- rep(0,dim(M)[1])
  for (i in 1:(dim(M)[1]-1)) {
    a1 <- angle.2vectors_best(a=c(0,0),b=tangents[i,])
    a2 <- angle.2vectors_best(a=c(0,0),b=tangents[i+1,])
    dA <- a1-a2
    ds <- sqrt(sum((M[i,]-M[i+1,])^2))
    ks[i] <- dA/ds
  }
  if (plots==T) {
    plot(ks, type='l')
  }
  return(ks)
}

compare.rtc <- function(P.source, P.target, w=0, selfun="Rs", method="rmsd.RTC", plots=F) {
  # Comparaison of profiles based on Gilboa et al. 2004 and Karasik et al. 2011
  # Function calculates different distances between P.source and P.target objects
  # expressed as Radius OR Tangent OR Curvature function
  # Created: 2019/02/03
  # Last Update: 2019/06/11
  # More details: Gilboa et al. 2004, Karasik et al. 2011
  
  # Arguments:
  #   P.source:   'pottery' object of source - list containing 10 elements  ('list')
  #   P.target:   'pottery' object of target - list containing 10 elements  ('list')
  #   w:          weights ('vector', logical/numeric, 1/length(s) )
  #               0 for no weights
  #   selfun:     representation to compare ('vector', char, 1)
  #               'Rs': Radius function
  #               'Ts': Tangent function
  #               'Ks': Curvature function
  #   method:    minfun:     minimised function ('vector', char, 1'):
  #               'rmsd.RTC' for rmsd (version from Karasik et al. 2011)
  #               'rms.RTC_PS' for rms calculated on source (version from Karasik et al. 2011)
  #               'rms.RTC_PT' for rms calculated on target (version from Karasik et al. 2011)
  #               'cor.wt' for weighted correlation
  #               'dist.gilbao' for weighted euclidean distance (version from Gilboa et al. 2004)
  #               'cosine_similarity' for cosine similarity
  #   plots:      visualisation of process/results ('vector', logical, 1)
  # Value:
  #   result:     result of calculation ('vector', num, 1)
  
  PS <- P.source
  PT <- P.target
  s_min <- which(PT$s==min(PS$s))
  s_max <- which(PT$s==max(PS$s))
  sel <- s_min:s_max
  # define weights
  if(w==0){
    w <- rep(1,length(PS$s))
  }
  if(selfun=="Rs") {
    PS.val <- PS$Rs
    PT.val <- PT$Rs[sel]
    PT.val_orig <- PT$Rs
  }
  if(selfun=="Ts") {
    PS.val <- PS$Ts
    PT.val <- PT$Ts[sel]
    PT.val_orig <- PT$Ts
  }
  if(selfun=="Ks") {
    PS.val <- PS$Ks
    PT.val <- PT$Ks[sel]
    PT.val_orig <- PT$Ks
  }
  # RMSD of Karasik et al. 2011
  if (method=="rmsd.RTC") {
    L <- sum (w)      # L < - sum ( w * rep(seg_dist,length(PS.val)) )
    result <-  sqrt ( 1/L *  sum (( PS.val-PT.val )^2 * w ) )
  }
  if (method=="rms.RTC_PS") {
    L <- sum (w)      # L < - sum ( w * rep(seg_dist,length(PS.val)) )
    result <-  sqrt ( 1/L *  sum (( PS.val )^2 * w ) )
    plots <- F
  }
  if (method=="rms.RTC_PT") {
    L <- sum (w)      # L < - sum ( w * rep(seg_dist,length(PS.val)) )
    result <-  sqrt ( 1/L *  sum (( PT.val )^2 * w ) )
    plots <- F
  }
  # Correlation
  if (method=="cor") {
    result <- cor(PS.val, PT.val)
  }
  # Weighted correlation
  if (method=="cor.wt") {
    result <- cov.wt(cbind(PS.val, PT.val), wt = w, cor = TRUE)$cor[2,1]
  }
  # Distance (Gilbao et al. 2004)
  if (method=="dist.gilbao") {
    result <- sqrt ( sum( (PS.val - PT.val)^2 * w ) )
  }
  # Cosine similarity
  if (method=="cosine_similarity") {
    result <- dot.prod(PS.val, PT.val)
  }
  if (plots==T) {
    plot(PT$s,PT.val_orig, type='l', col=colTarget, xlab="arc-length (cm)", ylab="")
    lines(PT$s[sel],PT.val, col=colTarget, lwd=2)
    lines(PS$s,PS.val, col=colSource, lwd=2)
    if (selfun=="Rs") {
      title(main=paste("Rs = ", round(result,2), sep=""))
      title(ylab="radius - R(s)")
    }
    if (selfun=="Ts") {
      title(main=paste("Ts = ", round(result,2), sep=""))
      title(ylab="tangent - T(s)")
    }
    if (selfun=="Ks") {
      title(main=paste("Ks = ", round(result,2), sep=""))
      title(ylab="curvature - k(s)")
    }
  }
  return(result)
}

compare.rtc.full <- function(P.source, P.target, w=0, W=c(1,1,1)/3 , plots=T) {
  # Comparaison of profiles based on Gilboa et al. 2004 and Karasik et al. 2011
  # Function calculates different distances between P.source and P.target objects
  # expressed as Radius AND Tangent AND Curvature function
  # Created: 2019/02/03
  # Last Update: 2019/06/11
  
  # Arguments:
  #   P.source:   'pottery' object of source - list containing 10 elements  ('list')
  #   P.target:   'pottery' object of target - list containing 10 elements  ('list')
  #   w:          weights of spline ('vector', logical/numeric, 1/length(s) )
  #               0 for no weights
  #   W:          weights of three functions: c(w_radius, w_tantent, w_curvature)
  #   plots:      visualisation of process/results ('vector', logical, 1)
  # Value:
  #   result:     result of calculation ('vector', num, 1)
  
  if(sum(W)!=1) { warning("Sum of function weights is not equal to one !!!") }
  WR <- W[1]
  WT <- W[2]
  WK <- W[3]
  dR <- compare.rtc(P.source, P.target, w=w, selfun="Rs", method="rmsd.RTC", plots=F)
  dT <- compare.rtc(P.source, P.target, w=w, selfun="Ts", method="rmsd.RTC", plots=F)
  dK <- compare.rtc(P.source, P.target, w=w, selfun="Ks", method="rmsd.RTC", plots=F)
  RR <- mean(compare.rtc(P.source, P.target, w=w, selfun="Rs", method="rms.RTC_PS"),
             compare.rtc(P.source, P.target, w=w, selfun="Rs", method="rms.RTC_PT"))
  TT <- mean(compare.rtc(P.source, P.target, w=w, selfun="Ts", method="rms.RTC_PS"),
             compare.rtc(P.source, P.target, w=w, selfun="Ts", method="rms.RTC_PT"))
  KK <- mean(compare.rtc(P.source, P.target, w=w, selfun="Ks", method="rms.RTC_PS"),
             compare.rtc(P.source, P.target, w=w, selfun="Ks", method="rms.RTC_PT"))
  d <- WR/RR*dR + WT/TT*dT + WK/KK*dK
  if(plots==T) {
    layout(matrix(1:4,2,2))
    s_min <- which(P.target$s==min(P.source$s))
    s_max <- which(P.target$s==max(P.source$s))
    sel <- s_min:s_max
    plot(P.target$profil, asp=1, type="n", xlab="r (cm)", ylab="z (cm)", main=paste("d = ", round(d,2)))
    polygon(P.target$profil, col="lightgrey", border = "darkgrey")
    polygon(P.source$profil, col="lightgrey", border = "darkgrey")
    lines(P.target$extint_seg[sel,], col="blue")
    lines(P.source$extint_seg, col="red")
    compare.rtc(P.source, P.target, w=w, selfun="Ts", method="rmsd.RTC", plots=T)
    compare.rtc(P.source, P.target, w=w, selfun="Rs", method="rmsd.RTC", plots=T)
    compare.rtc(P.source, P.target, w=w, selfun="Ks", method="rmsd.RTC", plots=T)
  }
  return(d)
}

# ICP (with optimisation)
profil.transform <- function (profil, method=c("rtsz"), is.rim=F, par=c(0,0,1,0), plots=F, xlim=c(0,0), ylim=c(0,0)) {
  # The function transforms 'profile' object
  # Created: 2015
  # Last Update: 2016/12/02
  # More details:
  # Notes: 'rtsz' transformation (translate_radius, rotate_theta, scale_size, translate_z)
  #        'strz' transformation (s=scale_by_size, t=rotate_by_theta, r_transl_by_r, z_transl_by_z)
  
  # Arguments:
  #   profil:     list containing 3 elements  ('list')
  #               'profil$profil' contains all profile coordinates ('matrix', num, Nx2)
  #               'profil$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #               'profil$intern' contains intern profile coordinates ('matrix', num, Nx2)
  #   method:     ordre of transformations:
  #               r: translation by r (in mm)
  #               t: rotation by theta (in rad)
  #               s: size (in scalar)
  #               z: translation by z (in mm)
  #   par:        4 transformation parameters ('vector'). see 'method' for more details
  #   plots:      visualisation of process/results ('logical')
  # Value:
  #   profil:     list containing 3 elements  ('list')
  #               'profil$profil' contains all profile coordinates ('matrix', num, Nx2)
  #               'profil$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #               'profil$intern' contains intern profile coordinates ('matrix', num, Nx2)
  #               

  if (plots==T) { profil_orig = profil  }
  pro <- profil$profil; extern <- profil$extern; intern <- profil$intern
  r.tr <- par[1]
  theta <- par[2]
  siz <- par[3]
  z.tr <- par[4]
  if (method=="rtsz") {
    pro[,1] <- pro[,1]+r.tr
    extern[,1] <- extern[,1]+r.tr
    intern[,1] <- intern[,1]+r.tr
    pro <- rot2dM(pro,theta)
    extern <- rot2dM(extern,theta)
    intern <- rot2dM(intern,theta)
    pro <- pro*siz
    extern <- extern*siz
    intern <- intern*siz
    pro[,2] <- pro[,2]+z.tr
    extern[,2] <- extern[,2]+z.tr
    intern[,2] <- intern[,2]+z.tr
  }
  if (method=="strz") {
    pro <- pro*siz
    extern <- extern*siz
    intern <- intern*siz
    pro <- rot2dM(pro,theta)
    extern <- rot2dM(extern,theta)
    intern <- rot2dM(intern,theta)
    pro[,1] <- pro[,1]+r.tr
    extern[,1] <- extern[,1]+r.tr
    intern[,1] <- intern[,1]+r.tr
    pro[,2] <- pro[,2]+z.tr
    extern[,2] <- extern[,2]+z.tr
    intern[,2] <- intern[,2]+z.tr
  }
  if (is.rim==T) {
    # max value is on inter
    if ( max(intern[,2]) > max(extern[,2]) ) {
      idx = which.max(intern[,2]):1             # select indexes of intern which should be extern
      extern = rbind(intern[idx,],extern)       # merge them with extern
      intern = intern[-idx,]                    # ...and take them out from intern
    } else {
      idx = which.max(extern[,2]):1             # select indexes of intern which should be extern
      intern = rbind(extern[idx,],intern)       # merge them with extern
      extern = extern[-idx,]                    # ...and take them out from intern
    }
  }
  profil <- list(profil=pro, extern=extern, intern=intern)
  if (plots==T) {
    pp <- rbind(profil$profil,profil_orig$profil)
    if (xlim[1]==0 & xlim[2]==0 ) { xlim=c(-max(abs(pp[,1])),max(abs(pp[,1]))) }
    if (ylim[1]==0 & ylim[2]==0 ) { ylim=c(min(pp[,2]),max(pp[,2])) }
    plot(profil_orig$profil, asp=1, type="n", xlab="r", ylab="z", xlim=xlim, ylim=ylim)
    polygon(profil_orig$profil, col="lightgrey")
    lines(profil_orig$extern, col="blue")
    lines(profil_orig$intern, col="red")
    polygon(profil$profil, col="grey")
    lines(profil$extern, col="blue")
    lines(profil$intern, col="red")
  }
  return(profil)
}

icp.p.all <- function (P.source, P.target, par=c(0,0,1,0), is.rim=F, seg=F, projection=F, EI="EI", minfun=c("rmsd.cv"), method=c("rtsz"), xlim=c(0,0), ylim=c(0,0), plots=F, plots.full=F, plots.add.res=F, main=NULL, cex.main=1, verbose=F) {
  # ICP performed on two 'profile' objects
  # Function calculates distance between P.source and P.target objects expressed as:
  # 'ssd', 'rmsd.cv' or 'rmsd.norm' (see 'minfun' for more information)
  # Created: 2015
  # Last Update: 2016/12/02
  # More details:
  # Notes:  (i) only 'rtsz' transformation available for the moment
  #         (ii) 'icp.p.f' used for optimisation
  
  # Arguments:
  #   P.source:   'pottery' object of source - list containing 3 elements  ('list')
  #   P.target:   'pottery' object of target - list containing 3 elements  ('list')
  #   par:        4 transformation parameters ('vector', num, 4). c(r.tr,theta,siz,z.tr)
  #   projection: calcualte projection along the normal (deprecated)
  #   EI:         sides used for calculation
  #               'E' for extern
  #               'I' for intern
  #               'EI' for both
  #   minfun:     minimised function ('vector', char, 1'):
  #               'rmsd.cv' for Coefficient of variation of the RMSD
  #               'rmsd.norm' for Normalized RMSD
  #               'ssd' for Summed squared distances
  #               'ssd.normBySurf' for Summed squared distances normalized by surface
  #
  #   method:     ordre of transformations performed on P.source ('vector', char, 1'):
  #               r: translation by r (in mm)
  #               t: rotation by theta (in rad)
  #               s: size (in scalar)
  #               z: translation by z (in mm)
  #   plots:      visualisation of process/results ('vector', logical, 1')
  #   plots_full: visualisation of process/results ('vector', logical, 1')
  #   plots_res:  visualisation of results ('vector', logical, 1')
  #   verbose:    print of process/results ('vector', logical, 1')
  # Value:
  #   res:        vector containing 5 elements  ('vector', num, 5)
  #               first four give parameter values, the fifth the residual distance
  
  P.source <- profil.transform(P.source, par=par, method=method, is.rim, plots=F)
  if (plots==T) {
    if (xlim[1]==0 & xlim[2]==0 ) { xlim=c(min(P.target$profil[,1]),0) }
    if (ylim[1]==0 & ylim[2]==0 ) { ylim=c(min(P.target$profil[,2]),0) }
    plot(P.target$profil, asp=1, type="n", xlim=xlim, ylim=ylim, xlab="r", ylab="z", main=main, cex.main=cex.main)
    polygon(P.target$profil, col=rgb(0,0,0,0.2), border=rgb(0,0,0,0.2))
    polygon(P.source$profil, col=rgb(1,0,0,0.2), border=rgb(1,0,0,0.2))
    if (plots.full==F) {
      if (EI=="EI" || EI=="E") { lines(P.target$extern, col=rgb(0,0,0,0.2)); lines(P.source$extern, col=rgb(1,0,0,0.2), lwd=2) }
      if (EI=="EI" || EI=="I") { lines(P.target$intern, col=rgb(0,0,0,0.2)); lines(P.source$intern, col=rgb(1,0,0,0.2), lwd=2) }
    } else { 
      if (EI=="EI" || EI=="E") { lines(P.target$extern, col=rgb(0,0,1,0.4)); lines(P.source$extern, col=rgb(0,0,1,0.4), lwd=2) }
      if (EI=="EI" || EI=="I") { lines(P.target$intern, col=rgb(1,0,0,0.4)); lines(P.source$intern, col=rgb(1,0,0,0.4), lwd=2) }
    }
  }
  # projection along normal
  if (projection==T) {
    if (EI=="EI" || EI=="E") {
      proj.extern <- projectionL2L_points(P.source$extern, P.target$extern)
      minDistExtern <- sqrt(apply((P.source$extern - proj.extern)^2,1,sum))
      if (plots.full==T) { for (i in 1:dim(P.source$extern)[1]){ lines(rbind(P.source$extern[i,], proj.extern[i,]), col="darkgrey") } }
    }
    if (EI=="EI" || EI=="I") {
      proj.intern <- projectionL2L_points(P.source$intern, P.target$intern)
      minDistIntern <- sqrt(apply((P.source$intern - proj.intern)^2,1,sum))
      if (plots.full==T) { for (i in 1:dim(P.source$intern)[1]){ lines(rbind(P.source$intern[i,], proj.intern[i,]), col="darkgrey") } }
    }
  }
  # only distance
  else {
    if (EI=="EI" || EI=="E") {
      minDistExtern <- rep(NA,dim(P.source$extern)[1])
      for (i in 1:dim(P.source$extern)[1]) {
        minDistExtern[i] <- min(sqrt((P.source$extern[i,1]-P.target$extern[,1])^2+(P.source$extern[i,2]-P.target$extern[,2])^2))
        if (plots.full==T) { lines(rbind(P.target$extern[which.min(sqrt((P.source$extern[i,1]-P.target$extern[,1])^2+(P.source$extern[i,2]-P.target$extern[,2])^2)),],P.source$extern[i,]), col="darkgrey") }
      }
    }
    if (EI=="EI" || EI=="I") {
      minDistIntern <- rep(NA,dim(P.source$intern)[1])
      for (i in 1:dim(P.source$intern)[1]) {
        minDistIntern[i] <- min(sqrt((P.source$intern[i,1]-P.target$intern[,1])^2+(P.source$intern[i,2]-P.target$intern[,2])^2))
        if (plots.full==T) { lines(rbind(P.target$intern[which.min(sqrt((P.source$intern[i,1]-P.target$intern[,1])^2+(P.source$intern[i,2]-P.target$intern[,2])^2)),],P.source$intern[i,]), col="darkgrey") }
      }
    }
  }
  if (EI=="E") { minDist <- minDistExtern }
  if (EI=="I") { minDist <- minDistIntern }
  if (EI=="EI") { minDist <- c(minDistExtern,minDistIntern) }
  minDist[is.na(minDist)] = median(minDist[!is.na(minDist)])  # !!!! what to do with NA
  if (minfun=="rmsd") { result <- rmsd(minDist) }
  if (minfun=="rmsd.cv") { result <- rmsd.cv(minDist) }
  if (minfun=="rmsd.norm") { result <- rmsd.norm(minDist) }
  if (minfun=="ssd") { result <- sum(minDist^2) }
  if (minfun=="ssd.normBySurf") { result <- sum(minDist^2)/surf(P.source$profil) }
  if (minfun=="ssd.normByS") { result <- sum(minDist^2)/par[3] }
  if (minfun=="var") { result <- var(minDist) }
  if (plots.add.res==T) {
    title(main=round(result,4))
  }
  if (verbose==T) {
    print(paste("r:", round(par[1],2), "; theta:", round(par[2],2), "; size:", round(par[3],2), "; z:", round(par[4], 2), "; res:", round(result,4),  sep=""))
  }
  res <- c(par[1],par[2],par[3],par[4],result)
  return(res)
}

icp.p.f <- function(P.source, P.target, par, is.rim, projection, EI, minfun, method, xlim=c(0,0), ylim=c(0,0), plots=F, plots.full=F, plots.add.res=F, verbose=F) {
  # Function used in optimization
  # returns only fifth element of the 'icp.p.all' result
  # See 'icp.p.all' for more details
  return(icp.p.all(P.source=P.source, P.target=P.target, par=par, is.rim=is.rim, projection=projection, EI=EI, minfun=minfun, method=method,
                   xlim=c(0,0), ylim=c(0,0), plots=plots, plots.full=F, plots.add.res=F, verbose=plots)[5])
}

hclimbing <- function (par, fn, change, lower, upper, control, type="min",...) {
  # simulated annealing optimisation (Cortez 2014)
  ### hill.R file ###
  # pure hill climbing:
  # par - initial solution
  # fn - evaluation function
  # change - function to generate the next candidate
  # lower - vector with lowest values for each dimension
  # upper - vector with highest values for each dimension
  # control - list with stopping and monitoring method:
  # $maxit - maximum number of iterations
  # $REPORT - frequency of monitoring information
  # type - "min" or "max"
  # ... - extra parameters for FUN
  
  fpar=fn(par,...)
  for(i in 1:control$maxit) {
    par1=change(par,lower,upper)
    fpar1=fn(par1,...)
    if(control$REPORT>0 &&(i==1||i%%control$REPORT==0))
      cat("i:",i,"s:",par,"f:",fpar,"s'",par1,"f:",fpar1,"\n")
    if( (type=="min" && fpar1<fpar) || (type=="max" && fpar1>fpar)) { par=par1;fpar=fpar1 }
  }
  if(control$REPORT>=1) cat("best:",par,"f:",fpar,"\n")
  
  return(list(sol=par,eval=fpar))
}

hchange <- function (par, lower, upper, dist, round=TRUE,...) {
  # simulated annealing optimisation (Cortez 2014)
  # slight random change of vector par:
  # par - initial solution
  # lower - vector with lowest values for each dimension
  # upper - vector with highest values for each dimension
  # dist - random distribution function
  # round - use integer (TRUE) or continuous (FALSE) search
  # ... - extra parameters for dist
  # examples: dist=rnorm, mean=0, sd=1; dist=runif, min=0,max=1D=length(par) # dimension
  step=dist(D,...) # slight step
  if(round) step=round(step)
  par1=par+step
  # return par1 within [lower,upper]:
  return(ifelse(par1<lower,lower,ifelse(par1>upper,upper,par1)))
}



# PROFILE PREPARATION, BREAKING AND STANDARDISATION
profile2d.xy <- profil.xy <- function(xy, r=xy[which.max(xy[,2]),1], plots=F) {
  # Creation of the 'profile' object from outline defined by xy-coordinates
  # Created: 2015
  # Last Update: 2016/10/30

  # Arguments:
  #   xy:       matrix of xy coordinates ('matrix', num, Nx2)
  #   r:        radius ('vector')
  #   plots:    visualisation of process/results ('logical')
  # Value:
  #   profil:   list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)
  
  profil <- outline(xy)
  profil <- profil
  profil[,1] <- profil[,1]-min(profil[,1])
  profil[,2] <- profil[,2]-min(profil[,2])
  ratio <- profil[which.max(profil[,2]),1]/r
  profil <- profil/ratio
  profil[,1] <- -profil[,1]
  profil[,2] <- profil[,2]-max(profil[,2])
  extern <- profil[1:which(profil[,1]==0)[1],]
  intern <- profil[dim(profil)[1]:which(profil[,1]==0)[length(which(profil[,1]==0))],]
  if (plots==T) {
    plot(profil, type="l", asp=1)
    lines(extern, col="blue", lwd=2)
    lines(intern, col="red", lwd=2)
  }
  return (list(profil=profil, extern=extern, intern=intern))
}

profil.break <- function(profil, breaks=6, only_rim=F, min.pts=20, seed=NULL, plots=F) {
  # Breaking the profile into fragments
  # Created: 2016/06/14
  # Last Update: 2017/01/28; 2019/06/05
  # More details:
  
  # Arguments:
  #   profil:   list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)
  #   breaks:   how many fragments
  #   only_rim: if only rim should be kept
  #   min.pts:  minimal nr of points in the fragment
  #   plots:    visualisation of process/results ('logical')
  # Value:
  #   profil.b: list containing 'breaks' number of elements  ('list')
  #             each element corresponds to fragment
  #             each element contains:
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)

  if (breaks*min.pts>dim(profil$profil)[1]) { warning("Too many breaks for not enough outline points") }
  D <- 1
  if (!is.null(seed)) { set.seed(seed) }
  while ( sum(abs(D)<min.pts)>0 ) {
    ind <- c(1,sort(sample(1:dim(profil$extern)[1],size=breaks-1)),dim(profil$extern)[1])
    D <- ind-ind[c(length(ind)-1, (1:length(ind)-1))]
  }
  if(only_rim==T) { breaks=1 }
  profil.b <- vector("list", breaks)
  for (i in 1:breaks) {
    ext <- profil$extern[ind[i]:ind[i+1],]
    int <- profil$intern[ind[i]:ind[i+1],]
    
    profil.b[[i]]$extern <- ext
    profil.b[[i]]$intern <- int
    profil.b[[i]]$profil <- rbind(ext, int[dim(int)[1]:1,])
  }
  if (plots==T) {
    plot(profil$profil, asp=1, type="n", xlab="r", ylab="z")
    polygon(profil$profil, col="lightgrey", border = "darkgrey")
    lines(profil$extern, col="blue")
    lines(profil$intern, col="red")
    for (i in 1:breaks) {
      polygon(profil.b[[i]]$profil, col="lightgrey", border = "darkgrey", lwd=2)
      lines(profil.b[[i]]$extern, col="blue", lwd=2)
      lines(profil.b[[i]]$intern, col="red", lwd=2)
    }
  }
  return(profil.b)
}

profil.rechant <- function(profil, pts) {
  # Resample extern and intern part of the profile on desired number of points
  # Created: 2019/06/14
  # Last Update: 2019/06/14

  # Arguments:
  #   profil:   list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)
  #   pts:      number of points ('vector', num, 1)
  # Value:
  #   profil:   list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, ptsx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, ptsx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, ptsx2)
  
  extern <- rechant(profil$extern, pts=pts, method="n")
  intern <- rechant(profil$intern, pts=pts, method="n")
  profil <- rechant(profil$profil, pts=pts*2)
  return (list(profil=profil, extern=extern, intern=intern))
}

profil.rechant.equidistant <- function(profil, segment_length=0.2, add=F, plots=F) {
  # Resample extern and intern part of the profile by equally spaced segments
  # Created: 2019/06/06
  # Last Update: 2019/06/06
  
  # Arguments:
  #   profil:   list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)
  #   segment_length:   distance between points ('vector', num, 1)
  #   add:      whatever merge the result with profil
  # Value:
  #   profil:   when add==F => list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, ptsx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, ptsx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, ptsx2)
  #   profil:   when add==T => list containing 6 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)
  #             '$extint_seg' contains all profile coordinates ('matrix', num, ptsx2)
  #             '$extern_seg' contains extern profile coordinates ('matrix', num, ptsx2)
  #             '$intern_seg' contains intern profile coordinates ('matrix', num, ptsx2)
  
  extern <- profil$extern
  intern <- profil$intern
  profil <- profil$profil
  # rechant by segment of the 'segment_length' length
  extern_seg <- rechant.equidistant_v2(extern,segment_length)
  intern_seg <- rechant.equidistant_v2(intern,segment_length)
  if (plots) {
    plot(rbind(extern,intern), asp=1)
    points(extern_seg,col="red", pch=19)
    points(intern_seg,col="blue", pch=19)
  }
  extern_seg <- matrix(extern_seg,dim(extern_seg)[1],dim(extern_seg)[2])
  intern_seg <- matrix(intern_seg,dim(intern_seg)[1],dim(intern_seg)[2])
  extint_seg <- rbind(extern_seg[dim(extern_seg)[1]:1,],intern_seg)
  if (add==F) {
    return (list(extern_seg=extern_seg, intern_seg=intern_seg, extint_seg=extint_seg))
  }
  if (add==T) {
    return (list(profil=profil, extern=extern, intern=intern, extern_seg=extern_seg, intern_seg=intern_seg, extint_seg=extint_seg))
  }
}

profil.plot <- function(profil, main=NULL, col.pro="lightgray", col.pro.border="darkgrey", col.ext="blue", col.int="red", add=F) {
  # Standardization of the 'profile' object according to Rim
  # Created: 2017/03/15
  # Last Update: 2017/03/15
  # More details:
  
  # Arguments:
  #   profil:   list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)
  #   main:     plot title ('char')
  #   col.pro:  color of the profile
  #   col.pro.border: color of the profile outline
  #   col.ext:  color of the extern profile outline
  #   col.int:  color of the outern profile outline
  #   add:      add to graphic ('logical')
  # Value:
  #   graphical
  
  if (add==F) {
    plot(profil$profil, asp=1, type="n", xlab="r", ylab="z", main=main, xlim=range(c(profil$profil[,1],profil$profil[,1])), ylim=range(c(profil$profil[,2],profil$profil[,2])))
  }
  polygon(profil$profil, col=col.pro, border=col.pro.border)
  lines(profil$extern, col=col.ext)
  lines(profil$intern, col=col.int)
}

profil.stand2RimRadius <- function(prof, r, plots=F) {
  # Standardization of 'profile' object accordint to a Rim radius
  # Created: 2019/06/12
  # Last Update: 2019/06/12
  # Note: ratio <- abs(prof$profil[which.max(prof$profil[,2]),1]/r)   # version not good for vessel bases !!!!
  
  # Arguments:
  #   prof:     list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)
  #   r:        radius ('vector', num, 1)         
  #   plots:    visualisation of process/results ('logical')
  # Value:
  #   profil:   list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)
  
  temp = which.max(c  (prof$extern[1,2],prof$intern[1,2]))
  if (temp==1) { px <- abs(prof$extern[1,1]) } else { px <- abs(prof$intern[1,1]) }
  ratio <- px/r
  profil <- prof$profil/ratio
  extern <- prof$extern/ratio
  intern <- prof$intern/ratio
  if (plots==T) {
    plot(profil, type="l", asp=1)
    lines(extern, col="blue", lwd=2)
    lines(intern, col="red", lwd=2)
  }
  return (list(profil=profil, extern=extern, intern=intern))
}

profil.stand2Rim <- function(prof, plots=F) {
  # Standardization of the 'profile' object according to Rim
  # Created: 2019/06/12
  # Last Update: 2019/06/19
  # Note:   # ratio <- abs(prof$profil[which.max(prof$profil[,2]),1])   # version not good for bottoms!!!!
  
  # Arguments:
  #   pro:      list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)
  #   plots:    visualisation of process/results ('logical')
  # Value:
  #   profil:   list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains extern profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains intern profile coordinates ('matrix', num, Nx2)
  
  temp = which.max(c  (prof$extern[1,2],prof$intern[1,2]))
  if (temp==1) { px <- abs(prof$extern[1,1]) } else { px <- abs(prof$intern[1,1]) }
  ratio <- px
  profil <- prof$profil/ratio
  extern <- prof$extern/ratio
  intern <- prof$intern/ratio
  if (plots==T) {
    plot(profil, type="l", asp=1)
    lines(extern, col="blue", lwd=2)
    lines(intern, col="red", lwd=2)
  }
  return (list(profil=profil, extern=extern, intern=intern))
}


# SUPPORT FUNCTIONS
funFun <-  function(LIST, selNormalisation, selSampling, Radius, PtsNb, PtsDist) {
  proc = 1/length(LIST)
  withProgress(message = 'Preparing profiles', value = 0, {
    for (i in 1:length(LIST)) {
      if (selNormalisation=="") {
        LIST[[i]] <- LIST[[i]]
      }
      else if (selNormalisation=="Normalisation_rim") {
        LIST[[i]] <- profil.stand2Rim(LIST[[i]])
      }
      else if (selNormalisation=="Normalisation_radius") {
        LIST[[i]] <- profil.stand2RimRadius(LIST[[i]], r=Radius[i])
      }
      if (selSampling=="Sampling_segments") {
        LIST[[i]] <- profil.rechant(LIST[[i]], pts=PtsNb)
      }
      else if(selSampling=="Sampling_equidistant") {
        LIST[[i]] <- profil.rechant.equidistant(LIST[[i]], segment_length = PtsDist, add=T)
      }
      setProgress(proc*i)
    }
  })
  return(LIST)
}

funFunRtc <- function(LIST, PtsDist) {
  for(i in 1:length(LIST)) { LIST[[i]] <- profil.prepare2rtc(LIST[[i]], PtsDist) }
  return(LIST)
}

funFunDct <- function(LIST, PtsDist) {
  for(i in 1:length(LIST)) { LIST[[i]] <- profil.prepare2dct(LIST[[i]], PtsDist) }
  return(LIST)
}

profil.prepare2dct <- profil.prepare2bez <- function(profil, PtsDist) {
  # Used for DCT, RDP and Bezier
  
  extern <- profil$extern
  intern <- profil$intern
  extern_seg <- profil$extern_seg
  intern_seg <- profil$intern_seg
  extint_seg <- profil$extint_seg
  profil <- profil$profil
  s <- 1:dim(extint_seg)[1]
  s <- s-dim(extern_seg)[1]
  s <- s*PtsDist
  return (list(profil=profil, extern=extern, intern=intern, extern_seg=extern_seg, intern_seg=intern_seg, extint_seg=extint_seg, s=s))
}

profil.prepare2rtc <- function(profil, PtsDist) {
  # Used RTC
  
  extern <- profil$extern
  intern <- profil$intern
  extern_seg <- profil$extern_seg
  intern_seg <- profil$intern_seg
  extint_seg <- profil$extint_seg
  profil <- profil$profil
  Rs <- radius_representation(extint_seg)
  Ts <- tangent_representation(extint_seg)
  Ks <- curvature_CurvatureOfTheCurve(extint_seg)
  s <- 1:dim(extint_seg)[1]
  s <- s-dim(extern_seg)[1]
  s <- s*PtsDist
  return (list(profil=profil, extern=extern, intern=intern, extern_seg=extern_seg, intern_seg=intern_seg, extint_seg=extint_seg, Rs=Rs, Ts=Ts, Ks=Ks, s=s))
}

# A quick function to fit a planar sky (should be able to use hyper.fit instead if needed,
# but SVD gives the least-squares solution right away and we don't need uncertainties)
.svdfitplane <- function(x, fitintercept=TRUE)
{
  ndim = dim(x)[[2]]
  stopifnot(!is.null(ndim) && ndim > 1)
  medians = vector(mode="numeric",ndim)
  x[,1] = -x[,1]
  for(i in 1:ndim)
  {
    medians[i] = median(x[,i])
    x[,i] = x[,i] - medians[i]
  }
  if(fitintercept) x[,4] = x[,1]*0 + 1
  ndim = dim(x)[[2]]
  svdfit = svd(x)$v
  # I don't recall what the purpose of this step is. Oops!
  svdfit = svdfit[,ndim-!fitintercept]/svdfit[1,ndim-!fitintercept]
  # Add the previously subtracted medians into the intercept
  svdfit[ndim] = svdfit[ndim] - sum(svdfit[1:(ndim-fitintercept)]*medians)
  # The first element is unity and not actually useful
  fitpar = svdfit[2:ndim]
  # Return the 3D scatter by projecting the orthognal vector to the plane
  # onto the scatter in the first dimension
  fitpar[ndim] = sd(as.matrix(x) %*% svdfit)
  # svdfit[1] is always unity, but just to be explicit...
  fitpar[ndim+1] = fitpar[ndim]*svdfit[1]/sqrt(sum(svdfit[1:(ndim-fitintercept)]^2))
  
  return(fitpar)
}

.kidsDefault <- function(img, band="r")
{
  gain_eff = 0
  if(band == "r")
  {
    gain_eff = 3.1e13
  } else if(band =="g") {
    gain_eff = 1.75e13
  } else if(band == "u") {
    gain_eff = 3.2e12
  } else if(band == "i") {
    gain_eff = 1.3e13
  }
  dimimg = dim(img)
  regions = list()
  for(dimi in 1:2)
  {
    dmin = 1
    dmax = dimimg[dimi]
    if(dimimg[1] > 200)
    {
      dmid = ceiling(dimimg[dimi]/2)
      dmin = dmid-100
      dmax = dmid+99
    }
    regions[[dimi]] = dmin:dmax
  }
  psfmodel=list(sersic=list(
    re=1.3+0.3*(band=="r"),
    nser=0.5, ang=0, axrat=1, mag=0))
  
  return(list(gain_eff=gain_eff,skylevel=0,xregion=regions[[1]],yregion=regions[[2]],
    psfmodel=psfmodel))
}

.kidsPSFmodel <-function(gid, band)
{
  re = 0
  ang = NA
  axrat = 0
  isr = band=="r"
  isg = band=="g"
  if(gid=="79635")
  {
    re=1.274*isr + 1.68*isg
    ang=130*isr + 44*isg
    axrat = 0.95*isr+0.965*isg
  } else if(gid == "77754")
  {
    re = 2.55*isr + 3.06*isg
    ang = 41*isr + 91*isg
    axrat = 0.94*isr + 0.991*isg
  } else if(gid == "238282") {
    re = 3.3*isg + 2.06*isr
    ang = 127.5*isg + 11*isr
    axrat = 0.961*isg + 0.85*isr
  }

  if((re>0) && (is.finite(ang)) && (axrat) > 0) {
    return(list(sersic=list(mag=0, nser=0.5,
    re=re, ang=ang, axrat=axrat)))
  } else {
    return(NULL)
  }
}

processKidsImage <- function(galaxyid=NULL, band="r", img=NULL, sigma=NULL, orig=NULL)
{
  rval = .kidsDefault(img, band)
  if(is.null(galaxyid)) return(rval)
  isr = band=="r"
  isg = band=="g"
  if(galaxyid == "79635")
  {
    xregion = (126:325)+19*isg
    yregion = (137:336)+8*isg
    gain_eff = 3.150171e+13*isr + 1.749256e+13*isg
    yoff = isg*15
    xbounds = c(121,244,360)
    xbound2 = 245
    xpatch1 = xbounds[1]:(xbounds[2]-1)
    ypatch1 = 15:85+yoff
    skypatch1 = input[xpatch1,ypatch1]
    xpatch2 = xpatch1
    ypatch2 = ((430+5*(band=="g")):473)+yoff
    skypatch2 = input[xpatch2,ypatch2]
    xpatch3 = (xbounds[2]+1):xbounds[3]
    ypatch3 = (31+10*isg):(90-7*isg)
    skypatch3 = input[xpatch3,ypatch3]
    xpatch4 = xpatch3
    ypatch4 = 420:(473+16*isg)
    skypatch4 = input[xpatch4,ypatch4]
    skypatches = list(
      list(list(x=xpatch1,y=ypatch1,z=skypatch1),
           list(x=xpatch2,y=ypatch2,z=skypatch2)),
      list(list(x=xpatch3,y=ypatch3,z=skypatch3),
           list(x=xpatch4,y=ypatch4,z=skypatch4))
    )
    skylevels = numeric(length(skypatches))
    patchskylevels = list(length(skypatches))
    tempcons = list()
    layout(rbind(c(2,3,5),c(6,1,4)))
    for(i in 1:length(skypatches))
    {
      x = numeric(0)
      z = numeric(0)
      y = numeric(0)
      for(j in 1:length(skypatches[[i]]))
      {
        x=c(x,as.vector(row(skypatches[[i]][[j]]$z))-1+skypatches[[i]][[j]]$x[1])
        y=c(y,as.vector(col(skypatches[[i]][[j]]$z))-1+skypatches[[i]][[j]]$y[1])
        z=c(z,as.vector(skypatches[[i]][[j]]$z)) 
      }
      tofit = data.frame(z=z, x=x, y=y)
      
      tempcon = 0*input
      #tmp = hyper.fit(tofit)
      fit = .svdfitplane(tofit)
      skymodel = fit[1]*x + fit[2]*y + fit[3]
      xi = 1
      for(j in 1:length(skypatches[[i]]))
      {
        nx = length(skypatches[[i]][[j]]$z)
        meanz = mean(skypatches[[i]][[j]]$z)
        print(paste0("Pre/post means: ",meanz,",",mean(
          skypatches[[i]][[j]]$z-skymodel[xi:(xi+nx-1)])))
        xi = xi+nx
        tempcon[skypatches[[i]][[j]]$x,skypatches[[i]][[j]]$y] = 1
        tempcons[[i]] = tempcon
      }  
      x = xbounds[i]:(xbounds[i+1]-1)
      y = 1:dim(input)[[2]]
      xa=as.vector(row(input[x,y])) + xbounds[i]-1
      ya=as.vector(col(input[x,y]))
      skymodel = fit[1]*xa + fit[2]*ya + fit[3]
      input[x,y] = input[x,y] - skymodel
      xi = 1
      skylevels[i] = gain_eff*var(as.vector(tofit$z - fit[1]*tofit$x - fit[2]*tofit$y - fit[3]))
      patchskylevels[[i]] = numeric(length(skypatches[[i]]))
      for(j in 1:length(skypatches[[i]]))
      {
        print(paste0("ivar[",i,",",j,"]=",1/var(as.vector(skypatches[[i]][[j]]$z))))
        nz = length(skypatches[[i]][[j]]$z)
        zs = xi:(xi+nz-1)
        skypatch = tofit$z[zs]  - (fit[1]*tofit$x[zs] + fit[2]*tofit$y[zs] + fit[3])
        xi = xi + nz
        hist(skypatch,breaks = 100,freq=FALSE,main="",xlab="mgy")
        title(paste0("skypatch[",i,",",j,"]"))
        xh = seq(-2e-11,3e-11,1e-13)
        skylevel = gain_eff*var(as.vector(skypatch))
        lines(xh,dnorm(xh,sd = mean(sigma[x,y])),col="red")
        lines(xh,dnorm(xh,sd = sqrt(skylevel/gain_eff)),col="blue")
        lines(xh,dnorm(xh,sd = sqrt(skylevels[i]/gain_eff)),col="green")
        patchskylevels[[i]][j] = skylevel
      }
    }
    magimage(input,col=cmap)
    for(i in 1:length(skypatches))
    {
      tempcon = tempcons[[i]]
      tempcon=magimage(tempcon,add=T,col=NA)#hsv(s=0,alpha=0.5)
      contour(tempcon,add=T,drawlabels = F,levels=1,col='cyan')
    }
    print(paste0("ivar top left: ",1/var(as.vector(orig[700:830,270:370]))))
    print(paste0("ivar top right: ",1/var(as.vector(orig[1750:1860,475:610]))))
    skylevel = var(as.vector(orig[1750:1860,475:610]))*gain_eff
    
    # Sanity check - a patch of sky not affected by the dither artifacts near the galaxy
    # Variance ~ counts * throughput^2 
    sigma = sqrt((input/throughput_atm + skylevel/throughput_sys)/gain_eff)
    
    # Remap 
    input = input[xregion,yregion] + skylevel
    sigma = sigma[xregion,yregion]
    magimage(input,col=cmap)
    par(mfrow=c(1,1))
  } else {
    if(galaxyid == "77754") {
      skylevel = var(as.vector(input[205:250,280:400]))*rval$gain_eff
      gain_eff = 3.077552e+13*isr + 1.736964e+13*isg
    } else if(galaxyid == "238282") {
      skylevel = var(as.vector(input[21:49,20:50]))*rval$gain_eff
      gain_eff = 3.047264e+13*isr + 1.731326e+13*isg
      rval$xregion = (36:155)-5*isg
      rval$yregion = (54:173)-6*isg
    }
    input = input[rval$xregion,rval$yregion]
    # Not needed since we don't correct the sky level
    #if(exists("gain_eff") && gain_eff > 0) sigma = sqrt(((input-skylevel)/throughput_atm + skylevel/throughput_sys)/gain_eff)
    sigma = sigma[rval$xregion,rval$yregion]
    input = input + skylevel
  }
  if(exists("skylevel")) rval$skylevel = skylevel
  if(exists("gain_eff")) rval$gain_eff = gain_eff
  if(exists("input")) rval$input = input
  if(exists("sigma")) rval$sigma = sigma
  if(exists("xregion")) rval$xregion = xregion
  if(exists("yregion")) rval$yregion = yregion
  psfmodel = .kidsPSFmodel(galaxyid,band)
  if(!is.null(psfmodel)) rval$psfmodel=psfmodel
  return(rval)
}
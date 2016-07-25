# Loads and fits a SAMI galaxy given by the string gid (which should be an integer)
stopifnot(exists("gid") && is.numeric(as.integer(gid)) && exists("band"))

require(RColorBrewer)
#cmap = rev(colorRampPalette(brewer.pal(9,'YlOrBr'))(200))
errcmap = rev(colorRampPalette(brewer.pal(9,'RdYlBu'))(200))
cmap = errcmap

datapath = "inst/extdata/"
source("inst/example/SAMIDataPrep.R")
gpath = paste0("~/raid/sami/kids/G",gid,"/",band,"/")

if(FALSE)
{
  input= readFITS(paste0(datapath,"G",gid,"_r_fitim.fits"))
  sigma = readFITS(paste0(datapath,"G",gid,"_r_sigma.fits"))
  mask = readFITS(paste0(datapath,"G",gid,"_r_mskim.fits"))$imDat
  segim = readFITS(paste0(datapath,"G",gid,"_r_segim.fits"))
  psfim = readFITS(paste0(datapath,"G",gid,"_r_psfim.fits"))$imDat
  psfims = readFITS(paste0(datapath,"G",gid,"_r_psfim.fits"))$imDat
  orig = readFITS(paste0("~/raid/sami/kids/G",gid,"/r/cutim.fits"))$imDat
} else {
  input= readFITS(paste0(gpath,"fitim.fits"))
  sigma = readFITS(paste0(gpath,"sigma.fits"))
  mask = readFITS(paste0(gpath,"M01_mskim.fits"))$imDat
  segim = readFITS(paste0(gpath,"segim.fits"))
  psfim = readFITS(paste0(gpath,"psfim.fits"))$imDat
  psfims = readFITS(paste0(gpath,"psfim.fits"))$imDat
  orig = readFITS(paste0(gpath,"cutim.fits"))$imDat
}

sigmahead = sigma$header
inputhead = input$header
segimhead = segim$header
sigma = sigma$imDat
input = input$imDat
segim = segim$imDat

# KiDS-specific info
# https://www.eso.org/sci/facilities/paranal/instruments/omegacam/inst.html
# See above: Should probably use QE*filter throughput
# http://www.e2v.com/resources/account/download-datasheet/1238
# https://www.eso.org/sci/facilities/paranal/instruments/omegacam/doc/VST-MAN-OCM-23100-3110-2_7_1.pdf
ccdgains = c(2.37,2.52,2.62,2.56,2.56,2.78,2.73,2.37,2.57,2.56,2.56,2.46,2.40,2.32,2.39,2.52,2.40,2.48,2.52,2.44,2.66,2.71,2.67,2.57,2.39,2.59,2.49,2.55,2.48,2.23,2.54,2.39)
gain_inv_e = mean(ccdgains)
if(band == "r") {
  throughput_sys = 0.42
  throughput_atm = 0.385
} else if(band == "g") {
  throughput_sys = 0.46
  throughput_atm = 0.4
}

# Overrides for this particular galaxy
datamod = processKidsImage(gid, band, input, sigma, orig)
datamodnames = names(datamod)
if("input" %in% datamodnames) input = datamod$input
if("sigma" %in% datamodnames) sigma = datamod$sigma
if("xregion" %in% datamodnames && "yregion" %in% datamodnames)
{
  mask = mask[datamod$xregion,datamod$yregion]
  segim = segim[datamod$xregion,datamod$yregion]
}
psfmodel = datamod$psfmodel
skylevel = datamod$skylevel
gain_eff = datamod$gain_eff

writecorrected = TRUE
if(writecorrected) {
  writeFITSim(sigma,file=paste0("~/raid/sami/models/",gid,"/kids_sigma_",
    band,".fits"),header=sigmahead)
  writeFITSim(input-skylevel,file=paste0("~/raid/sami/models/",gid,"/kids_fitim_",
    band,".fits"),header=inputhead)
  writeFITSim(segim==segim[dim(segim)[1]/2,dim(segim)[2]/2],file=paste0("~/raid/sami/models/",
    gid,"/kids_mask_",band,".fits"),header=inputhead)
}

skymin = log10(skylevel/4)
skymax = log10(skylevel*4)

xc = dim(input)[1]/2
yc = dim(input)[2]/2
rc = sqrt(xc^2 + yc^2)

#Very rough model (not meant to look too good yet):
model=list(
  sersic=list(
    xcen= list(xc, xc),
    ycen= list(yc, yc),
    mag= list(16, 16),
    re= list(rc/16, rc/4),
    nser= list(4.0, 1.0),
    ang= list(0.0, 135), #theta/deg: 0= |, 45= \, 90= -, 135= /, 180= |
    axrat= list(1.0, 0.6) #min/maj: 1= o, 0= |
  ),
  magzero=list(0),
  sky=list(bg=skylevel)
)

finesample = 3L
psfarea = pi*(psfmodel$sersic$re)^2*psfmodel$sersic$axrat
psfdim = ceiling(psfmodel$sersic$re*3)
psfdim = as.integer(psfdim + !(psfdim %% 2))
psf = profitMakePointSource(model=psfmodel, image = matrix(0,psfdim,psfdim))
psff = profitMakePointSource(model=psfmodel, image = matrix(0,psfdim*finesample,psfdim*finesample))

# The pure model (no PSF):
magimage(profitMakeModel(model,dim=dim(input)),magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# The original image:
magimage(input,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# The convolved model (with PSF):
modelconv = profitMakeModel(model,psf=psf,dim=dim(input))
magimage(modelconv,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# The convolved model (with fine-sampled PSF):
modelconvfd = profitMakeModel(model,psf=psff,dim=dim(input),finesample=finesample)
magimage(modelconvfd,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# The convolved model (with fine-sampled PSF), prior to downsampling:
modelconvf = profitMakeModel(model,psf=psff,dim=dim(input),finesample=finesample,returnfine = TRUE)
magimage(modelconvf,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# The fine-sampled model, pre-convolution, prior to downsampling:
modelf = profitMakeModel(model,dim=dim(input),finesample=finesample,returnfine = TRUE)
magimage(modelf,magmap=T,stretch='asinh',stretchscale=1/median(abs(input)),col=cmap)

# Difference between fine and non-finesampled model
diff = (modelconvfd$z-modelconv$z)/modelconvfd$z
magimage(diff,col=errcmap,stretch='asinh',stretchscale=1/median(abs(diff)))

# What should we be fitting:

tofit=list(
  sersic=list(
    xcen= list(TRUE,NA), #We will trust that the x and y positions are okay already
    ycen= list(TRUE,NA), #We will trust that the x and y positions are okay already
    mag= list(TRUE,TRUE),
    re= list(TRUE,TRUE),
    nser= list(TRUE,FALSE), #The second sersic is our disk- we will fix this for our first fit
    ang= list(FALSE,TRUE), #The bulge will be fixed to have axrat=1, so no need to fir for the orientation
    axrat= list(FALSE,TRUE) #The bulge has axrat=1 for our first fit
  ),
  magzero=list(FALSE),
  sky=list(bg=TRUE)
)

# What parameters should be fitted in log space:

tolog=list(
  sersic=list(
    xcen= list(F,F),
    ycen= list(F,F),
    mag= list(F,F),
    re= list(T,T), #re is best fit in log space
    nser= list(T,T), #nser is best fit in log space
    ang= list(F,F),
    axrat= list(T,T) #axrat is best fit in log space
  ),
  magzero=list(FALSE),
  sky=list(bg=T)
)

# The priors. If the parameters are to be sampled in log space (above) then the priors will refer to dex not linear standard deviations. Priors should be specified in their unlogged state- the logging is done internally.

priorsd=list(
  sersic=list(
    xcen= list(2.5,2.5),   # should have tight constraints on x and y
    ycen= list(2.5,2.5),   # should have tight constraints on x and y
    mag= list(5.0, 5.0),   # 5 mag SD
    re= list(1.0, 1.0),    # i.e. 1 dex in re is the SD
    nser= list(1.0, 1.0),  # i.e. 1 dex in nser is the SD
    ang= list(30.0, 30.0), # very broad 30 deg ang SD
    axrat= list(1.0, 0.6)  # i.e. 1 dex in axrat is the SD
  ),
  magzero=list(5),
  sky=list(bg = 0.05)
)

priors = as.list(unlist(priorsd))
npriors = length(priors)
for(p in 1:npriors)
{
  priors[[p]] = eval(bquote(function(x){dnorm(x,0,.(priors[[p]]))}))
}
priors = relist(priors,priorsd)

#the hard intervals should also be specified in log space if relevant:

intervals=list(
  sersic=list(
    xcen=list(function(x){interval(x,-Inf,Inf,reflect=F)},function(x){interval(x,-Inf,Inf,reflect=F)}),
    ycen=list(function(x){interval(x,-Inf,Inf,reflect=F)},function(x){interval(x,-Inf,Inf,reflect=F)}),
    mag=list(function(x){interval(x,10,30,reflect=F)},function(x){interval(x,10,30,reflect=F)}),
    re=list(function(x){interval(x,-1,2.5,reflect=F)},function(x){interval(x,-1,2.5,reflect=F)}), # i.e. 1 dex in re is the SD
    nser=list(function(x){interval(x,log10(0.5),1.3,reflect=F)},function(x){interval(x,-0.3,1.3,reflect=F)}), # i.e. 1 dex in nser is the SD
    ang=list(function(x){x = ((x+180) %% 360) - 180},function(x){interval(x,-Inf,Inf,reflect=F)}),
    axrat=list(function(x){interval(x,-2,0,reflect=F)},function(x){interval(x,-2,0,reflect=F)}) # i.e. 1 dex in axrat is the SD
  ),
  magzero=list(function(x){interval(x,-Inf,Inf,reflect=F)}),
  sky=list(bg=list(function(x){interval(x,skymin,skymax,reflect=F)}))
)

#Setup the data structure we need for optimisation:

psf = psff
if(!exists("DataG")) DataG=profitSetupData(image=input,mask=mask,sigma=sigma,segim = segim,psf = psf,model = model, 
  tofit = tofit, tolog=tolog, priors = priors, intervals=intervals,algo.func = "", finesample=finesample, verbose=TRUE) 

DataG$psfarea = psfarea
DataG$gain = gain_eff

# This produces a fairly complex R object, but with all the bits we need for fitting, e.g. (notice the tolog parameteres are now logged):

DataG$init

#These are the parameters we wish to fit for, and we take the initial guesses from the model list we provided before.

#We can test how things currently look (we get an output because we set verbose=TRUE earlier):

rv = profitLikeModel(DataG$init,DataG,makeplots=T,cmap = cmap, errcmap=errcmap)

# Now with covariance estimated from the model + sky + gain:
#profitLikeModel(DataG$init,DataG,estcovar = TRUE, image=TRUE)

best = c(100.21085,99.19477,19.70171,14.50953,0.44658,1.8276456,-0.208384,136.62,-0.2181,-9.7101,109.707)
names(best) = DataG$parm.names
bestf = c(100.212,99.179,19.66,14.513,0.454,1.829,-0.167,136.07,-0.21936,-9.69377)
if(band == "g") bestf = c(100.6264050,99.4205854,20.4668126,14.8644415,0.7693278,1.9200132,0.6111524,135.7312155,-0.2196236,-10.3740616)
if(gid == "77754") bestf = c(1.004505e+02,9.966245e+01,1.845005e+01,1.531247e+01,5.531417e-01,1.573627e+00,-9.399536e-03,1.708520e+02,-2.577917e-01,-9.976029e+00)
if(gid == "238282") bestf = c(59.97907394,59.02058617,18.32457835,16.33379224,0.01615005,1.26221494,0.16354749,32.64407153,-0.20673562,-9.77888680)
names(bestf) = DataG$parm.names
init = bestf

# bestldg
# sersic.xcen1  sersic.ycen1   sersic.mag1   sersic.mag2    sersic.re1    sersic.re2  sersic.nser1   sersic.ang2 sersic.axrat2        sky.bg 
# 100.5814694    99.3961693    20.9342004    14.8734759     0.4865213     1.9156031    -0.3003859   135.7662751    -0.2187757   -10.3727678

rv = profitLikeModel(bestf,DataG,makeplots=T,cmap = cmap, errcmap=errcmap)

dola = FALSE
docma = FALSE
dold = TRUE

if(dola)
{
  DataG$algo.func = "LA"
  LAfit=LaplaceApproximation(profitLikeModel,parm=DataG$init,Data=DataG,Iterations= 1e4,Method='BFGS',CovEst='Identity',sir=FALSE)
  #The best BFGS fit is given by:
  
  #LAfit_best = LAfit$Summary1[,1]
}

# Fit using CMA-ES, an evolutionary algorithm that adapts based on the covariance matrix
# Should be more robust to being trapped in local minima and reasonably fast

if(docma)
{
  require(cmaeshpc)
  # The priors are a bit too broad, perhaps...
  DataG$algo.func = "CMA"
  cma_sigma = unlist(priorsd)[which(unlist(tofit))]/4.0
  cmafit = cmaeshpc(DataG$init, profitLikeModel, Data=DataG, # lowerlims=lowerlim, upperlims=upperlim, lower=lowerlim, upper=upperlim,
    control=list(maxit=2e3,diag.sigma=TRUE,diag.eigen=TRUE,diag.pop=TRUE,diag.value=TRUE,
    fnscale=-1.0,sigma=cma_sigma,maxwalltime=Inf, trace=TRUE, stopfitness = 0, stop.tolx=1e-2*cma_sigma))
  
  # Best CMA fit:
  profitLikeModel(cmafit$par, DataG, makeplots = T)
}

if(dold)
{
  #Now we can try a LaplacesDemon fit:
  DataG$algo.func="LD"
  LDfits = list()
  niters = 5e3
  # This algorithm takes forever and returns hideous posteriors
  #t1 = proc.time()[['elapsed']]
  #LDfits$ADMG=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataG,Iterations=niters,Algorithm="ADMG",
  #  Thinning=1,Specs=list(n=0,Periodicity=2*length(init)))
  t2 = proc.time()[['elapsed']]
  
  # Found very low acceptance rates with this algorithm
  LDfits$AMM=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataG,Iterations=niters,Algorithm="AMM",
    Thinning=1,Specs=list(Adaptive=niters/10,B=NULL,n=0,Periodicity=2*length(init), w=0.05))
  t3 = proc.time()[['elapsed']]
  
  # Haven't tried yet
  #LDfits$HARM=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataG,Iterations=niters,Algorithm="HARM",
  #  Thinning=1,Specs=list(alpha.star=0.44))
  t4 = proc.time()[['elapsed']]
  
  # Still best
  LDfits$CHARM=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataG,Iterations=niters,Algorithm="CHARM",
    Thinning=1,Specs=list(alpha.star=0.44))
  t5 = proc.time()[['elapsed']]
  
  dosummary = FALSE
  if(dosummary)
  {
    #If it has converged well you will have a Summary2 structure using the ESS:
    LDfit$Summary2[,1]
    
    #If not you can still check Summary1:
    LDfit$Summary1[,1]
    
    #The global fit is very close to the initial LA fit on this occassion.
    magtri(LDfit$Posterior2)
    
    #Otherwise try:
    magtri(LDfit$Posterior1)
    
    #Similarly we can try:
    profitLikeModel(LDfit$Summary2[,1],DataG,makeplots=T)
    
    #Or if that doesn't work:
    profitLikeModel(LDfit$Summary1[,1],DataG,makeplots=T)
  }
}

DataG$algo.func=""
rv = profitLikeModel(bestf,DataG,makeplots=T,cmap = cmap, errcmap=errcmap)

# Now fit the *pre-convolution* best-fit model: hopefully you get the best-fit back!
nopsfmodel = model
nopsfmodel$psf = NULL
nopsftofit = tofit
nopsftofit$psf = NULL
nopsftolog = tolog
nopsftolog$psf = NULL
nopsfpriors = priors
nopsfpriors$psf = NULL
nopsfintervals = intervals
nopsfintervals$psf = NULL

nx = dim(rv$model$z)[1]
ny = dim(rv$model$z)[2]
padx = floor(dim(psf)[1]/2)
pady = floor(dim(psf)[2]/2)
nxpad = nx + 2*padx
nypad = ny + 2*pady
cropx = (1+padx):(nx+padx)
cropy = (1+pady):(ny+pady)

DataM=profitSetupData(image=rv$model$z, mask=mask,sigma=sigma,segim = segim,psf = NULL,
  model = nopsfmodel, tofit = nopsftofit, tolog=nopsftolog, priors = nopsfpriors, 
  intervals=nopsfintervals,algo.func = "", finesample = finesample, verbose=TRUE)
DataM$gain = gain_eff

initp = init
initp['sersic.xcen1'] = init['sersic.xcen1'] + padx
initp['sersic.ycen1'] = init['sersic.ycen1'] + pady

segimpad = matrix(0,nxpad,nypad)
sigmapad = segimpad
maskpad = segimpad
segimpad[(1:nx)+padx,(1:ny)+pady] = segim
sigmapad[(1:nx)+padx,(1:ny)+pady] = sigma
maskpad[(1:nx)+padx,(1:ny)+pady] = mask
DataP=profitSetupData(image=matrix(0,nxpad,nypad), mask=maskpad,sigma=sigmapad,segim = segimpad,psf = NULL,
  model = nopsfmodel, tofit = nopsftofit, tolog=nopsftolog, priors = nopsfpriors, 
  intervals=nopsfintervals, algo.func = "", finesample = finesample, verbose=TRUE)

# Gain to convert from photoelectrons to source photons
gain_atm = 1.0/(gain_inv_e*throughput_atm)
gain_sys = 1.0/(gain_inv_e*throughput_sys)

skylevel = 10^initp['sky.bg']

DataP$usecalcregion=FALSE
bestmodelcounts = (profitLikeModel(initp,DataP)$model$z-skylevel)*gain_eff
bestsky = matrix(skylevel*gain_eff, nrow=nx, ncol=ny)
randmodelcounts = profitPoissonMC(bestmodelcounts,1,throughput_atm,1)
# I don't know why this is faster than just calling rpois but it is
randmodelsky = profitPoissonMC(bestsky,3,throughput_sys,1)
randmodel = (round(randmodelcounts[cropx,cropy]*gain_inv_e)*gain_atm +
  round(randmodelsky*gain_inv_e)*gain_sys)/gain_eff

sigma = sqrt(((randmodel-skylevel)/throughput_atm + skylevel/throughput_sys)/gain_eff)

DataM=profitSetupData(image=randmodel,mask=mask,sigma=sigma,segim = segim,psf = NULL, 
  model = nopsfmodel, tofit = nopsftofit, tolog=nopsftolog, priors = nopsfpriors, 
  intervals=nopsfintervals,algo.func = "", verbose=TRUE)
DataM$gain = gain_eff
DataM$psfarea = NULL

rv = profitLikeModel(init,DataM,makeplots=T,cmap = cmap, errcmap=errcmap)
#LAfitm=LaplaceApproximation(profitLikeModel,parm=init,Data=DataM,Iterations= 1e4,Method='BFGS',CovEst='Identity',sir=FALSE)
#LDfitm=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataM,Iterations=1e4,Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44))

randmodelc = (round(profitBruteConvMC(randmodelcounts,psf,2)[cropx,cropy]*gain_inv_e)*gain_atm + 
  round(randmodelsky*gain_inv_e)*gain_sys)/gain_eff
sigma = sqrt(((randmodelc-skylevel)/throughput_atm + skylevel/throughput_sys)/gain_eff)
DataMC=profitSetupData(image=randmodelc, mask=mask, sigma=sigma, segim = segim, psf = psff,
  model = model, tofit = tofit, tolog=tolog, priors = priors, intervals=intervals, 
  algo.func = "", finesample = finesample, verbose=TRUE)
DataMC$gain = gain_eff
rv = profitLikeModel(init,DataMC,makeplots=T,cmap = cmap, errcmap=errcmap)

makerandc = FALSE
if(makerandc)
{
  nmodels = 1e3
  nxc = nx/4
  nxcpad = nxc + 2*padx
  nyc = ny/4
  nycpad = nyc + 2*pady
  cropmcs = array(dim=c(nxc,nyc,nmodels))
  
  padxci = (floor(nx/2)-floor(nxc/2)-padx+1):(floor(nx/2)+floor(nxc/2)+padx)
  padyci = (floor(ny/2)-floor(nyc/2)-pady+1):(floor(ny/2)+floor(nyc/2)+pady)
  
  cropxc = (1+padx):(nxc+padx)
  cropyc = (1+pady):(nyc+pady)
  
  bmcs = bestmodelcounts[padxci,padyci]
  bsc = bestsky[(floor(nx/2)-floor(nxc/2)+1):(floor(nx/2)+floor(nxc/2)),
    (floor(ny/2)-floor(nyc/2)+1):(floor(ny/2)+floor(nyc/2))]
  for(i in 1:nmodels)
  {
    cropmc[,,i] = (profitBruteConvMC(profitPoissonMC(bmcs,3*i,throughput_atm,1),psf,3*i+1,1,gain_inv_e)[cropxc,cropyc]*gain_atm +
      profitPoissonMC(bsc,3*i+2,throughput_sys,gain_inv_e)*gain_sys)/gain_eff
    if(i %% 50) print(i)
  }
  randvars = matrix(0,nxc,nyc)
  for(j in 1:nyc)
  {
    for(i in 1:nxc)
    {
      randvars[i,j] = var(cropmcs[i,j,])
    }
  }
}

domultimc = TRUE
LDfitms = list()
LDfitmcs = list()
if(domultimc)
{ 
  DataM$algo.func="LD"
  DataMC$algo.func="LD"
  nchains = 20
  niter = 5e3
  tdeconv = 0
  tconv = 0
  for(i in 1:nchains)
  {
    modelcounts = profitPoissonMC(bestmodelcounts,3*i,throughput_atm,1)
    modelcountsc = round(profitBruteConvMC(modelcounts,psf,3*i+1,gain_inv_e))[cropx,cropy]*gain_atm
    modelcounts = round(modelcounts[cropx,cropy]*gain_inv_e)*gain_atm
    skycounts = profitPoissonMC(bestsky,3*i+2,throughput_sys,gain_inv_e)*gain_sys
    DataM$image = (modelcounts + skycounts)/gain_eff
    DataM$sigma = sqrt(((modelcounts+skycounts-skylevel)/throughput_atm + skylevel/throughput_sys)/gain_eff^2)
    t1mc = proc.time()[['elapsed']]
    LDfitms[[i]]=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataM,Iterations=niter,
      Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44))
    t2mc = proc.time()[['elapsed']]
    tdeconv = tdeconv + t2mc - t1mc
    DataMC$image = (modelcounts + skycounts)/gain_eff
    DataMC$sigma = sqrt(((modelcountsc++skycounts-skylevel)/throughput_atm + skylevel/throughput_sys)/gain_eff^2)
    t1mc = proc.time()[['elapsed']]
    LDfitmcs[[i]]=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataMC,Iterations=niter,
      Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44))
    t2mc = proc.time()[['elapsed']]
    tconv = tdeconv + t2mc - t1mc
  }
  
  postms = LDfitms[[1]]$Posterior1
  postmcs = LDfitmcs[[1]]$Posterior1
  for(i in 2:nchains)
  {
    postms = rbind(postms,LDfitms[[1]]$Posterior1)
    postmcs = rbind(postmcs,LDfitms[[1]]$Posterior1)
  }
  
  parnames = c("sersic.mag1","sersic.mag2","sersic.re1","sersic.re2","sersic.nser1","sersic.axrat2","sky.bg")
  profitMagtri(postms[,parnames],inputpar=best[parnames])
  profitMagtri(postmcs[,parnames],inputpar=best[parnames])
}

DataMC2 = DataMC

plotrandc = FALSE
if(plotrandc)
{
  #bestmodelc = profitBruteConv(bestmodel$model$z,psf,matrix(1))
  bmcmean = bestmodelc[100,100]*gain_eff
  tmp = values[100,100,]*gain_eff
  bmcerr = 5*sqrt(bmcmean)
  minx = floor(bmcmean - bmcerr)
  maxx = ceiling(bmcmean + bmcerr)
  dx = 10
  breaks = seq(minx-dx/2,maxx-dx/2,dx)
  y = hist(tmp,breaks=breaks,plot=FALSE)$count
  x = breaks[1:(length(breaks)-1)]
  magplot(x,y/sum(y)/dx,type="s",log="y")
  lines(x+0.5,dchisq(x+0.5,bmcmean),type="s")
  
  DataMC2$sigma = sqrt(randvars)
  rv = profitLikeModel(best,DataMC2,makeplots=T,cmap = cmap, errcmap=errcmap)
}

DataMC2$usecovar = TRUE
DataMC2$algo.func = "estcovar"
rvc = profitLikeModel(best,DataMC2,makeplots=T,cmap = cmap, errcmap=errcmap)
DataMC2$covarinv = rvc$covarinv
DataMC2$algo.func = "LD"

# Fit the single, 
#LAfitmc2=LaplaceApproximation(profitLikeModel,parm=init,Data=DataM,Iterations= 1e4,Method='BFGS',CovEst='Identity',sir=FALSE)
LDfitmc2=LaplacesDemon(profitLikeModel,Initial.Values=init,Data=DataMC2,Iterations=2e3,Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44))

# Make a test grid of values to see how well-behaved LP is with fine changes to input params:
testparams = LAfit_best
gridparams = list()
paramstogrid = c("sersic.mag1","sersic.mag2")
for(param in paramstogrid)
{
  gridparams[[param]] = seq(testparams[[param]]-0.5,testparams[[param]]+0.5,0.025)
}
modelgrid = profitMakeModelGrid(DataG, gridparams, LAfit_best)
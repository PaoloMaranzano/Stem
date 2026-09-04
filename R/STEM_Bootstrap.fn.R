#' @keywords internal
#' @noRd

`STEM_Bootstrap.fn` <-
function(x, StemModel, seed.list.out, output.kriging=NULL,distance='euclidean',precision=0.01,regularization=regularization, verbose = FALSE){

###########
#1) fix a seed x
#2) simulate a new data set z
#3) estimate the parameter vector phi
#4) predict in new spatial locations
##########
	##STEP 1: FIX THE SEED
	#x = seed
	set.seed(x)
	position   = which(unlist(seed.list.out)== x)
	if (isTRUE(verbose)) message("=== Bootstrap iteration n. ", position, " ===")
	#write(position , file = "position.out", append=T)


        ##STEP 2: SIMULATIOM
        StemModel$skeleton$phi = StemModel$estimates$phi.hat   #it follows that the initial values are given by the ML Estimates
	na.pattern	= is.na(StemModel$data$z)
	simulated.z	= STEM_Simulation(StemModel = StemModel,distance=distance)
	### The replicate reproduces the observed design, missing values included:
	### a bootstrap sample with a complete response would understate the
	### uncertainty of a fit obtained from an incomplete one.
	if (any(na.pattern)) simulated.z[na.pattern] = NA_real_
	StemModel$data$z = simulated.z

	###STEP 3: PARAMETER ESTIMATION
	MLE 	= STEM_Estimation(StemModel = StemModel,distance=distance,precision=precision,regularization = regularization, verbose = verbose)

    ##STEP4: SPATIAL PREDICTION
	#data.newlocations	= output.kriging$data.newlocations
	#time.point 			= output.kriging$time.point
	#~ spat.predictions  =  STEM_Kriging(StemModel 		= StemModel,
	#~ output.estimation		= MLE,
	#~ coord.newlocations 	= data.newlocations$coord.newlocations,
	#~ covariates.newlocations= data.newlocations$covariates.newlocations,
	#~ K.newlocations 		= data.newlocations$K.newlocations,
	#~ time.point 			= time.point,
	#~ regular.grid 		= output.kriging$regular.grid)

	return(MLE = MLE) #, spat.predictions=spat.predictions))
}


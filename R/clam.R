
# do for future versions: add greyscale proxy graphs (add dat$model, dat$smooth & calculate required depths...), model acc.rates. Add depths as alternative to depthsfile (but keep depthsfile). 

#' clam
#' 
#' clam is R code for classical (non-Bayesian) age-depth modelling.
#' 
#' @docType package
#' @author Maarten Blaauw <maarten.blaauw@qub.ac.uk>
#' @importFrom grDevices dev.off grey rgb pdf png
#' @importFrom graphics abline image layout legend lines par plot points polygon rect
#' @importFrom stats approx density dnorm lm loess pnorm predict qnorm quantile rnorm runif smooth.spline spline weighted.mean
#' @importFrom utils read.csv read.table write.table packageName
#' @name clam
NULL  


#' @name clam
#' @title The main age-depth modelling function
#' @description Produce age-depth models for cores with dated depths.
#' @details 
#' Cores containing several 14C and/or other dates can be processed semi-automatically in order to obtain age-depth models. 
#' In the process, any 14C dates are calibrated, and age-depth curves are repeatedly drawn through point estimates sampled from the dates.
#'  Age-depth models can be based on linear interpolation, linear/polynomial regression, or cubic, smooth or locally weighted splines. 
#'  For each date, the probability of a calendar year being sampled is proportionate to its calibrated probability (see Blaauw, 2010). 
#'  Uncertainty ranges as well as a 'best' age-model are calculated.
#'
#' Additional cores should be put in a comma-separated file in a sub-folder of the directory where the cores are stored. 
#' By default this parent folder is called \code{coredir="clam_runs"} (if no folder called \code{"Cores"} already exists). If your core is called MyCore1, save MyCore1.csv as \code{clam_runs/MyCore1/MyCore1.csv}.
#' Ensure that the names of the core's folder and filename's root (the part before .csv) match, e.g., using exactly similar upper- and lower case letters.
#' 
#' Avoid the use of spaces or non-standard (non-ASCII) characters within the file or in folder or file names. 
#' The plain text file should consist of 6 or 7 columns (also called fields), containing in the following exact order (see the example below):
#'  \enumerate{
#'   \item Identification labels (e.g. 14C lab codes)
#'   \item 14C ages for 14C-dated depths; leave empty for non-14C dated depths
#'   \item cal BP ages (for any non-14C dates such as the core surface; leave empty for levels with 14C dates)
#'   \item errors (reported 1 standard deviation errors. This column should never be left empty. Errors should always be larger than 0)
#'   \item age offsets if known (otherwise leave empty)
#'   \item depths (depths in the sequence were the dated samples were taken, default unit depth="cm"; this column should never be left empty)
#'   \item thicknesses of the sampled slices (optional column; leave empty for default of 1)
#'  }
#' Add a final empty line to your core's .csv file by pressing 'Enter' after the file's last value. 
#' 
#' These files can be made in spreadsheet software such as MS-Excel, but it is always a good idea to check the file's formatting in a plain-text editor such as WordPad. Remove any lines which contain only commas, and it is also recommended to remove quotes ()\code{\" or \'}) in the headers or elsewhere. 
#'
#' Age-models for the core can then be produced by typing, e.g., \code{clam("MyCore1")}.
#'
#' By default the northern hemisphere terrestrial calibration curve is used (\code{cc=1}, \code{cc1="IntCal13.14C"}). 
#' To use alternative curves, change \code{cc} to \code{cc=2 (cc2="Marine13.14C")}, \code{cc=3 (cc3="SHCal13.14C")}, \code{cc=4 (cc4="mixed.14C")}. 
#' You can also provide custom-built calibration curves, indicating its location using \code{ccdir}.
#'
#' The provided example (default \code{core="Example"}) is core Quilichao-1 which was sampled from a Colombian lake (Berrio et al., 2002). 
#' This core was chosen because it was dated at a rather high resolution, and appears to contain a hiatus (e.g., try \code{hiatus=450}
#' for a hiatus at 450 cm depth).
#'
#' Each clam run will produce a range of files within the core's folder. One, ending with \code{"_calibrated.txt"} contains the calibrated 
#' age ranges of the 14C and other dates. The others will be named according to the core's name followed by the model type, 
#' and contain the age estimates for all depths (files ending with \code{"_ages.txt"}), settings (files ending with \code{"_settings.txt"}) 
#' and graphs (files ending with \code{".pdf"} and \code{".png"}). 
#' The file containing the age estimates has 5 columns; first the depths, then the minima and maxima of the confidence intervals, 
#' then a "best" estimate, and finally the reconstructed accumulation rates. The reported values are rounded to 0 decimals by default
#' (\code{decimals=0}). Accumulation rates are in yr/cm ("deposition time") by default (\code{cmyr=FALSE}), but can be reported in cm/yr (\code{cmyr=TRUE}).
#'
#' see accompanying webpage \url{http://www.chrono.qub.ac.uk/blaauw/clam.html} and Blaauw 2010 (Quaternary Geochronology 5: 512-518).
#'
#' @param core Name of the core, given using quotes. Defaults to the core provided with clam, \code{core="Example"}.
#' @param type The type of age-depth model. Five different types are provided:
#'  \enumerate{
#'   \item linear interpolation between neighbouring levels (1, "int", "inter" or "interp") 
#'   \item linear or higher polynomial regression (2, "reg", "regr", "poly" or "polyn", default linear) 
#'   \item cubic spline (3, "spl" or "spline") 
#'   \item smooth spline (4, "sm" or "smooth", default smoothing 0.3) 
#'   \item locally weighted spline (5, "loess" or "lowess", default smoothing 0.75, cannot extrapolate)
#'  }
#' @param smooth Degree of smoothing. Gives polynomial degree for model type 2. Not relevant for \code{type=1 or type=3}.
#'  \itemize{
#'   \item for type=2: \code{smooth=1} (linear), \code{smooth=2} second-order polynomial, \code{smooth=3} for third-order polynomial, etc. 
#'   \item for type=4: \code{smooth=0.3} 
#'   \item for type=5: \code{smooth=0.75} 
#'  }
#' @param prob Confidence intervals (between 0 and 1), default \code{prob=0.95} or 95\%.
#' @param its Amount of age-model iterations; defaults to \code{its=1000}.
#' @param coredir The directory where core runs are stored (each core in its own directory named after the core's name). 
#' Defaults to \code{coredir="clam_runs"}, or to \code{coredir="Cores"} if this folder exists where R is working. 
#' @param wghts Weights can be applied to dated depths as follows:
#'  \itemize{
#'   \item 0 no weighting 
#'   \item 1 weighted to calibrated probabilities of sampled calendar years (default, \code{wghts=1}). 
#'   \item 2 weighted to (inverse squared) errors of the dates.
#'  }
#' @param cc calibration curve for C14 dates (1, 2 or 3).
#' @param cc1 For terrestrial, northern hemisphere C14 dates.
#' @param cc2 For marine C14 dates.
#' @param cc3 For southern hemisphere C14 dates.
#' @param cc4 For mixed terrestrial/marine C14 dates.
#' @param postbomb Use a postbomb curve for negative (i.e. postbomb) 14C ages. \code{0 = none, 1 = NH1, 2 = NH2, 3 = NH3, 4 = SH1-2, 5 = SH3}. See \url{http://calib.org/CALIBomb/}.
#' @param pb1 For Northern hemisphere region 1 postbomb C14 dates.
#' @param pb2 For Northern hemisphere region 2 postbomb C14 dates.
#' @param pb3 For Northern hemisphere region 3 postbomb C14 dates.
#' @param pb4 For Southern hemisphere regions 1-2 postbomb C14 dates.
#' @param pb5 For Southern hemisphere region 3 postbomb C14 dates.
#' @param ccdir Directory where the calibration curves for C14 dates \code{cc} are located. By default \code{ccdir=""}. 
#' For example, use \code{ccdir="."} to choose current working directory, or \code{ccdir="Curves/"} to choose sub-folder \code{Curves/}.
#' @param outliers The number of any dates to be considered outlying, e.g. \code{c(5,6)} for the fifth and sixth dated depth counting from the top of a core.
#' @param ignore The number of any dates that should be ignored, e.g., \code{c(5,6)} for the fifth and sixth date counting from the top of a core.
#' @param youngest The age beyond which dates should be truncated (e.g., \code{youngest=-60} if the core was sampled in -60 cal BP or AD 2010).
#' @param extradates Depths of any additional dates with their files of ages and probabilities.
#' @param slump Upper and lower depths of sections of abrupt accumulation that should be excised, e.g., \code{c(600, 550, 120, 100)} for two sections of 600-550 and 120-100 cm depth.
#' @param est Which point estimate to use as 'best' age. It is highly recommended to not only use these 'best' point estimates, as chronological uncertainties are often considerable and should not be ignored.  
#'  \enumerate{
#'   \item averages of age-depth model derived ages (default, \code{est=1})
#'   \item midpoints of age-depth model derived age estimates
#'   \item midpoints of calibrated ranges
#'   \item weighted means of calibrated ranges 
#'   \item medians of calibrated distributions 
#'   \item maximum densities of calibrated distributions 
#'   \item midpoints of entire calibrated distributions (including years outside the calibrated ranges)
#'  }
#' @param calibt Off by default; provide two parameters such as \code{c(3,4)}.
#' @param mixed.effect Set to \code{TRUE} to activate mixed-effect modelling.
#' @param dmin Minimum depth of age-depth model (e.g., extrapolate).
#' @param dmax Maximum depth of age-depth model (e.g., extrapolate).
#' @param every Resolution at which (ages for) depths are calculated.
#' @param yrmin Minimum of calendar axis of age-depth plot (calculate automatically by default).
#' @param yrmax Maximum of calendar axis of age-depth plot (calculated automatically by default).
#' @param yrsteps Temporal resolution at which calibrated ages are calculated (in calendar years).
#' @param pbsteps Temporal resolution at which postbomb C14 ages are calibrated (in calendar years).
#' @param hpdsteps Temporal resolution at which highest posterior density ranges are calibrated (in calendar years).
#' @param BCAD Use BC/AD or cal BP scale.
#' @param decimals Amount of decimals for rounding.
#' @param cmyr Accumulation rates can be provided as yr/cm (default, \code{cmyr=TRUE}, more accurately named deposition times) or cm/yr (\code{cmyr=FALSE}).
#' @param ageofdepth Calculate age estimates of a specific depth.
#' @param depth Depth units.
#' @param depthseq Sequence of depths for which age estimates are to be calculated (default: from \code{dmin} to \code{dmax} with steps of size every)
#' @param depths.file Use a file with depths for depthseq.
#' @param thickness Thickness of the dated samples.
#' @param hiatus Depths of any hiatuses, e.g., \code{c(500, 300)}. Each sub-section must have at least 2 dates (\code{4} for smoothing spline; does not work with loess as it cannot extrapolate).
#' @param remove.reverse Proportion of age-models with reversals that can be removed before prompting a warning.
#' Set at \code{FALSE} to avoid removing models with reversals.
#' @param times Half-range of calibration curve used to calibrate dates (multiplication factor for the dates' errors).
#' @param sep Separator between the fields of the plain text file containing the dating information.
#' @param ext Extension of the file containing the dating information.
#' @param runname Text to add to the core name for specific runs, e.g., "MyCore_Test1"
#' @param storedat Store the dates and age-model within R after a \code{clam} run.
#' @param threshold Below which value should probabilities be excluded from calculations.
#' @param proxies Set to \code{TRUE} to plot proxies against age after the run.
#' @param revaxes Set to \code{TRUE} to plot ages on the vertical axis and depth on the horizontal axis.
#' @param revd Plot depth axis in reverse.
#' @param revyr Plot age axis in reverse.
#' @param calhght Heights of the calibrated distributions in the age-depth plot.
#' @param maxhght Maximum height of age probability distributions.
#' @param mirror Plot the age distributions in "mirror" style (above and below depth).
#' @param plotrange Plot the confidence ranges of the age-model.
#' @param mar Plot margins (amount of white space along edges of axes 1-4).
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted).
#' @param bty Type of box to be drawn around plots. Draw a box around the graph (\code{"n"} for none, 
#' and \code{"l"}, \code{"7"}, \code{"c"}, \code{"u"}, "]" or \code{"o"} for correspondingly shaped boxes).
#' @param plotpdf Produce a pdf file of the age-depth plot.
#' @param plotpng Produce a png file of the age-depth plot.
#' @param greyscale Produce a grey-scale representation of all age-models (number gives resolution, e.g., 500 bins; will cancel plotting of the confidence intervals).
#' @param yrlab Alternative names can be provided.
#' @param dlab Alternative names can be provided.
#' @param calcol Colour of the calibrated distributions in the age-depth plot.
#' @param C14col Colour of the calibrated ranges of the dates.
#' @param outcol Colour of outlying dates.
#' @param outlsize Size of symbols outlying dates.
#' @param bestcol Colour of the "best" age-depth model (based on chosen value for est).
#' @param rangecol Colour of plotted confidence ranges.
#' @param slumpcol Colour of slump.
#' @param plotname Print the core name on the graph.
#' @param ash Plot all distributions at the same height.
#' 
#' @author Maarten Blaauw
#' @return Age model construction together with a text output and files saved to a folder in the \code{coredir/core} directory.
#' @examples 
#'  clam(, coredir=tempdir()) # Create the example in Cores/Example folder
#'  clam(, coredir=tempdir(), extradates=470) 
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/clam.html}
#' \link{calibrate}
#' \link{mix.calibrationcurves}
#' \link{pMC.age}
#' \link{age.pMC}
#' \link{student.t}
#' \link{deptime.depth}
#' \link{deptime.age}
#' \link{plot_proxies}
#' 
#' @references
#' Berrio, J.C., Hooghiemstra, H., Marchant, R., Rangel, O., 2002. Late-glacial and Holocene history of the dry forest area 
#' in the south Colombian Cauca Valley. _Journal of Quaternary Science_ 17, 667-682
#'
#' Blaauw, M., 2010. Methods and code for 'classical' age-modelling of radiocarbon sequences. Quaternary Geochronology 5, 512-518
#' \url{http://dx.doi.org/10.1016/j.quageo.2010.01.002}
#' @export
clam <- function(core="Example", type=1, smooth=c(), prob=0.95, its=1000, coredir=c(), wghts=1, cc=1, cc1="IntCal13.14C", cc2="Marine13.14C", cc3="SHCal13.14C", cc4="mixed.14C",  postbomb=FALSE, pb1="postbomb_NH1.14C", pb2="postbomb_NH2.14C", pb3="postbomb_NH3.14C", pb4="postbomb_SH1-2.14C",pb5="postbomb_SH3.14C", ccdir="", outliers=c(), ignore=c(), youngest=c(), extradates=c(), slump=c(), est=1, calibt=FALSE, mixed.effect=FALSE, dmin=c(), dmax=c(), every=1, yrmin=c(), yrmax=c(), yrsteps=1, pbsteps=0.01, hpdsteps=1, BCAD=FALSE, decimals=0, cmyr=FALSE, ageofdepth=c(), depth="cm", depthseq=c(), depths.file=FALSE, thickness=1, hiatus=c(), remove.reverse=0.5, times=5, sep=",", ext=".csv", runname=c(), storedat=TRUE, threshold=1e-6, proxies=FALSE, revaxes=FALSE, revd=TRUE, revyr=TRUE, calhght=0.3, maxhght=0.01, mirror=TRUE, plotrange=TRUE, bty="l", mar=c(3.5,3,2,1), mgp=c(2,1,0), plotpdf=TRUE, plotpng=TRUE, greyscale=c(), yrlab=c(), dlab=c(), calcol=rgb(0,0.5,0.5,0.5), C14col=rgb(0,0,1,0.5), outcol="red", outlsize=1, bestcol="black", rangecol=rgb(0,0,0,0.3), slumpcol=grey(0.75), plotname=TRUE, ash=FALSE) {
  # If coredir is left empty, check for a folder named Cores in the current working directory, and if this doesn't exist, for a folder called Bacon_runs (make this folder if it doesn't exist yet).
  # Check if we have write access. If not, tell the user to provide a different, writeable location for coredir. 
  if(length(coredir) == 0) {
    if(dir.exists("Cores")) 
      coredir <- "Cores" else
        if(dir.exists("clam_runs"))
          coredir <- "clam_runs" else {
            coredir <- "clam_runs"			  
            wdir <- dir.create(coredir, FALSE)
            if(!wdir)
              stop("Cannot write into the current directory.\nPlease set coredir to somewhere where you have writing access, e.g. Desktop or ~.")
		    }    
  } else {
    if(!dir.exists(coredir))
      wdir <- dir.create(coredir, FALSE)
	  if(!dir.exists(coredir)) # if it still doesn't exist, we probably don't have enough permissions
        stop("Cannot write into the current directory.\nPlease set coredir to somewhere where you have writing access, e.g. Desktop or ~.")
	}
  coredir <- .validateDirectoryName(coredir)
  cat("The run's files will be put in this folder: ", coredir, core, "\n", sep="")
	  
  # Copy Example's files into core's directory
  if(core == "Example") {
	if(!dir.exists(paste(coredir, "Example/", sep="")))  
      dir.create(paste(coredir, "Example/", sep=""), showWarnings = FALSE, recursive = TRUE)
    fileCopy <- system.file("extdata/Example/", package="clam")
    file.copy(fileCopy, coredir, recursive = TRUE, overwrite=FALSE)
  } 
    
  # warn and stop if unexpected settings are provided
  if(type > 5 || type < 1 || prob < 0 || prob > 1 || its < 100 || wghts < 0 || wghts > 1 || est < 1 || est > 7 || yrsteps <= 0 || hpdsteps <= 0 || every <= 0 || decimals < 0 || cmyr < 0 || cmyr > 1 || thickness < 0 || times < 1 || calhght < 0 || (type==5 && length(hiatus)>0))
    stop("\n Warning, clam cannot run with these settings! Please check the manual.\n\n", call.=FALSE)

  dets <- suppressWarnings(read.csv(paste(coredir, core, "/", core, ext, sep=""), sep=sep))
  d <- dets[,6]
  if(min(diff(d)) < 0)
    cat("\n Warning, depths not in ascending order (top ones should come first).\n\n")

  # avoid confusing warning when using sample for first time in session. Probably not necessary any more
  # tmp <- suppressWarnings(sample(1:1e3, 1, prob=rep(.001,1e3), replace=TRUE))

  # avoid Windows/Mac habit of silently adding .txt extension to plain text files
  Gates <- list.files(paste(coredir, core, sep=""), pattern=".csv.txt")
  if(length(Gates) > 0) {
    cat("\nRemoving unnecessary .txt extension from .csv file", Gates[1], "\n")
    file.rename(paste(coredir, core, "/", core, ".csv.txt", sep=""),
      paste(coredir, core, "/", core, ".csv", sep=""))
  }

  # set the calibration curve
  ccdir <- .validateDirectoryName(ccdir)
  if(ccdir == "") # so, if no alternative folder provided, use clam's calibration curves
    ccdir = paste(system.file("extdata", package=packageName()), "/", sep="")
    
  if(cc==1) calcurve <- read.table(paste(ccdir, cc1,  sep="")) else
    if(cc==2) calcurve <- read.table(paste(ccdir, cc2,  sep="")) else
      if(cc==3) calcurve <- read.table(paste(ccdir, cc3,  sep="")) else
        if(cc==4) calcurve <- read.table(paste(ccdir, cc4,  sep="")) else
          stop("I do not understand which calibration curve you mean, check the manual", call.=FALSE)
  if(cc==1) ccname <- cc1 else
    if(cc==2) ccname <- cc2 else
      if(cc==3) ccname <- cc3 else
        if(cc==4) ccname <- cc4 

  # negative C14 ages need a postbomb curve
  pbnames <- c(pb1, pb2, pb3, pb4, pb5)
  cdat <- dets[,2]
  if(length(cdat[!is.na(cdat)]) > 0)
    if(min(cdat[!is.na(cdat)]) < 0)
      if(postbomb==FALSE)
        cat("Warning, negative 14C ages, should I use a postbomb curve?") else {
          if(postbomb>5)
            stop("I do not understand which postbomb curve you mean, check the manual", call.=FALSE)
          yrsteps <- min(pbsteps, yrsteps)
          pb <- read.table(system.file("extdata", pbnames[postbomb], package=packageName()))
          pb.x <- seq(min(pb[,1]), max(pb[,1]), by=yrsteps)
          pb.y <- approx(pb[,1], pb[,2], pb.x)$y
          pb.sd <- approx(pb[,1], pb[,3], pb.x)$y
          calcurve <- cbind(c(pb.x, calcurve[,1]), c(pb.y, calcurve[,2]), c(pb.sd, calcurve[,3]))
        }

  # work in BC/AD if needed, and prepare for calculations in f14C
  if(BCAD) {
    theta <- 1950-calcurve[,1]
    border <- max(which(theta > 0))
    theta <- c(theta[1:(border-1)], theta[border]:theta[border+2], theta[(border+3):length(theta)])
    mu <- approx(1950-calcurve[,1], calcurve[,2], theta)$y
    sigma <- approx(1950-calcurve[,1], calcurve[,3], theta)$y
    theta[theta <=0] <- theta[theta <= 0]-1
    calcurve <- cbind(theta, mu, sigma)
  } else theta <- calcurve[,1]
  if(length(yrlab) == 0)
    yrlab <- ifelse(BCAD, "cal BC/AD", "cal BP")
  f.mu <- exp(-calcurve[,2]/8033)
  f.sigma <- exp(-(calcurve[,2]-calcurve[,3])/8033) - f.mu

  # prepare for slumps and hiatuses
  #cat(paste("Core name Test:", core))
  if(length(greyscale) > 0) storedat <- TRUE
  if(length(slump) > 0) {
    if(length(slump) %% 2 == 1)
      stop("\n Warning, slumps need both upper and lower depths. Please check the manual", call.=FALSE)
    slump <- matrix(sort(slump), ncol=2, byrow=TRUE)
    if(length(dmax)==0)
      dmax <- max(dets[,6])
    if(length(extradates) > 0)
      dmax <- max(dmax, extradates)
    for(i in 1:nrow(slump)) {
      d[d > min(slump[i,])] <- d[d > min(slump[i,])] - (max(slump[i,]) - min(slump[i,]))
      dmax <- dmax - (max(slump[i,])-min(slump[i,]))
    }
  if(length(hiatus) > 0)
    for(i in 1:nrow(slump)) {
       below.slump <- which(hiatus > max(slump[i,]))
       above.slump <- which(hiatus < min(slump[i,]))
       hiatus[below.slump] <- hiatus[below.slump] - (max(slump[i,])-min(slump[i,]))
       hiatus <- hiatus[c(above.slump, below.slump)]
    }
  }

  # read in the data
  dat <- .read.clam(core, coredir, ext, hpdsteps, yrsteps, prob, times, sep, BCAD, storedat, ignore, thickness, youngest, slump, threshold, theta, f.mu, f.sigma, calibt, extradates, calcurve, postbomb)
  cat("\n Calibrating dates... ")

  # calculate the depths to be used, based on the ranges and resolution
  if(length(dmin)==0)
    dmin <- floor(min(dat$depth))
  if(length(dmax)==0)
    dmax <- ceiling(max(dat$depth))
  if(depths.file)
    if(file.exists(dd <- paste(coredir, core, "/", core, "_depths.txt", sep=""))) {
      if(length(depthseq) == 0)
        depthseq <- seq(dmin, dmax, by=every)
      depthseq <- sort(unique(c(depthseq, suppressWarnings(read.table(dd))[,1])))
      dmin <- min(depthseq)#, read.table(dd)[,1])
      dmax <- max(depthseq)#, read.table(dd)[,1])
    } else
      stop(paste("\nCannot find file ", dat$core, "_depths.txt!\n", sep=""), call.=FALSE)
  if(length(depthseq) == 0)
    depthseq <- seq(dmin, dmax, by=every)
  if(proxies) {
    storedat <- TRUE
    if(file.exists(dd <- paste(coredir, core, "/", core, "_proxies.csv", sep="")))
      dat$proxies <- suppressWarnings(read.csv(dd, sep=sep)) else
        stop(paste("\nCannot find file ", dat$core, " _proxies.csv!\n", sep=""), call.=FALSE)
    dmin <- min(depthseq, dat$proxies[,1])
    dmax <- max(depthseq, dat$proxies[,1])
    depthseq <- sort(unique(c(depthseq, dat$proxies[,1])))
  }

  if(length(ageofdepth) > 0)
    depthseq <- sort(unique(c(ageofdepth, depthseq)))

  # decide which models and point estimates should be used
  if(any(type==c(1, "int", "inter", "interp"))) type <- 1 else
    if(any(type==c(2, "reg", "regr", "poly", "polyn"))) type <- 2 else
      if(any(type==c(3, "spline", "spl"))) type <- 3 else
        if(any(type==c(4, "smooth", "sm"))) type <- 4 else
          if(any(type==c(5, "loess", "lowess"))) type <- 5
  best <- cbind(dat$mid1, dat$mid1, dat$mid1, dat$wmn, dat$med, dat$mode, dat$mid2)
  Est <- best[,est]  

  # remove outliers from the age-depth modelling
  if(length(outliers) > 0) {
    depths <- dat$depth[-outliers]
    errors <- dat$error[-outliers]
    calibs <- dat$calib[-outliers]
    Est <- Est[-outliers]
  } else {
      depths <- dat$depth
      errors <- dat$error
      calibs <- dat$calib
    }

  # age-depth modelling with curves through sampled age estimates
  # in sections if one or more hiatuses are present
  if(length(hiatus) > 0) {
    allrange <- c(0,0,0,0)
    hiatusseq <- sort(c(range(depthseq), hiatus))
    for(i in 2:length(hiatusseq)) {
      cat(paste("\n section ", i-1, ",", sep=""))
      section <- depthseq[min(which(depthseq >= hiatusseq[i-1])) : max(which(depthseq <= hiatusseq[i]))]
      if(i>2) section <- section[-1]
      sel <- min(which(depths >= min(section))):max(which(depths <= max(section)))
      if(mixed.effect)
        if(length(outliers) > 0)
          smp <- .mixed.effect(its, depths, dat$cal[-outliers], dat$cage[-outliers], errors, calibs, Est, theta, f.mu, f.sigma, yrsteps, calibt) else
            smp <- .mixed.effect(its, depths, dat$cal, dat$cage, errors, calibs, est, theta, f.mu, f.sigma, yrsteps, calibt) else
              smp <- .smpl(its, depths[sel], calibs[sel], Est[sel])
      calrange <- .model.clam(type, smooth, its, wghts, depths[sel], errors[sel], section, prob, est, dat, smp, greyscale, remove.reverse, storedat, ageofdepth, BCAD)
      allrange <- rbind(allrange, calrange)
    }
  calrange <- allrange[2:nrow(allrange),]
  } else {
    if(mixed.effect)
	  if(length(outliers) > 0)
	    smp <- .mixed.effect(its, depths, dat$cal[-outliers], dat$cage[-outliers], errors, calibs, est, theta, f.mu, f.sigma, yrsteps, calibt) else
	      smp <- .mixed.effect(its, depths, dat$cal, dat$cage, errors, calibs, Est, theta, f.mu, f.sigma, yrsteps, calibt) else
	        smp <- .smpl(its, depths, calibs, Est)
     calrange <- .model.clam(type, smooth, its, wghts, depths, errors, depthseq, prob, est, dat, smp, greyscale, remove.reverse, storedat, ageofdepth, BCAD)
    }
  dat$model <- approx(calrange[,1], (calrange[,2]+calrange[,3])/2, dat$depth)$y

  if(est==2) calrange[,4] <- (calrange[,2]+calrange[,3])/2
  if(!BCAD && any(diff(calrange[,4]) < 0) || BCAD && any(diff(calrange[,4]) > 0))
    reversal <- TRUE else reversal <- FALSE
  gfit <- round(.gfit(theta, f.mu, f.sigma, dat, calrange, outliers), 2)

  # re-correct the depths if slumps were applied
  if(length(slump) > 0) {
    dat <- .read.clam(core, coredir, ext, hpdsteps, yrsteps, prob, times, sep, BCAD, storedat, ignore, thickness, youngest, slump=c(), threshold, theta, f.mu, f.sigma, calibt, extradates, calcurve, postbomb) # read in the original dates again
    calrange <- calrange[which(calrange[,1] <= dmax),]
    d <- calrange[,1]
    for(i in 1:nrow(slump)) {
      d[d > min(slump[i,])] <- d[d > min(slump[i,])] + (max(slump[i,]) - min(slump[i,]))
      dmax <- dmax + (max(slump[i,]) - min(slump[i,]))
      calrange[,1] <- d
      hiatus[hiatus > min(slump[i,])] <- hiatus[hiatus > min(slump[i,])] + (max(slump[i,]) - min(slump[i,]))
    }
  }

  # produce the age-depth plot, and a pdf copy if desired
  if(length(yrmin)==0)
    yrmin <- min(dat$mid1, calrange[,2])
  if(length(yrmax)==0)
    yrmax <- max(dat$mid1, calrange[,3])
  if(length(ageofdepth > 0))
    layout(matrix(c(1,2,1,3), nrow=2), heights=c(.7,.3))
  .ageplot(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, yrlab, dlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, if(length(greyscale)>0) get('chron') else c(), C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, core, bty, mar, mgp, ash)

  # write files providing calibrated dates, age-model and settings
  colnames(calrange) <- c("Depth", paste("min.", 100*prob, "%range", sep=""), paste("max.", 100*prob, "%range", sep=""), "point")
  .write.clam(dat, coredir, runname, calrange, core, prob, type, remove.reverse, smooth, wghts, its, outliers, ignore, est, BCAD, yrsteps, every, decimals, cmyr, depth, depthseq, hiatus, gfit, reversal, plotpdf, plotpng, yrmin, yrmax, dmin, dmax, dlab, yrlab, plotrange, greyscale, if(length(greyscale)>0) get('chron') else c(), C14col, outcol, outlsize, bestcol, rangecol, calhght, maxhght, mirror, calcol, slump, slumpcol, revaxes, revyr, revd, calibt, youngest, extradates, plotname, calcurve, ccname, postbomb, pbnames, depths.file, bty, mar, mgp, ash)
  closeAllConnections()

  if(storedat) {
    calrange <<- calrange
    dat <<- dat
    smp <<- smp
  }

  # plot the age distribution of a provided depth
  if(length(ageofdepth) > 0) {
    if(revaxes)
      abline(v=ageofdepth, lty=2) else
        abline(h=ageofdepth, lty=2)
    xlim <- range(.ageofdepth)
    if(!BCAD) xlim <- xlim[2:1]
    .ageofdepth <- get(".ageofdepth")
    hst <- density(.ageofdepth, n=max(1, max(xlim)-min(xlim)))
    yr <- seq(min(xlim), max(xlim), by=yrsteps)
    hst <- cbind(c(yr, max(xlim), min(xlim)), c(approx(hst$x, hst$y, yr)$y, 0, 0))
    plot(hst, type="n", main="", xlim=xlim, xlab=yrlab, ylab="")
    polygon(hst, col="grey")
    legend("topleft", paste(ageofdepth, depth), bty="n")
    layout(matrix(1))
    rng <- round(calrange[max(which(calrange[,1] <= ageofdepth)),])
    cat("\n  Age range of ", ageofdepth, " ", depth, ": ", rng[3], " to ", rng[2], ifelse(BCAD, " cal BC/AD", " cal BP"), " (", rng[3]-rng[2], " yr, ", prob, " % range)",  sep="")
  }

  # report the confidence ranges, the goodness-of-fit, and whether any age-reversals occurred
  rng <- round(calrange[,3]-calrange[,2])
  cat("\n  ", core, "'s ", 100*prob, "% confidence ranges span from ", min(rng), " to ", max(rng), " yr (average ", round(mean(rng)), " yr)", sep="")
  cat("\n  Fit (-log, lower is better):", gfit, "\n")
  if(reversal) 
    cat("  Age reversals occurred. Try other model?\n")
}



#' @name mix.calibrationcurves
#' @title Build a custom-made, mixed calibration curve. 
#' @description If two curves need to be 'mixed' to calibrate, e.g. for dates of mixed terrestrial and marine carbon sources, then this function can be used. 
#' @details The proportional contribution of each of both calibration curves has to be set. 
#'
#' @param proportion Proportion of the first calibration curve required. e.g., change to \code{proportion=0.7} if \code{cc1} should contribute 70\% (and \code{cc2} 30\%) to the mixed curve.
#' @param cc1 The first calibration curve to be mixed. Defaults to the northern hemisphere terrestrial curve IntCal13.
#' @param cc2 The second calibration curve to be mixed. Defaults to the marine curve IntCal13.
#' @param name Name of the new calibration curve.
#' @param dirname Directory where the file will be written. If using the default \code{dirname="."}, 
#' the new curve will be saved in current working directory. 
#' @param offset Any offset and error to be applied to \code{cc2} (default 0 +- 0).
#' @author Maarten Blaauw, J. Andres Christen
#' @return A file containing the custom-made calibration curve, based on calibration curves \code{cc1} and \code{cc2}.
#' @examples
#'   mix.calibrationcurves(dirname=tempdir())
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
mix.calibrationcurves <- function(proportion=.5, cc1="IntCal13.14C", cc2="Marine13.14C", name="mixed.14C", dirname=".", offset=c(0,0)) {
  ccloc <- normalizePath(system.file("extdata/", package='clam')) 
  dirname <- .validateDirectoryName(dirname)
  
  cc1 <- read.table(normalizePath(paste(ccloc, "/", cc1,  sep="")))
  cc2 <- read.table(normalizePath(paste(ccloc, "/", cc2,  sep="")))
  cc2.mu <- approx(cc2[,1], cc2[,2], cc1[,1], rule=2)$y + offset[1] # interpolate cc2 to the calendar years of cc1
  cc2.error <- approx(cc2[,1], cc2[,3], cc1[,1], rule=2)$y
  cc2.error <- sqrt(cc2.error^2 + offset[2]^2)
  mu <- proportion * cc1[,2] + (1-proportion) * cc2.mu
  error <- proportion * cc1[,3] + (1-proportion) * cc2.error
  write.table(cbind(cc1[,1], mu, error), paste(dirname, name,  sep="") , row.names=FALSE, col.names=FALSE, sep="\t")
}



#' @name pMC.age
#' @title Calculate C14 ages from pmC values.
#' @description Calculate C14 ages from pmC values of radiocarbon dates.
#' @details Post-bomb dates are often reported as pMC or percent modern carbon. Since Bacon expects radiocarbon ages,
#'  this function can be used to calculate radiocarbon ages from pMC values. The reverse function is \link{age.pMC}.
#' @param mn Reported mean of the pMC.
#' @param sdev Reported error of the pMC.
#' @param ratio Most modern-date values are reported against \code{100}. If it is against \code{1} instead, use \code{1} here.
#' @param decimals Amount of decimals required for the radiocarbon age.
#' @author Maarten Blaauw, J. Andres Christen
#' @return Radiocarbon ages from pMC values. If pMC values are above 100\%, the resulting radiocarbon ages will be negative.
#' @examples
#'   pMC.age(110, 0.5) # a postbomb date, so with a negative 14C age
#'   pMC.age(80, 0.5) # prebomb dates can also be calculated
#'   pMC.age(.8, 0.005, 1) # pMC expressed against 1 (not against 100\%)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#'  \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
pMC.age <- function(mn, sdev, ratio=100, decimals=0) {
  y <- -8033 * log(mn/ratio)
  sdev <- y - -8033 * log((mn+sdev)/ratio)
  round(c(y, sdev), decimals)
}



#' @name age.pMC
#' @title Calculate pMC values from C14 ages
#' @description Calculate pMC values from radiocarbon ages
#' @details Post-bomb dates are often reported as pMC or percent modern carbon. Since Bacon expects radiocarbon ages, 
#' this function can be used to calculate pMC values from radiocarbon ages. The reverse function of \link{pMC.age}.
#' @param mn Reported mean of the 14C age.
#' @param sdev Reported error of the 14C age.
#' @param ratio Most modern-date values are reported against \code{100}. If it is against \code{1} instead, use \code{1} here.
#' @param decimals Amount of decimals required for the pMC value.
#' @author Maarten Blaauw, J. Andres Christen
#' @return pMC values from C14 ages.
#' @examples
#'   age.pMC(-2000, 20)
#'   age.pMC(-2000, 20, 1)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
age.pMC <- function(mn, sdev, ratio=100, decimals=3) {
  y <- exp(-mn / 8033)
  sdev <- y - exp(-(mn + sdev) / 8033)
  signif(ratio*c(y, sdev), decimals)
}



#' @name add.dates
#' @title Add dates to age-depth plots
#' @description Add dated depths to plots, e.g. to show dates that weren't used in the age-depth model
#' @details Sometimes it is useful to add additional dating information to age-depth plots, e.g., to show outliers or how dates calibrate with different estimated offsets. 
#' @param mn Reported mean of the date. Can be multiple dates. 
#' @param sdev Reported error of the date. Can be multiple dates. 
#' @param depth Depth of the date. 
#' @param cc The calibration curve to use: \code{cc=1} for IntCal13 (northern hemisphere terrestrial), \code{cc=2} for Marine13 (marine), \code{cc=0} for none (dates that are already on the cal BP scale).
#' @param above Threshold for plotting of probability values. Defaults to \code{above=1e-3}.
#' @param exx Exaggeration of probability distribution plots. Defaults to \code{exx=50}.
#' @param normal By default, Bacon uses the student's t-distribution to treat the dates. Use \code{normal=TRUE} to use the normal/Gaussian distribution. This will generally give higher weight to the dates.
#' @param normalise By default, the date is normalised to an area of 1 (\code{normalise=TRUE}). 
#' @param t.a The dates are treated using the student's t distribution by default (\code{normal=FALSE}). 
#' The student's t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010). 
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file). 
#' For symmetry reasons, t.a must always be equal to t.b-1. 
#' @param t.b The dates are treated using the student's t distribution by default (\code{normal=FALSE}). 
#' The student's t-distribution has two parameters, t.a and t.b, set at 3 and 4 by default (see Christen and Perez, 2010). 
#' If you want to assign narrower error distributions (more closely resembling the normal distribution), set t.a and t.b at for example 33 and 34 respectively (e.g., for specific dates in your .csv file). 
#' For symmetry reasons, t.a must always be equal to t.b-1. 
#' @param age.res Resolution of the date's distribution. Defaults to \code{date.res=100}.
#' @param times The extent of the range to be calculated for each date. Defaults to \code{times=20}.
#' @param col The colour of the ranges of the date. Default is semi-transparent red: \code{col=rgb(1,0,0,.5)}.
#' @param border The colours of the borders of the date. Default is semi-transparent red: \code{border=rgb(1,0,0,0.5)}.
#' @param rotate.axes The default of plotting age on the horizontal axis and event probability on the vertical one can be changed with \code{rotate.axes=TRUE}.
#' @param mirror Plot the dates as 'blobs'. Set to \code{mirror=FALSE} to plot simple distributions.
#' @param up Directions of distributions if they are plotted non-mirrored. Default \code{up=TRUE}.
#' @param BCAD The calendar scale of graphs is in \code{cal BP} by default, but can be changed to BC/AD using \code{BCAD=TRUE}. 
#' @author Maarten Blaauw, J. Andres Christen
#' @return A date's distribution, added to an age-depth plot.
#' @examples
#'   clam(coredir=tempfile())
#'   add.dates(5000, 100, 60)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @references
#' Blaauw, M. and Christen, J.A., Flexible paleoclimate age-depth models using an autoregressive 
#' gamma process. Bayesian Anal. 6 (2011), no. 3, 457--474. 
#' \url{https://projecteuclid.org/download/pdf_1/euclid.ba/1339616472}
#' @export
add.dates <- function(mn, sdev, depth, cc=1, above=1e-3, exx=50, normal=TRUE, normalise=TRUE, t.a=3, t.b=4, age.res=100, times=20, col=rgb(1,0,0,.5), border=rgb(1,0,0,.5), rotate.axes=FALSE, mirror=TRUE, up=TRUE, BCAD=FALSE) {
  if(cc > 0)
    cc = copyCalibrationCurve(cc)
  
  for(i in 1:length(mn)) {
	yrs <- seq(mn[i]-times*sdev[i], mn[i]+times*sdev[i], length=age.res)
	if(length(cc) < 2)
	  cc <- cbind(yrs, yrs, rep(0, length(yrs)))	
	ages <- approx(cc[,1], cc[,2], yrs)$y
	errors <- approx(cc[,1], cc[,3], yrs)$y
	
	if(normal)
	  probs <- dnorm(ages, mn[i], sqrt(sdev[i] + errors)^2) else
	    probs <- (t.b + (mn[i]-ages)^2  / (2*(sdev[i]^2 + errors^2))) ^ (-1*(t.a+0.5))
	if(normalise)
	  probs <- probs / sum(probs)	
	these <- which(probs >= above)
	if(length(these) > 0) {
	  yrs <- yrs[these]
	  probs <- probs[these]	
	}
	  	
	if(!up)
	  up <- -1  
	if(BCAD)
      yrs <- 1950 - yrs
	if(mirror)
  	  pol <- cbind(c(yrs, rev(yrs)), depth[i] + exx*c(probs, -rev(probs))) else
        pol <- cbind(c(min(yrs), yrs, max(yrs)), depth[i] - up*exx*c(0, probs,  0))
	if(rotate.axes)
	  pol <- pol[,2:1]		
	polygon(pol, col=col, border=border)
  }
}



#' @name copyCalibrationCurve
#' @title Copy a calibration curve.
#' @description Copy one of the the calibration curves into memory. 
#' @details Copy the radiocarbon calibration curve defined by cc into memory. 
#' @return The calibration curve (invisible).
#' @param cc Calibration curve for 14C dates: \code{cc=1} for IntCal13 (northern hemisphere terrestrial), \code{cc=2} for Marine13 (marine), 
#' \code{cc=3} for SHCal13 (southern hemisphere terrestrial). 
#' @param postbomb Use \code{postbomb=TRUE} to get a postbomb calibration curve (default \code{postbbomb=FALSE}).
#' @author Maarten Blaauw, J. Andres Christen
#' @examples 
#' intcal13 <- copyCalibrationCurve(1)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/manualBacon_2.3.pdf}
#' @export
copyCalibrationCurve <- function(cc=1, postbomb=FALSE) {
  if(postbomb) {
    if(cc==1) fl <- "postbomb_NH1.14C" else
      if(cc==2) fl <- "postbomb_NH2.14C" else
        if(cc==3) fl <- "postbomb_NH3.14C" else
          if(cc==4) fl <- "postbomb_SH1-2.14C" else 
            if(cc==5) fl <- "postbomb_SH3.14C" else  
              stop("Calibration curve doesn't exist\n")		       
  } else
  if(cc==1) fl <- "IntCal13.14C" else
    if(cc==2) fl <- "Marine13.14C" else
      if(cc==3) fl <- "SHCal13.14C" else
        stop("Calibration curve doesn't exist\n") 		 
  cc <- system.file("extdata", fl, package='clam')
  cc <- read.table(cc)
  invisible(cc)
}



# See Christen and Perez 2009, Radiocarbon 51:1047-1059. Instead of assuming the standard Gaussian model (default in clam), a student t distribution can be used with two parameters. Christen and Perez 2009 suggest t.a = 3 and t.b = 4; this can be put as clam( calibt=c(3,4) )
.calibt <- function(t.a, t.b, f.cage, f.error, theta, f.mu, f.sigma)
  (t.b + ((f.cage-f.mu)^2) / (2*(f.sigma^2 + f.error^2))) ^ (-1*(t.a+0.5))


#' @name student.t 
#' @title Comparison dates calibrated using both the student-t distribution and the the normal distribution.
#' @description Visualise how a date calibrates using the student-t distribution and the the normal distribution.
#' @details Radiocarbon and other dates are usually modelled using the normal distribution (red curve). The student-t approach (grey distribution) however allows for wider tails and thus tends to better accommodate outlying dates. This distribution requires two parameters, called 'a' and 'b'.
#' @param y The reported mean of the date.
#' @param error The reported error of the date.
#' @param t.a Value for the student-t parameter \code{a}.
#' @param t.b Value for the student-t parameter \code{b}.
#' @param postbomb Which postbomb curve to use for negative 14C dates
#' @param cc calibration curve for C14 dates (1, 2 or 3).
#' @param cc1 For northern hemisphere terrestrial C14 dates.
#' @param cc2 For marine C14 dates.
#' @param cc3 For southern hemisphere C14 dates.
#' @param cc4 A custom calibration curve
#' @param ccdir Directory where the calibration curves for C14 dates \code{cc} are allocated. By default \code{ccdir=""}. 
#' Use \code{ccdir="."} to choose current working directory. Use \code{ccdir="Curves/"} to choose sub-folder \code{Curves/}.
#' @param Cutoff Threshold above which calibrated probabilities are plotted
#' @param times 8 by default.
#' @author Maarten Blaauw
#' @examples 
#' student.t() 
#' 
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/clam.html}
#' @references
#' Christen, J.A., P{\'e}rez, S. E., 2009. A new robust statistical model for radiocarbon data. Radiocarbon 51 (3), 1047-1059
#' \url{https://journals.uair.arizona.edu/index.php/radiocarbon/article/view/3562}
#' @export
student.t <- function(y=2450, error=50, t.a=3, t.b=4, cc=1,  postbomb=c(), cc1="IntCal13", cc2="Marine13", cc3="SHCal13", cc4="mixed", ccdir="",Cutoff=1e-5, times=8)
{
  ccdir <-.validateDirectoryName(ccdir)
  # set the calibration curve
  if(ccdir=="")
    ccdir = paste(system.file("extdata", package=packageName()), "/", sep="")

  
  if(cc == 0)
  {
    x <- seq(y-(times*error), y+(times*error), length=500)
    norm.cal <- dnorm(x, y, error)
    norm.cal <- norm.cal/sum(norm.cal)
    t.cal <- (t.b + ((y-x)^2) / (2*(error^2))) ^ (-1*(t.a+0.5))
    t.cal <- t.cal/sum(t.cal)
    t.cal <- cbind(c(min(x), x, max(x)), c(0, t.cal, 0))
    plot(x, norm.cal, type="l", xlab="cal BP", xlim=range(x)[2:1], ylab="", ylim=c(0, max(t.cal[,2], norm.cal)), col=2, lwd=1.5)
    polygon(t.cal, col=rgb(0,0,0,.25), border=rgb(0,0,0,.5))
    legend("topleft", "Gaussian", text.col=2, bty="n")
    legend("topright", paste("student-t (a=", t.a, ", b=", t.b, ")", sep=""), bty="n", text.col=grey(.4))
  } else
  {
    
    if(cc1=="IntCal13") cc1 <- read.table(paste(ccdir, "IntCal13.14C",  sep="")) else
      cc1 <- read.csv(paste(ccdir, cc1,  sep=""))[,1:3]
    if(cc2=="Marine13") cc2 <- read.table(paste(ccdir, "Marine13.14C",  sep="")) else
      cc2 <- read.csv(paste(ccdir, cc2,  sep=""))[,1:3]
    if(cc3=="SHCal13") cc3 <- read.table(paste(ccdir, "SHCal13.14C",  sep="")) else
      cc3 <- read.table(paste(ccdir, cc3,  sep=""))[,1:3]
    if(cc4=="mixed") cc4 <- read.table(paste(ccdir, "mixed.14C",  sep="")) else
      cc4 <- read.table(paste(ccdir, cc4,  sep=""))[,1:3]
    if(cc==1) cc <- cc1 else if(cc==2) cc <- cc2 else if(cc==3) cc <- cc3 else cc <- cc4

    if (y < 0)
      if (length(postbomb) == 0) 
        stop("Warning, negative ages require a postbomb curve. Provide value for postbomb")
    else 
    {
      
      if(postbomb==1) bomb <- read.table(system.file("extdata","postbomb_NH1.14C", package=packageName()))[,1:3] else
        if(postbomb==2) bomb <- read.table(system.file("extdata","postbomb_NH2.14C", package=packageName()))[,1:3] else
          if(postbomb==3) bomb <- read.table(system.file("extdata","postbomb_NH3.14C", package=packageName()))[,1:3] else
            if(postbomb==4) bomb <- read.table(system.file("extdata","postbomb_SH1-2.14C", package=packageName()))[,1:3] else
              if(postbomb==5) bomb <- read.table(system.file("extdata","postbomb_SH3.14C", package=packageName()))[,1:3] else
                stop("Warning, cannot find postbomb curve #", postbomb, " (use values of 1 to 5 only)")
              
      bomb.x <- seq(max(bomb[, 1]), min(bomb[, 1]), length = 500)
      bomb.y <- approx(bomb[, 1], bomb[, 2], bomb.x)$y
      bomb.z <- approx(bomb[, 1], bomb[, 3], bomb.x)$y
      bomb <- cbind(bomb.x, bomb.y, bomb.z, deparse.level = 0)
      #if (info$postbomb < 4) #JEV warning
      if (postbomb < 4) 
        cc <- rbind(bomb, cc1, deparse.level = 0)
      else cc <- rbind(bomb, cc3, deparse.level = 0)
    }
    
    norm.cal <- dnorm(cc[, 2], y, sqrt(cc[, 3]^2 + error^2))
    norm.cal <- cbind(cc[, 1], norm.cal/sum(norm.cal))
    acc <- which(norm.cal[, 2] >= Cutoff)
    if(y < 0) acc <- 1:200 # feo pero funciona
    norm.cal <- norm.cal[acc, ]
    
    t.cal <- (t.b + ((y - cc[, 2])^2)/(2 * (cc[, 3]^2 + error^2)))^
      (-1 * (t.a + 0.5))
    t.cal <- cbind(cc[, 1], t.cal/sum(t.cal))
    acc <- which(t.cal[, 2] >= Cutoff)
    if(y < 0) acc <- 1:200 # feo pero funciona
    t.cal <- t.cal[acc, ]
    
    plot(norm.cal, type = "l", xlab = "cal BP", xlim = range(c(t.cal[,1], norm.cal[, 1]))[2:1], ylab = "",
	  ylim = c(0, max(t.cal[,2], norm.cal[, 2])), col = 2, lwd = 1.5)
    polygon(t.cal, col = rgb(0, 0, 0, 0.25), border = rgb(0, 0, 0, 0.5))
    legend("topleft", "Gaussian", text.col = 2, bty = "n")
    legend("topright", paste("student-t (a=", t.a, ", b=", t.b, ")", sep = ""), bty = "n", text.col = grey(0.4))
  }     

}
# find the calibrated distributions of 14C dates
.caldist <- function(f.cage, f.error, theta, f.mu, f.sigma, yrsteps, threshold, calibt, BCAD, normalise=FALSE)
  {
    if(f.cage > 1)
      {
        if(f.cage > 1) yrsteps <- min(yrsteps, .1)
        pb <- theta[which(f.mu > 1)]
        if(length(pb)==0)
          stop("help, something exploded with a postbomb date")
        x <- approx(theta, f.mu, seq(min(pb), max(pb), by=yrsteps))
        xsd <- approx(theta, f.sigma, x$x)$y
        theta <- c(x$x, theta[which(f.mu <= 0)])
        f.mu <- c(x$y, f.mu[which(f.mu <= 0)])
        f.sigma <- c(xsd, f.sigma[which(f.mu <= 0)])
        threshold <- 0
      }

    # calibrate; find how far f.cage (measurement) is from f.mu (calibration curve)
    if(length(calibt) < 2)
      cal <- cbind(theta, dnorm(f.mu, f.cage, sqrt(f.error^2+f.sigma^2))) else
        cal <- cbind(theta, .calibt(calibt[1], calibt[2], f.cage, f.error, theta, f.mu, f.sigma))

    # interpolate and normalise calibrated distribution to 1
    cal <- cal[min(which(cal[,2] > 0)):max(which(cal[,2] > 0)),] # remove unnecessary data
    cal <- approx(cal[,1], cal[,2], seq(min(cal[,1]), max(cal[,1]), by=yrsteps))
    cal <- cbind(cal$x, cal$y/sum(cal$y))
    if(BCAD && (0 %in% cal[,1]))
			   cal <- cal[-which(cal[,1]==0),] # 0 BC/AD does not exist
    # only report those normalised calibrated probabilities beyond a threshold
    cal[cal[,2] > threshold,]
  }


# find the highest posterior density (hpd) of the calibrated distribution
.hpd <- function(dat, prob, hpdsteps, yrsteps)
  {
    # interpolate and rank the ages according to their calibrated distribution probabilities
    dat <- approx(dat[,1], dat[,2], seq(min(dat[,1]), max(dat[,1]), by=yrsteps))
    o <- order(dat$y, decreasing=TRUE)
    dat <- cbind(dat$x[o], dat$y[o]/sum(dat$y))

    # only retain those ages with cumulative normalised probabilities within required percentage
    dat <- dat[which(cumsum(dat[,2]) <= prob),]
    dat <- dat[order(dat[,1]),]

    # identify any individual ranges within the hpd range and calculate their probability
    dif <- which(diff(dat[,1]) > hpdsteps)
    if(length(dif)==0)
      hpds <- cbind(min(dat[,1]), max(dat[,1]), 100*prob) else
        {
          dif <- c(dat[1,1], sort(c(dat[dif,1], dat[dif+1,1])), dat[nrow(dat),1])
          dif <- matrix(dif, ncol=2, byrow=TRUE)
          probs <- c()
          for(i in 1:nrow(dif))
            probs[i] <- round(100*sum(dat[which(dat[,1]==dif[i,1]):which(dat[,1]==dif[i,2]),2]), 1)
          hpds <- cbind(dif, probs)
        }
    hpds
  }


# calculate the age-depth model and its uncertainty
.model.clam <- function(type, smooth, its, wghts, depths, errors, depthseq, prob, est, dat, smp, greyscale, remove.reverse, storedat, ageofdepth, BCAD)
  {
    # warn for extrapolation, refuse to do so for loess
    if(min(depthseq) < min(dat$depth) || max(depthseq) > max(dat$depth))
      if(type==5)
        stop(" cannot extrapolate using loess! Change settings.\n ", call.=FALSE) else
          cat(" extrapolating beyond dated levels, dangerous!\n ")

    # choose model: interpolation, (polynomial) regression, spline, smooth spline or loess
    chron <- array(0, dim=c(length(depthseq), its))
    if(type==1) chron <- .interp(depthseq, depths, its, chron, smp) else
      if(type==2) chron <- .poly(depthseq, smooth, wghts, errors, depths, its, chron, smp) else
        if(type==3) chron <- .spline(depthseq, smooth, depths, its, chron, smp) else
          if(type==4) chron <- .smooth(depthseq, smooth, wghts, errors, depths, its, chron, smp) else
            if(type==5) chron <- .loess(depthseq, smooth, wghts, errors, depths, its, chron, smp)

    # test against age reversals
    warp <- c()
    if(remove.reverse!=FALSE)
      for(i in 1:ncol(chron))
        if(!BCAD && min(diff(chron[,i])) <= 0 || BCAD && max(diff(chron[,i])) >= 0)
          warp <- c(warp, i)
    if(length(warp) > 0)
      if(length(warp) > remove.reverse*its)
        cat("\n\n !!! Too many models with age reversals!!!\n") else
          {
            cat("\n Removing", length(warp), "models with age reversals,", its-length(warp), "models left...")
            chron <- chron[,-warp]
            smp <- smp[,-warp,]
          }

    if(length(ageofdepth) > 0)
      if(ageofdepth %in% depthseq)
        .assign_to_global(".ageofdepth", chron[which(depthseq==ageofdepth),]) 

    if(storedat)
      {
        chron <<- chron
        smp <<- smp
      }

    # find uncertainty ranges of calendar age for each depth of the core
    calrange <- array(0, dim=c(nrow(chron), 2))
    mn <- c()
    for(i in 1:nrow(chron))
      {
        x <- chron[i,2:ncol(chron)]
        qp <- (1-prob)/2
        calrange[i,] <- quantile(x, c(qp, 1-qp))
        if(est==1) mn[i] <- mean(x)
      }
    if(est==1)
      cbind(depthseq, cbind(calrange, mn)) else
      cbind(depthseq, cbind(calrange, chron[,1]))
  }


# sample point age estimates from the calibrated distributions ('its' times)
# the probability of a year being sampled is proportional to its calibrated probability
.smpl <- function(its, depths, calibs, Est)
  {
    smp <- array(1, dim=c(length(depths), 1+its, 2))
    smp[,1,1] <- Est
    for(i in 1:length(calibs))
      smp[i,(1:its)+1,] <-
        calibs[[i]][sample(1:length(calibs[[i]][,1]), its, prob=calibs[[i]][,2], TRUE),]
    smp
  }

# akin to Heegaard et al.'s mixed effect modelling, but using calibrated dates
.mixed.effect <- function(its, depths, cals, cages, errors, calibs, Est, theta, f.mu, f.sigma, yrsteps, calibt)
  {
    cat("\n Mixed effect modelling, this will take some time")
    smp <- array(1, dim=c(length(depths), 1+its, 2))
    smp[,1,1] <- Est
    for(i in 1:length(cals))
      if(!is.na(cals[i]))
        if(length(calibt)==0)
          {
            x <- rnorm(its, cals[i], errors[i])
            smp[i,(1:its)+1,] <- c(x, dnorm(x, cals[i], errors[i]))
          } else
            {
              x <- (cals[i]-10*errors[i]) : (cals[i]+10*errors[i])
              x <- cbind(x, .calibt(calibt[1], calibt[2], cals[i], errors[i], x, x, 0))
              o <- order(x[,2], decreasing=TRUE)
              x <- cbind(x[o,1], cumsum(x[o,2])/sum(x[,2]))
              sampled.x <- max(which(x[,2] <= runif(1, 0, max(x[,2]))))
              smp[i,(1:its)+1,] <- x[sampled.x,]
            } else
             for(j in 1:its)
               {
                 if(j/(its/3) == round(j/(its/3))) cat(".")
                 yr <- rnorm(1, cages[i], errors[i])
                 f.yr <- exp(-yr/8033)
                 f.error <- f.yr - exp(-(yr+errors[i])/8033)
                 yr <- cbind(theta, dnorm(f.mu, f.yr, sqrt(f.error^2+f.sigma^2)))
                 yr <- yr[yr[,2]>0,]
                 yr <- approx(yr[,1], yr[,2], seq(min(yr[,1]), max(yr[,1]), by=yrsteps))
                 smp.yr <- sample(length(yr$x), 1, prob=yr$y)
                 smp[i,j+1,] <- c(yr$x[smp.yr], yr$y[smp.yr])
               }
    smp
  }


# interpolate linearly between the data (default)
.interp <- function(depthseq, depths, its, chron, smp)
  {
    cat(" Interpolating, sampling")
    for(i in 1:its)
      {
        temp <- approx(depths, smp[,i,1], depthseq, ties=mean)$y

        # allow for extrapolation... dangerous!
        if(min(depthseq) < min(depths))
          {
            minus <- which(depthseq < min(depths))
            slope <- diff(temp)[max(minus)+1]/diff(depthseq)[max(minus)+1]
            temp[minus] <- temp[max(minus)+1] + slope * (depthseq[minus] - min(depths))
          }
        if(max(depthseq) > max(depths))
          {
            maxim <- which(depthseq > max(depths))
            slope <- diff(temp)[min(maxim)-2]/diff(depthseq)[min(maxim)-2]
            temp[maxim] <- temp[min(maxim)-1] + slope * (depthseq[maxim] - max(depths))
          }
        chron[,i] <- temp
        if(i/(its/5) == round(i/(its/5))) cat(".")
      }
    chron
  }


# polynomial regressions of certain order through the data (default linear, y=ax+b)
.poly <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp)
  {
    if(length(smooth)==0)
      cat(" Using linear regression, sampling") else
      cat(paste(" Using polynomial regression (degree ", smooth, "), sampling", sep=""))
    if(wghts==0) w <- c() else w <- 1/errors^2
    for(i in 1:its)
      {
        if(wghts==1) w <- smp[,i,2]
        chron[,i] <- predict(lm(smp[,i,1] ~ poly(depths, max(1, smooth)), weights=w), data.frame(depths=depthseq))
        if(i/(its/5) == round(i/(its/5))) cat(".")
      }
    chron
  }


# fit cubic spline interpolations through the data
.spline <- function(depthseq, smooth, depths, its, chron, smp)
  {
    if(length(smooth) < 1) smooth <- .3
    cat(paste(" Using cubic spline sampling", sep=""))
    for(i in 1:its)
      {
        chron[,i] <- spline(depths, smp[,i,1], xout=depthseq)$y
        if(i/(its/5) == round(i/(its/5))) cat(".")
      }
    chron
  }


# fit cubic smoothed splines through the data, with smoothing factor
.smooth <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp)
  {
    if(length(smooth) < 1) smooth <- .3
    cat(paste(" Using smoothing spline (smoothing ", smooth, "), sampling", sep=""))
    if(wghts==0) w <- c() else w <- 1/errors^2
    for(i in 1:its)
      {
        if(wghts==1) w <- smp[,i,2]
        chron[,i] <- predict(smooth.spline(depths, smp[,i,1], w=w, spar=smooth), depthseq)$y
        if(i/(its/5) == round(i/(its/5))) cat(".")
      }
    chron
  }


# fit locally weighted (1/errors^2) splines through the data, with smoothing factor
.loess <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp)
  {
    if(length(smooth) < 1) smooth <- .75
    cat(paste(" Using loess (smoothing ", smooth, "), sampling", sep=""))
    if(wghts==0) w <- c() else w <- 1/errors^2
    for(i in 1:its)
      {
        if(wghts==1) w <- smp[,i,2]
        chron[,i] <- predict(loess(smp[,i,1] ~ depths, weights=w, span=smooth), depthseq)
        if(i/(its/5) == round(i/(its/5))) cat(".")
      }
    chron
  }


# read the data and perform first calculations incl. calibrations
.read.clam <- function(name, namedir,ext, hpdsteps, yrsteps, prob, times, sep, BCAD, storedat, ignore, thickness, youngest, slump, threshold, theta, f.mu, f.sigma, calibt, extradates, calcurve, postbomb)
  {
    # read the file with the dating information
    dat <- list(coredir=paste(namedir, name, "/", sep=""), name=name)
    if(!file.exists(paste(namedir, name, sep="")))
      stop(paste("\n\n Warning, cannot find a folder within",namedir," named ", name, ". Have you saved it in the right place and with the right name? Please check the manual\n\n", sep=""), call.=FALSE)
    if(!file.exists(paste(dat$coredir, name, ext, sep="")))
      stop(paste(" \n\n Warning, cannot find file ", name, ".csv in folder",namedir, name, ". Have you saved it in the right place and named it correctly? Please check the manual\n\n", sep=""), call.=FALSE)
    dets <- suppressWarnings(read.table(paste(dat$coredir, name, ext, sep=""), comment.char="", header=TRUE, sep=sep, na.strings = c("#N/A!", "NA", "@NA")))

    # ignore dates if required, add thickness column if it was left out
    if(length(ignore) > 0)
      {
        dat$ignore <- as.character(dets[ignore,1])
        dets <- dets[-ignore,]
      }
    if(ncol(dets) < 7)
      dets <- cbind(dets, thickness) else
      dets[is.na(dets[,7]),7] <- thickness

    # should slumps be taken into account?
    if(length(slump) > 0)
      {
        d.adapt <- dets[,6]
        d.lost <- c()
        for(i in 1:nrow(slump))
          {
            below.slump <- which(dets[,6] > max(slump[i,]))
            above.slump <- which(dets[,6] < min(slump[i,]))
            d.lost <- c(d.lost, which(!(1:nrow(dets) %in% c(above.slump, below.slump))))
            d.adapt[below.slump] <- d.adapt[below.slump] - (max(slump[i,])-min(slump[i,]))
          }
        dets[,6] <- d.adapt
        if(length(d.lost) > 0)
          dets <- dets[-d.lost,]
      }
    # check for common errors
    dets <- dets[,1:7]
    x <- 0
    for(i in 2:7) if(is.factor(dets[,i])) x <- 1
    if(x == 1)
		stop(paste("\n Some value fields in ", name, ".csv contain letters, please adapt", sep=""), call.=FALSE)
    if(length(dets[is.na(dets[,2]),2])+length(dets[is.na(dets[,3]),3]) != nrow(dets))
      stop(paste("\n Remove duplicate entries within the C14 and calendar fields in ", name, ".csv", sep=""), call.=FALSE)
    if(min(dets[,4]) <= 0)
      stop(paste("\n Errors of dates should be larger than zero. Please adapt ", name, ".csv", sep=""), call.=FALSE)
    dat$ID <- as.character(dets[,1])

    # correct for any reservoir effect
    dets[is.na(dets[,5]),5] <- 0
    dat$cage <- dets[,2] - dets[,5]
    dat$error <- dets[,4]

    # work in F14C for calibration
    dat$f.cage <- exp(-dat$cage/8033)
    dat$f.error <- dat$f.cage - exp(-(dat$cage+dat$error)/8033)

    # check if any 14C dates are (entirely or partly) beyond the calibration curve
    outside <- which(!is.na(dat$cage))
    rangecc <- c(min(calcurve[,2]-calcurve[,3]),max(calcurve[,2]+calcurve[,3]))
    outside <- outside[c(which(dat$cage[outside]-times*dat$error[outside] < rangecc[1]), which(dat$cage[outside]+times*dat$error[outside] > rangecc[2]))]
    if(length(outside) > 0)
      {
        truncate <- 0
        for(i in 1:length(outside)) # check if date lies only partly beyond the curve limits
          if((dat$cage[outside[i]]-times*dat$error[outside[i]] < rangecc[1] &&
            dat$cage[outside[i]]+times*dat$error[outside[i]] > rangecc[1]) ||
            (dat$cage[outside[i]]-times*dat$error[outside[i]] < rangecc[2] &&
            dat$cage[outside[i]]+times*dat$error[outside[i]] > rangecc[2]))
              truncate <- truncate + 1
        if(truncate > 0)
          cat("\n Warning, dates spanning beyond the calibration curve will be truncated! ")

        # remove dates which lie entirely outside the limits of the calibration curve
        outside <- outside[c(which(dat$cage[outside]+qnorm(1-(1-prob)/2)*dat$error[outside] < rangecc[1]), which(dat$cage[outside]-qnorm(1-(1-prob)/2)*dat$error[outside] > rangecc[2]))]
        if(length(outside) > 0)
          {
            cat("\n Warning, dates older than the calibration curve will be ignored! ")
            dets <- dets[-outside,]
            dat$cage <- dat$cage[-outside]
            dat$error <- dat$error[-outside]
            dat$f.cage <- dat$f.cage[-outside]
            dat$f.error <- dat$f.error[-outside]
            dat$outside <- dat$ID[outside]
            dat$ID <- dat$ID[-outside]
          }
      }

    # fill the 'dat' list with additional information
    dat$cal <- c(dets[,3], extradates)
    dat$res <- c(dets[,5], extradates)
    dat$depth <- c(dets[,6], extradates)
    dat$thick <- c(dets[,7], rep(thickness, length(extradates)))
    dat$BCAD <- BCAD

    # find distribution (calibrated if 14C) and point estimates for each date
    for(i in 1:length(dat$depth))
      {
        if(length(extradates) > 0 && i > nrow(dets))
          {
            tmp <- read.table(paste(dat$coredir, name, "_", extradates[i-nrow(dets)], ".txt", sep=""))
            calib <- cbind(tmp[,1], tmp[,2]/sum(tmp[,2]))
          } else
            if(is.na(dat$cage[[i]]))
              {
                age <- dat$cal[[i]]
                error <- dat$error[[i]]
                ageseq <- seq(age-(times*error), age+(times*error), by=yrsteps)
                calib <- cbind(ageseq, dnorm(ageseq, age, error))
              } else
                calib <- .caldist(dat$f.cage[[i]], dat$f.error[[i]], theta, f.mu, f.sigma, yrsteps, threshold, calibt, BCAD)
        if(length(youngest) > 0) # truncate ages younger than a limit
          {
            if(BCAD) calib <- calib[which(calib[,1] <= youngest),] else
              calib <- calib[which(calib[,1] >= youngest),]
            if(length(calib) == 0)
              if(BCAD)
                calib <- cbind(seq(youngest-(3*yrsteps), youngest+yrsteps, length=5), c(0:3,0)/3) else
                calib <- cbind(seq(youngest-yrsteps, youngest+(3*yrsteps), length=5), c(0,3:0)/3)
          }
        dat$calib[[i]] <- calib
        dat$hpd[[i]] <- .hpd(calib, prob=prob, hpdsteps=hpdsteps, yrsteps=yrsteps)
        dat$mid1[[i]] <- (dat$hpd[[i]][1] + dat$hpd[[i]][2*nrow(dat$hpd[[i]])])/2
        yrs <- calib[,1]
        dat$mid2[[i]] <- mean(c(max(yrs), min(yrs)))
        dat$wmn[[i]] <- weighted.mean(calib[,1], 1/calib[,2])
        dat$med[[i]] <- calib[max(which(cumsum(calib[,2]) <= .5)),1]
        dat$mode[[i]] <- calib[which(calib[,2] == max(calib[,2])),1][1]
      }

    if(storedat)
	  dets <<- dets
    dat
  }


# calculate goodness-of-fit (small number, so calculate its -log)
.gfit <- function(theta, f.mu, f.sigma, dat, calrange, outliers)
  {
    gfit <- c()
    if(length(outliers) > 0)
      {
        dat$cage <- dat$cage[-outliers]
        dat$error <- dat$error[-outliers]
        dat$cal <- dat$cal[-outliers]
        dat$model <- dat$model[-outliers]
      }
    gfit <- pnorm(dat$cal, dat$model, dat$error^2)
    if(length(c14 <- which(!is.na(dat$cage))) > 0) # if there are radiocarbon dates
      {
        gfit.c <- approx(theta, f.mu, dat$model[c14])$y # C14 age at cc of modelled cal date
        f.cage <- exp(-dat$cage[c14]/8033)
        f.error <- exp(-(dat$cage[c14]-dat$error[c14])/8033) - f.cage
        gfit.var <- f.error^2 + approx(theta, f.sigma, dat$model[c14])$y^2
        gfit[c14] <- pnorm(f.cage, gfit.c, sqrt(gfit.var)) # deviation between measured and cc ages
      }
    dat$gfit <- -sum(log(gfit[!is.na(gfit)]))
  }


# write files of the age-depth model, calibrated ranges, and settings
.write.clam <- function(dat, namedir,runname, calrange, name, prob, type, remove.reverse, smooth, wghts, its, outliers, ignore, est, BCAD, yrsteps, every, decimals, cmyr, depth, depthseq, hiatus, gfit, reversal, plotpdf, plotpng, yrmin, yrmax, dmin, dmax, yrlab, dlab, plotrange, greyscale, chron, C14col, outcol, outlsize, bestcol, rangecol, calhght, maxhght, mirror, calcol, slump, slumpcol, revaxes, revyr, revd, calibt, youngest, extradates, plotname, calcurve, ccname, postbomb, pbnames, depths.file, bty, mar, mgp, ash)
  {
    # age-depth model; age estimates, accumulation rates and ranges for every analysed depth
    runnames <- c("_interpolated", "_polyn_regr", "_cubic_spline", "_smooth_spline", "_loess")
    calrange <- cbind(calrange, round(c(diff(calrange[,4])/diff(calrange[,1]), NA), decimals+2))
    if(cmyr)
      calrange[,5] <- 1/calrange[,5]
    calrange[,2:4] <- round(calrange[,2:4], decimals)
    ifelse(length(runname)==0, runname <- runnames[type], runname)
    if(depths.file && file.exists(dd <- paste(namedir, name, "/", name, "_depths.txt", sep="")))
      {
        dd <- read.table(dd)[,1]
        this <- c()
        for(i in 1:length(dd))
          this[i] <- which(calrange[,1]==dd[i])[1] # find where the relevant ages are
	    write.table(calrange[this,], paste(dat$coredir, name, runname, "_ages.txt", sep=""), row.names=FALSE, col.names=c("depth", paste("min", 100*prob, "%", sep=""), paste("max", 100*prob, "%", sep=""), "best", "acc.rate"), quote=FALSE, sep="\t")
      } else
       write.table(calrange, paste(dat$coredir, name, runname, "_ages.txt", sep=""), row.names=FALSE, col.names=c("depth", paste("min", 100*prob, "%", sep=""), paste("max", 100*prob, "%", sep=""), "best", "accrate"), quote=FALSE, sep="\t")

    # calibrated ranges of all dates
    hpd.file <- file(paste(dat$coredir, name, "_calibrated.txt", sep=""), "w")
    cat(paste("Calibrated age ranges at ", 100*prob, "% confidence intervals\n", sep=""), file=hpd.file)
    for(i in 1:length(dat$depth))
      {
        cat(paste("\n\nDepth: ", dat$depth[[i]], "\nyrmin\tyrmax\tprobability\n"), file=hpd.file)
        hpds <- dat$hpd[[i]]
        for(j in 1:nrow(hpds))
          {
            for(k in 1:3) cat(hpds[j,k], "\t", file=hpd.file)
            cat("\n", file=hpd.file)
          }
      }
    close(hpd.file)

    # relevant settings and results
    set.file <- file(paste(dat$coredir, name, runnames[type], "_settings.txt", sep=""), "w")
    cat(paste("Settings (square brackets give names of the constants)\n\n",
      "Calibration curve: ", ccname,
      if(postbomb!=FALSE)
		paste(",", pbnames[postbomb], "for postbomb dates"),
      "\nAge-depth model: ",
      if(type==1) "linear interpolation between dated levels [type=1]" else
      if(type==2) ifelse(length(smooth)==0, "linear regression [type=2, smooth=c()]",
      paste("polynomial regression [type=2] of order", smooth, "[smooth]")) else
      if(type==3) "cubic spline [type=3]" else
      if(type==4) paste("smooth spline [type=4] with spar =", ifelse(length(smooth)<1, 0.3, smooth), "[smooth]") else
      if(type==5) paste("locally weighted spline [type=5] with span =", ifelse(length(smooth)<1, 0.75, smooth), "[smooth]"),
      if(wghts==1) "\nWeighted by the calibrated probabilities [wghts=1]",
      if(wghts==2) "\nWeighted by the errors (1/sdev^2) [wghts=2]",
      "\nCalculations at ", 100*prob, "% confidence ranges [prob=", prob, "]",
      "\nAmount of iterations: ", its, " [its]",
      "\nCalendar age point estimates for depths based on ",
      if(est==1) "weighted average of all age-depth curves [est=1]" else
      if(est==2) "midpoints of the hpd ranges of the age-depth curves [est=2]" else
      if(est==3) "midpoints of the hpd ranges of the dated levels [est=3]" else
      if(est==4) "weighted means of the dated levels [est=4]" else
      if(est==5) "medians of the dated levels [est=5]" else
      if(est==6) "modes/maxima/intercepts of the dated levels [est=6]",
      "\nCalendar scale used: ", if(BCAD) "cal BC/AD" else "cal BP",
      " [BCAD=", BCAD, "] at a resolution of ", yrsteps, " yr [yrsteps]",
      "\nAges were calculated every ", every, " [every] ", depth,
      " [depth], from ", min(depthseq), " [dmin] to ", max(depthseq), " [dmax] ", depth, sep=""), file=set.file)
      if(length(youngest) > 0) cat("\n\nDates with ages younger than", youngest, ifelse(BCAD, "BC/AD", "cal BP"), "were truncated", file=set.file)
      if(length(calibt)> 1) cat("\n\nInstead of assuming the standard Gaussian model, a student t distribution was used with t.a =", calibt[1], "and t.b =", calibt[2], "(see Christen and Perez 2009, Radiocarbon 51:1047-1059)", file=set.file)
    if(length(slump) == 2) cat("\n\nA slump was excised between", max(slump), "and", min(slump), depth, file=set.file)
    if(length(slump) > 2)
      {
        cat("\n\nSlumps were excised from ", file=set.file)
        sl <- array(sort(slump), dim=c(2, length(slump)/2))
        for(i in 1:ncol(sl))
          cat(sl[1,i], "to", sl[2,i], depth, if(i<ncol(sl)) "and ", file=set.file)
      }
    if(length(outliers) > 0)
      {
        cat("\n\nDates assumed outlying [outliers]: ", file=set.file)
        for(i in outliers) cat(i, " (", dat$ID[i], ") ", sep="", file=set.file)
      }
    if(length(ignore) > 0)
      {
        cat("\n\nDates ignored [ignore]: ", file=set.file)
        for(i in 1:length(ignore)) cat(ignore[i], " (", dat$ignore[i], ") ", sep="", file=set.file)
      }
    if(length(dat$outside) > 0)
      {
        cat("\n\nDates outside calibration curve and ignored: ", file=set.file)
        for(i in 1:length(dat$outside)) cat(dat$outside[i], " ", sep="", file=set.file)
      }
    cat(paste(
      if(length(hiatus) > 0)
        paste("\nA hiatus was inferred at", hiatus, depth, "[hiatus]"),
        "\n\nGoodness-of-fit (-log, lower is better): ", gfit,
      if(reversal) "\nSome age-depth reversals occurred"),
      if(remove.reverse) "\nAny models with age-depth reversals were removed",
      "\n\nProduced ", date(), sep="", file=set.file)
    close(set.file)

    if(plotpdf)
      {
        pdf(file=paste(dat$coredir, name, runname, ".pdf", sep=""))
        .ageplot(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, dlab, yrlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, if(length(greyscale)>0) chron else c(), C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, name, bty, mar, mgp, ash)
        dev.off()
      }
    if(plotpng)
      {
        png(filename = paste(dat$coredir, name, runname, ".png", sep=""))
        .ageplot(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, dlab, yrlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, if(length(greyscale)>0) chron else c(), C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, name, bty, mar, mgp, ash)
        dev.off()
      }
  }


.ageplot <- function(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, yrlab, dlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, chron, C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, name, bty="l", mar, mgp, ash=FALSE)
  {
    # set up initial parameters
    if(length(dlab)==0) dlab <- paste("Depth (", depth, ")", sep="")
    ifelse(BCAD || !revyr, yr.lim <- c(yrmin, yrmax), yr.lim <- c(yrmax, yrmin))
    if(revd) d.lim <- c(dmax, dmin) else d.lim <- c(dmin, dmax)

    par(xaxt="s", xaxs="r", yaxt="s", yaxs="r", bty=bty, mar=mar, mgp=mgp, font=2)
    if(revaxes) plot(0, type="n", ylim=yr.lim, xlim=d.lim, xlab=dlab, ylab=yrlab) else
      plot(0, type="n", xlim=yr.lim, ylim=d.lim, xlab=yrlab, ylab=dlab)
    if(plotname) legend("topleft", name, bty="n")

    # draw histograms of all age-depth models. Off by default, time-consuming!
    if(length(greyscale)==1)
      {
       	plotrange <- FALSE
        depgr=seq(dmin, dmax, length=greyscale)
        for(i in 2:greyscale)
          {
            temp <- density(chron[max(which(calrange[,1]<=depgr[i])),2:ncol(chron)], n=greyscale)
            if(revaxes)
              image(c(depgr[i-1], depgr[i]), temp$x, matrix(temp$y), col=grey(1-(0:100)/100), add=TRUE) else
                image(temp$x, c(depgr[i-1], depgr[i]), matrix(temp$y), col=grey(1-(0:100)/100), add=TRUE)
          }
      }

    # draw the age-depth models, per section if hiatuses were inferred
    if(length(hiatus) > 0)
      {
        if(length(slump) == 0)
          hiatusseq <- sort(c(range(depthseq), hiatus)) else
            hiatusseq <- sort(c(range(depthseq, depthseq+sum(slump[,2]-slump[,1])), hiatus))
        for(i in 2:length(hiatusseq))
          {
            sec <- calrange[min(which(calrange[,1] > hiatusseq[i-1])):max(which(calrange[,1] < hiatusseq[i])),]
            pol <- cbind(c(sec[,2], rev(sec[,3])), c(sec[,1], rev(sec[,1])))
            if(plotrange)
              if(revaxes)
                polygon(pol[,2], pol[,1], col=rangecol, border=rangecol) else
                  polygon(pol, col=rangecol, border=rangecol)
            if(revaxes)
              lines(sec[,1], sec[,4], lwd=2, col=bestcol) else
                lines(sec[,4], sec[,1], lwd=2, col=bestcol)
            if(revaxes)
              abline(v=hiatus, col="grey", lty="dashed") else
                abline(h=hiatus, col="grey", lty="dashed")
          }
      } else
      {
        pol <- cbind(c(calrange[,2], rev(calrange[,3])), c(calrange[,1], rev(calrange[,1])))
        if(plotrange)
          if(revaxes)
            polygon(pol[,2], pol[,1], col=rangecol, border=rangecol) else
              polygon(pol, col=rangecol, border=rangecol)
        if(revaxes)
          lines(calrange[,1], calrange[,4], lwd=2, col=bestcol) else
            lines(calrange[,4], calrange[,1], lwd=2, col=bestcol)
      }

    # draw slumps if these were given
    if(length(slump) > 0)
      for(i in 1:nrow(slump))
        if(revaxes)
          rect(min(slump[i,]), min(yr.lim)-1e4, max(slump[i,]), max(yr.lim)+1e4, col=slumpcol, border=slumpcol) else
            rect(min(yr.lim)-1e4, min(slump[i,]), max(yr.lim)+1e4, max(slump[i,]), col=slumpcol, border=slumpcol)

    # draw the calibrated distributions of the dates
    top <- 1
    for(i in 1:length(dat$depth))
      top <- min(top, max(dat$calib[[i]][,2])) # find the lowest peak

    if(calhght > 0)
      for(i in 1:length(dat$depth))
        {
          if(is.na(dat$cal[[i]])) col <- C14col else col <- calcol
          pol <- dat$calib[[i]] # already normalised to 1
          if(ash) pol[,2] <- pol[,2]/max(pol[,2])/1e3 # draw all same height
          pol[pol[,2] > maxhght,2] <- maxhght
          pol[,2] <- calhght*(dmax-dmin)*pol[,2]/(top*100)
          pol <- cbind(c(pol[,1], rev(pol[,1])),
            c(dat$depth[[i]]-pol[,2], dat$depth[[i]]+mirror*rev(pol[,2])))
          if(revaxes) polygon(pol[,2], pol[,1], col=col, border=col) else
            polygon(pol, col=col, border=col)
        }

    # draw the calibrated ranges of the dates
    for(i in 1:length(dat$depth))
      {
        if(is.na(dat$cal[[i]])) col <- C14col else col <- calcol
        for(j in 1:nrow(dat$hpd[[i]]))
          if(revaxes)
            rect(dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,1], dat$depth[i]+dat$thick[i]/2, dat$hpd[[i]][j,2], lwd=1, lend=2, col=col, border=NA) else
              rect(dat$hpd[[i]][j,1], dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,2], dat$depth[i]+dat$thick[i]/2, lwd=1, lend=2, col=col, border=NA)
      }
    if(length(outliers)>0) # any outliers?
      {
        for(i in outliers)
          for(j in 1:nrow(dat$hpd[[i]]))
            if(revaxes)
              rect(dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,1], dat$depth[i]+dat$thick[i]/2, dat$hpd[[i]][j,2], col=outcol, border=outcol, lwd=1, lend=2) else
                rect(dat$hpd[[i]][j,1], dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,2], dat$depth[i]+dat$thick[i]/2, col=outcol, border=outcol, lwd=1, lend=2)
        if(revaxes)
          points(dat$depth[outliers], dat$mid1[outliers], cex=outlsize, pch=4, col=outcol) else
            points(dat$mid1[outliers], dat$depth[outliers], cex=outlsize, pch=4, col=outcol)
      }
  }


#' @name deptime.depth
#' @title  Calculates *for each iteration* the slope of a straight curve between depths
#'  just above and below the desired point.
#' @description Calculates *for each iteration* the slope of a straight curve between depths
#'  above and below the desired point. Requires sufficiently dense density of depths, e.g. \code{yrsteps=1}.
#' @details 
#' To calculate sedimentation times at a depth. Before running this, run your core in clam and store the data, 
#' so, make sure to set \code{storedat=TRUE}. 
#' Renamed from previous accrate.depth function to avoid confusion with accrate.depth function of rbacon.
#' @param depth The depth for which accumulation rate estimates should be calculated.
#' @param yrcm Calculate in years per cm, or alternatively in cm per yr.
#' @param prob  Probability level at which to calculate the ranges.
#' @author Maarten Blaauw
#' @return The slope of a straight curve between depths above and below the desired point.
#' @examples
#'   clam(coredir=tempdir(), storedat=TRUE) 
#'   deptime.depth(20)
#'   deptime.depth(20, FALSE) # to calculate accumulation rates in cm/yr
#' 
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/clam.html}
#' @references
#' Blaauw, M., 2010. Methods and code for 'classical' age-modelling of radiocarbon sequences. Quaternary Geochronology 5, 512-518
#' \url{http://dx.doi.org/10.1016/j.quageo.2010.01.002}
#' @export
deptime.depth <- function(depth, yrcm=TRUE, prob=.95)
  {
  chron <- get('chron')
  calrange <- get('calrange')
    if(depth <= min(calrange[,1]) || depth >= max(calrange))
      stop("Deposition times cannot be calculated for the top or bottom of the core. Please check the manual", call.=FALSE)
    d <- max(which(calrange[,1] <= depth))
    if(yrcm)
      accrate <- (chron[d+1,]-chron[d-1,]) / (calrange[d+1,1]-calrange[d-1,1]) else
      accrate <- (calrange[d+1,1]-calrange[d-1,1]) / (chron[d+1,]-chron[d-1,])
    acc <- density(accrate)
    plot(acc, main="", xlab=if(yrcm) "yr/cm" else "cm/yr")
    abline(h=0)
    o <- order(acc$y, decreasing=TRUE)
    acc <- cbind(acc$x[o], cumsum(acc$y[o])/sum(acc$y))
    acc <- range(acc[acc[,2] <= prob,1])
    rect(acc[1], 0, acc[2], -999, col=grey(.5), border=grey(.5))
    cat(100*prob, "% ranges: ", acc[1], " to ", acc[2], if(yrcm) " yr/cm\n" else " cm/yr\n", sep="")
  }



#' @name deptime.age 
#' @title Calculates the slope of a straight curve at the desired age.
#' @description Calculates *for each iteration* the slope of a straight curve between 
#' depths above and below the desired age. Requires sufficiently dense density of depths, e.g. \code{steps=1}.
#' @details 
#' To calculate deposition times at an age. Before doing this, run your core in clam and store the data, 
#' so, make sure the option \code{storedat=TRUE}.
#' Renamed from previous accrate.age function to avoid confusion with accrate.age function of rbacon.
#' @param age Age to calculate deposition time (years per cm).
#' @param yrcm Calculate in years per cm, or alternatively in cm per yr.
#' @param prob Probability level at which to calculate the ranges.
#' @author Maarten Blaauw
#' @return The slope of a straight curve between depths above and below the desired point
#' @examples 
#'   clam(coredir=tempdir(), storedat=TRUE)
#'   deptime.age (5000)
#' deptime.age(5000, yrcm=FALSE) # to calculate sedimentation times in cm/yr, so accumulation rates
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/clam.html}
#' @references
#' Blaauw, M., 2010. Methods and code for 'classical' age-modelling of radiocarbon sequences. Quaternary Geochronology 5, 512-518
#' \url{http://dx.doi.org/10.1016/j.quageo.2010.01.002}
#' @export
deptime.age <- function(age, yrcm=TRUE, prob=.95)
  {
    chron <- get('chron') 
    calrange <- get('calrange')
	
    accrate <- c()
    for(i in 1:ncol(chron))
      {
        a <- max(which(chron[,i] <= age))
        if(yrcm)
          accrate <- c(accrate, (chron[a+1,i]-chron[a-1,i]) / (calrange[a+1,1]-calrange[a-1,1])) else
          accrate <- c( accrate, (calrange[a+1,1]-calrange[a-1,1]) / (chron[a+1,i]-chron[a-1,i]))
      }
    acc <- density(accrate)
    plot(acc, main="", xlab=if(yrcm) "yr/cm" else "cm/yr")
    abline(h=0)
    o <- order(acc$y, decreasing=TRUE)
    acc <- cbind(acc$x[o], cumsum(acc$y[o])/sum(acc$y))
    acc <- range(acc[acc[,2] <= prob,1])
    rect(acc[1], 0, acc[2], -999, col=grey(.5), border=grey(.5))
    cat(100*prob, "% ranges: ", acc[1], " to ", acc[2], if(yrcm) " yr/cm\n" else " cm/yr\n", sep="")
  }



#' @name plot_proxies 
#' @title Produce a plot of proxy values against calendar age.
#' @description  Produce a plot of proxy values against calendar age.
#' @details 
#' Only works after running clam on the core using \code{proxies=TRUE}. Requires a file containing the core depths as the first column, 
#' and any proxy values on subsequent columns. Values should be separated by comma's. The file should be stored as a .csv file in the core's directory.
#' @param prox  Position of the proxy that should be plotted, e.g. \code{1} for the first proxy in the file.
#' @param errors Plot an error envelope.
#' @param proxcol Colour of the error envelope.
#' @param revyr Direction of the calendar scale (\code{revyr=TRUE} will reverse the calendar scale from the default \code{FALSE}).
#' @author Maarten Blaauw
#' @return A plot of the age model function with proxies.
#' @examples 
#' clam(coredir=tempdir(), proxies=TRUE)
#' plot_proxies(3)
#' plot_proxies(3, revyr=FALSE)
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/clam.html}
#' @references
#' Blaauw, M., 2010. Methods and code for 'classical' age-modelling of radiocarbon sequences. Quaternary Geochronology 5, 512-518
#' \url{http://dx.doi.org/10.1016/j.quageo.2010.01.002}
#' @export
#' 
plot_proxies <- function(prox, errors=TRUE, proxcol=grey(0.5), revyr=TRUE)
{
  dat <- get('dat') #JEV warning
  calrange <- get('calrange') #JEV warning
  
  prx <- dat$proxies
  if(length(prox)>1) layout(matrix(1:length(prox), ncol=1))
  for(j in 1:length(prox))
  {
    pr <- prx[which(!is.na(prx[,prox+1])),]
    ages <- array(0, dim=c(nrow(pr),3))
    for(i in 1:nrow(pr))
      ages[i,] <- calrange[which(calrange[,1]==pr[i,1]),c(2,3,4)]
    xlim <- range(ages)
    if(!dat$BCAD) xlim <- rev(xlim)
    if(revyr) xlim <- rev(xlim)
    plot(ages[,3], pr[,prox+1], type="n", xlim=xlim, xlab=ifelse(dat$BCAD, "cal BC/AD", "cal BP"), ylab=names(pr)[prox+1])
    if(errors)
      for(i in 2:nrow(pr))
        polygon(c(ages[(i-1):i,1], ages[i:(i-1),2]), c(pr[c((i-1):i, i:(i-1)),prox+1]), col=proxcol, border=proxcol)
    lines(ages[,3], pr[,prox+1])
  }
  layout(1)
}

#' @name calibrate 
#' @title Calibrate individual 14C dates.
#' @description Calibrate individual 14C dates, plot them and report calibrated ranges.
#' @details 
#' Type \code{calibrate()} to see how a date of 2450 +- 50 14C BP gets calibrated (the calibration curve happens to show
#' a plateau around this 14C age). To calibrate a different date, provide its reported mean and error (1 
#' standard deviation error as reported by the radiocarbon laboratory) as follows: \code{calibrate(mean, error)},
#' e.g., for a date of 130 +- 20 14C BP, type calibrate\code{(cage=130, error=20)} or, shorter, \code{calibrate(130,20)}. 
#' As this date will fall partly beyond the younger extreme of the calibration curve, a warning will be given
#' (similar warnings will be given for too old dates).
#'   
#' In case the date has a reservoir effect or age offset, e.g. of 100 14C years, provide this as follows: 
#' \code{calibrate(130, 20, reservoir=100)}. If you want to include an uncertainty for this offset, provide this as follows,
#' e.g., for an uncertainty of 50yr, \code{calibrate(130,20,reservoir=c(100, 50))}. 
#' The uncertainty for the age offset will then be added to the error (by taking the square root of the sum 
#' of the squared error and the squared offset uncertainty). If the carbon of your sample has mixed marine/terrestrial sources,
#' instead apply the marine offset using \code{mix.calibrationcurves} \link{mix.calibrationcurves}, and calibrate the date using that custom-built curve.
#'
#' If you prefer to work with, e.g., 68 \% as opposed to the default 95 \% confidence intervals, 
#' type: \code{calibrate(130, 20, prob=0.68)} or \code{calibrate(130, 20,, 0.68)} (the commas between the brackets indicate the position of the option;
#' the standard deviation is the fourth option of the \code{calibrate} function). Clam calculates the calibrated distribution 
#' for every single calendar year (\code{yrsteps=1}) within a wide range of the 14C date (default but adaptable \code{times=5}
#' standard deviations or 99.999999 \% of its probability distribution). This range can also be adapted by 
#' changing the option expand (default \code{expand=0.1}). Probabilities below a threshold (default \code{threshold=1e-6}) will be neglected.
#' 
#' By default the northern hemisphere terrestrial calibration curve is used (\code{cc=1, cc1="IntCal13.14C"}). 
#' To use alternative curves, use \code{cc=2} (\code{cc2="Marine13.14C"}), \code{cc=3} (\code{cc3="SHCal13.14C"}), 
#' \code{cc=4} (\code{cc4="mixed.14C"}), or change the file names of \code{cc1, cc2, cc3 or cc4}.
#'  
#' Clam works in cal BP (calendar years before AD 1950) by default, but can work with cal BC/AD through the option \code{BCAD=TRUE}. 
#' 
#' By default the Gaussian distribution is used to calibrate dates. For use of the student-t distribution instead, 
#' provide two sensible values, e.g., \code{calibt=c(3,4)}.
#'
#' Calibrated distributions are usually reduced to their 68\% or 95\% calibrated ranges, taking into account the asymmetric 
#' and multi-peaked shape of these distributions. In clam, this is done by calculating the highest posterior density (hpd) ranges: 
#' \itemize{
#' \item i) the probability distribution (see above) is normalised to 100\%
#' \item ii) the calendar years are ranked according to their probabilities
#' \item iii) those calendar ages with a cumulative sum at or above the desired probability threshold (default 95\%) are retained, and 
#' \item iv) the extremes and probabilities of any sub-ranges within these calendar ages are reported. 
#' }
#' Calibrated ranges at 68\% will obviously result in narrower confidence intervals, and a perceived higher precision, than 95\% ranges. However, given the often
#' asymmetric and multi-modal nature of calibrated distributions, the probability that the 'true' calendar date 
#' lies outside the 1 standard deviation hpd ranges is considerable (c. 32\%). Therefore the use of 95\% calibrated ranges is preferable, 
#' and default in clam. The hpd ranges are calculated at yearly resolution by default (\code{hpdsteps=1}).
#'
#' Negative radiocarbon ages are calibrated with postbomb curves, but the user needs to tell clam which curve to use. 
#' For example, to use the first of the three northern hemisphere curves, provide the option \code{postbomb=1}, 
#' while for southern hemisphere samples, use \code{postbomb=4} or \code{postbomb=5}. Default curves can be changed; 
#' currently they are \code{pb1="postbomb_NH1.14C"}, \code{pb2="postbomb_NH2.14C"}, \code{pb3="postbomb_NH3.14C"}, 
#' \code{pb4="postbomb_SH1-2.14C"} and \code{pb5="postbomb_SH3.14C"}; see \url{http://calib.org/CALIBomb/}. 
#' If no \code{postbomb} option is provided 
#' for negative radiocarbon ages, clam will report an error and refuse to calibrate the date. Given the sub-year resolution of postbomb-curves, 
#' hpd ranges are calculated at high resolution by default (\code{pbsteps=0.01}). Choose alternative values with care as 
#' they may cause unexpected results. 
#' 
#' Generally the calculations are removed from memory after calibration; 
#' if you want to have them stored (say for subsequent manipulations), provide the option \code{storedat=TRUE}.
#'
#' A graph of the calibration is produced by default (\code{graph=TRUE}), and it can be adapted in several ways.
#' The limits of the horizontal (calendar scale) and vertical (14C scale) axes are calculated automatically 
#' but can be changed by providing alternative values for the options \code{yrmin, yrmax, minC14} and \code{maxC14}, respectively.
#' The titles of both axis can be changed by providing alternative titles to \code{xlab} and/or \code{ylab}, and 
#' also the top title can be adapted using title. The heights of the distributions of the 14C and calibrated 
#' ages can be set to alternative values using \code{calheight} (default \code{0.3} which plots the distribution up to 30\% of the height of the entire graph).
#' Parameters for white space around the 
#' graph can be changed (default \code{mar=c(3.5, 2, 2, 1}) for spacing below, to the left, above and to the right respectively), 
#' as can the spacing for the axis labels (\code{mgp=c(2,1,0)}). By default, the axes are connected at the lower left, \code{bty="l"}.
#' Check the R documentation of \code{par()} for more options.
#'   
#' The colours of the 14C date, the calibration curve, the entire distributions, as well as of the highest posterior density (\code{hpd}) 
#' ranges, can be changed by providing an alternative colour in \code{date.col}, \code{cc.col}, \code{dist.col}, and/or \code{sd.col}, respectively.
#' The default colours are transparent grey for the dates probability distributions (\code{dist.col=rgb(0,0,0, 0.3)} and \code{sd.col=rgb(0,0,0, 0.5)};
#' change the last value of rgb for different greyscale values), red for the uncalibrated mean and error bars (\code{date.col="red"}), 
#' and transparent green for the calibration curve (\code{cc.col=rgb(0, 0.5, 0, 0.7)}). R's rgb() function expects values between \code{0} and \code{1}
#' for red, green and blue, respectively, followed by a value for the semi-transparency (also between 0 and 1). Some graphic devices 
#' such as postscript are unable to use transparency; in that case provide different colours or leave the fourth value empty.
#' @param cage Mean of the uncalibrated C-14 age.
#' @param error	Error of the uncalibrated C-14 age.
#' @param reservoir Reservoir age, or reservoir age and age offset.
#' @param prob Probability confidence intervals (between 0 and 1).
#' @param cc Calibration curve for C-14 dates (1, 2, 3, or 4).
#' @param cc1 For northern hemisphere terrestrial C-14 dates.
#' @param cc2 For marine C-14 dates.
#' @param cc3 For southern hemisphere C-14 dates.
#' @param cc4 For mixed marine/terrestrial C-14 dates.
#' @param ccdir Directory where the calibration curves for C-14 dates \code{cc} are located. By default \code{ccdir=""}. 
#' Use \code{ccdir="."} to choose current working directory. Use \code{ccdir="/Curves"} to choose sub-folder \code{/Curves}.
#' @param postbomb Calibration curve for postbomb dates.
#' @param pb1 For Northern hemisphere region 1 postbomb C-14 dates.
#' @param pb2 For Northern hemisphere region 2 postbomb C-14 dates.
#' @param pb3 For Northern hemisphere region 3 postbomb C-14 dates.
#' @param pb4 For Southern hemisphere regions 1-2 postbomb C-14 dates.
#' @param pb5 For Southern hemisphere region 3 postbomb C-14 dates.
#' @param yrsteps Temporal resolution at which C-14 ages are calibrated (in calendar years).
#' @param pbsteps Temporal resolution at which postbomb C-14 ages are calibrated (in calendar years).
#' @param hpdsteps Temporal resolution at which highest posterior density ranges are calibrated (in calendar years).
#' @param calibt Off by default; provide two parameters such as c(3,4).
#' @param yrmin Minimum of calendar axis (default calculated automatically).
#' @param yrmax Maximum of calendar axis (default calculated automatically).
#' @param minC14 Minimum age of the C-14 age axis (default calculated automatically).
#' @param maxC14 Maximum of the C-14 age axis (default calculated automatically).
#' @param times Half-range of calibration curve used to calibrate dates (multiplication factor for the date's errors).
#' @param calheight Maximum height of the C14 and calibrated distributions (as proportion of the invisible secondary axes).
#' @param expand By which ratio should the calendar axis be expanded to fit the calibrated distribution.
#' @param threshold Below which value should probabilities be excluded from calculations.
#' @param graph	Plot a graph of the calibrated date. If set to FALSE, only the hpd ranges will be given.
#' @param storedat Store the dates within the R session after a clam run.
#' @param xlab Alternative names can be provided.
#' @param ylab Alternative names can be provided.
#' @param BCAD Use BC/AD or cal BP scale (default cal BP).
#' @param mar Plot margins (amount of white space along edges of axes 1-4).
#' @param mgp Axis text margins (where should titles, labels and tick marks be plotted).
#' @param bty Draw a box around the graph ("n" for none, and "l", "7", "c", "u", "]" or "o" for correspondingly shaped boxes).
#' @param xaxs Whether or not to extend the limits of the horizontal axis. Defaults to \code{xaxs="i"} which does not extend the limits.
#' @param yaxs Whether or not to extend the limits of the vertical axis. Defaults to \code{yaxs="i"} which does not extend the limits.
#' @param title	Title of the graph. Defaults to the values of the uncalibrated date.
#' @param date.col Colour of the "dot-bar" plot of the C14 date. Defaults to \code{date.col="red"}.
#' @param cc.col Colour of the calibration curve. Defaults to semi-transparent dark green; \code{cc.col=rgb(0,.5,0,0.7)}.
#' @param dist.col Colour of the calibrated distribution.
#' @param sd.col Colour of calibrated range.
#' @author Maarten Blaauw
#' @return A graph of the raw and calibrated C-14 date, and the calibrated ranges.
#' @examples 
#' calibrate()
#' calibrate(130, 20)
#' calibrate(130, 20, reservoir=100)
#' calibrate(130, 20, prob=0.68)
#' calibrate(cage=130, error=20)
#' calibrate(130, 20, reservoir=c(100, 50))
#' 
#' @seealso \url{http://www.chrono.qub.ac.uk/blaauw/clam.html}
#' @references
#' Blaauw, M., 2010. Methods and code for 'classical' age-modelling of radiocarbon sequences. Quaternary Geochronology 5, 512-518
#' \url{http://dx.doi.org/10.1016/j.quageo.2010.01.002}
#' @export
calibrate <- function(cage=2450, error=50, reservoir=0, prob=0.95, cc=1, cc1="IntCal13.14C", cc2="Marine13.14C", cc3="SHCal13.14C", cc4="mixed.14C", ccdir="", postbomb=FALSE, pb1="postbomb_NH1.14C", pb2="postbomb_NH2.14C", pb3="postbomb_NH3.14C", pb4="postbomb_SH1-2.14C", pb5="postbomb_SH3.14C", yrsteps=1, pbsteps=0.01, hpdsteps=1, calibt=FALSE, yrmin=c(), yrmax=c(), minC14=c(), maxC14=c(), times=5, calheight=0.3, expand=0.1, threshold=1e-6, storedat=FALSE, graph=TRUE, xlab=c(), ylab=c(), BCAD=FALSE, mar=c(3.5,3,2,1), mgp=c(1.7,.8,0), bty="l", xaxs="i", yaxs="i", title=c(), date.col="red", cc.col=rgb(0,.5,0,0.7), dist.col=rgb(0,0,0,0.3), sd.col=rgb(0,0,0,0.5)) {
  # set the calibration curve
  ccdir <- .validateDirectoryName(ccdir)
  if(ccdir == "")
    ccdir = paste(system.file("extdata", package=packageName()), "/", sep="")

  # set calibration curve
  if(cc==1) calcurve <- read.table(paste(ccdir, cc1,  sep="")) else
    if(cc==2) calcurve <- read.table(paste(ccdir, cc2,  sep="")) else
      if(cc==3) calcurve <- read.table(paste(ccdir, cc3,  sep="")) else
        if(cc==4) calcurve <- read.table(paste(ccdir, cc4,  sep="")) else
            stop("I do not understand which calibration curve you mean, please check the manual", call.=FALSE)     

  # include postbomb curve if required
  if(cage < 0) {
    pb <- 0
    if(postbomb==FALSE)
      stop("\n  Negative 14C age, should I use a postbomb curve?\n", call.=FALSE)
    if(postbomb==1) pb <- pb1 else
       if(postbomb==2) pb <- pb2 else
         if(postbomb==3) pb <- pb3 else
           if(postbomb==4) pb <- pb4 else
             if(postbomb==5) pb <- pb5 else
               stop("I do not understand which postbomb curve you mean, check the manual", call.=FALSE)
	yrsteps <- min(pbsteps, yrsteps)
    if(length(pb) > 0) {
      pb <- read.table(system.file("extdata", pb, package=packageName()))
      pb.x <- seq(min(pb[,1]), max(pb[,1]), by=yrsteps)
      pb.y <- approx(pb[,1], pb[,2], pb.x)$y
      pb.sd <- approx(pb[,1], pb[,3], pb.x)$y
      calcurve <- cbind(c(pb.x, calcurve[,1]), c(pb.y, calcurve[,2]), c(pb.sd, calcurve[,3]))
    }
    cat("  postbomb date, interpolating to every", pbsteps, "yr.")
  }

  # check whether date lies partly or entirely beyond the calibration curve
  if(length(reservoir) == 2) { # assuming that first value is mean offset, second is error
    error <- sqrt(error^2 + reservoir[2]^2)
    reservoir <- reservoir[1]
  }
  border <- 0
  if(cage-reservoir-error < min(calcurve[,2]+calcurve[,3]))
    if(cage-reservoir+error > min(calcurve[,2]-calcurve[,3]))
      border <- 1 else border <- 2
  if(cage-reservoir+error > max(calcurve[,2]-calcurve[,3]))
    if(cage-reservoir-error < max(calcurve[,2]+calcurve[,3]))
      border <- 1 else border <- 2
  if(border == 1)
    cat("\nDate falls partly beyond calibration curve and will be truncated!")
  if(border == 2)
    stop("\nCannot calibrate dates beyond calibration curve!\n\n")

  # work in BC/AD if needed, and prepare for calculations in f14C
  if(BCAD) {
    theta <- 1950-calcurve[,1]
    ad <- max(which(theta > 0)) # one side of the border between AD and BC
    theta <- c(theta[1:(ad-1)], theta[ad]:theta[ad+2], theta[(ad+3):length(theta)])
   	mu <- approx(1950-calcurve[,1], calcurve[,2], theta)$y
    sigma <- approx(1950-calcurve[,1], calcurve[,3], theta)$y
    theta[theta <= 0] <- theta[theta <= 0] - 1
    calcurve <- cbind(theta, mu, sigma)
  } else 
      theta <- calcurve[,1]
  f.mu <- exp(-calcurve[,2]/8033)
  f.sigma <- exp(-(calcurve[,2]-calcurve[,3])/8033) - f.mu
  f.cage <- exp(-(cage-reservoir)/8033)
  f.error <- f.cage - exp(-(cage-reservoir+error)/8033)

  # calibrate the date and report its highest posterior density (hpd) range
  if(length(xlab) == 0)
  xlab <- ifelse(BCAD, "cal BC/AD", "cal BP")
  calib <- .caldist(f.cage, f.error, theta, f.mu, f.sigma, yrsteps, threshold, calibt, BCAD)
  hpd <- .hpd(calib, prob, hpdsteps, yrsteps)
  colnames(hpd) <- c("yrmin", "yrmax", "prob")
  dat <- list(calib=calib, hpd=hpd)
  if(storedat)
    dat <<- dat
  cat("\nmin\tmax\tprob\n")
  for(i in 1:nrow(hpd)) {
    for(j in 1:3)
	  cat(hpd[i,j], "\t")
    cat("\n")
  }
  cat("\n")

  # produce a graph of the calibrated distribution (default)
  if(graph) {
    ifelse(BCAD,
      xrange <- 1950+c((1+expand)*(min(calib[,1])-1950), (1-expand)*(max(calib[,1])-1950)),
      xrange <- c((1+expand)*max(calib[,1]), (1-expand)*min(calib[,1])))
    if(length(yrmin) > 0)
	  xrange[2] <- yrmin
    if(length(yrmax) > 0)
	  xrange[1] <- yrmax
    ifelse(BCAD,
      cc <- calcurve[max(which(theta >= min(xrange))):min(which(theta <= max(xrange))),],
      cc <- calcurve[min(which(theta >= min(xrange))):max(which(theta <= max(xrange))),])

  # first plot the calibrated distribution, and its hpd ranges
  par(mar=mar, mgp=mgp, bty=bty, xaxs=xaxs, xaxt="s", yaxs=yaxs, yaxt="n", new=FALSE)
  pol <- cbind(c(calib[,1], rev(calib[,1])), c(calib[,2]/max(calib[,2]), rep(0, length=nrow(calib))))
  plot(0, type="n", xlim=xrange, ylim=c(0,1/calheight), xlab="", ylab="")
  polygon(pol, col=dist.col, border=NA)
  for(i in 1:nrow(hpd)) {
    if(hpd[i,1]==hpd[i,2]) {
      probs <- calib[which(calib[,1]==hpd[i,1]),]
      lines(rep(probs[1], 2), c(0, probs[2]/max(calib[,2])), col=grey(.5))
    } else {
        probs <- calib[max(which(calib[,1]<=hpd[i,1])):max(which(calib[,1]<=hpd[i,2])),]
        pol <- cbind(c(probs[,1], rev(probs[,1])), c(probs[,2]/max(calib[,2]), rep(0, length=nrow(probs))))
        polygon(pol, col=sd.col, border=NA)
      }
  }
  lines(calib[,1], calib[,2]/max(calib[,2]))
  abline(h=0)

  # now draw the 14C distribution (normal distribution, on vertical axis)
  par(new=TRUE, yaxt="s", yaxs="r", xaxt="n")
  if(length(cc) == 3)
	cc <- cbind(cc[1], cc[2], cc[3])
  if(reservoir != 0)
    main <- substitute(cage-res %+-% er, list(cage=cage, er=error, res=reservoir)) else
      main <- substitute(cage %+-% er, list(cage=cage, er=error))
  if(length(title)>0)
	main <- title
  if(length(minC14) == 0)
    minC14 <- min(cc[,2]-qnorm(1-(1-prob)/2)*cc[,3], cage-reservoir-qnorm(1-(1-prob)/2)*error)
  if(length(maxC14) == 0)
    maxC14 <- max(cc[,2]+qnorm(1-(1-prob)/2)*cc[,3], cage-reservoir+qnorm(1-(1-prob)/2)*error)
  if(length(ylab) == 0)
    ylab <- expression(paste(""^14, "C BP"))
  plot(0, type="n", xlim=xrange, ylim=c(minC14, maxC14), xlab=xlab, ylab=ylab, main=main)
  if(length(calibt) > 0)
    times <- 5*times
  yage <- (cage-reservoir-times*error):(cage-reservoir+times*error) # must not be on F14C for plot
  if(length(calibt) < 2)
    xage <- dnorm(exp(-yage/8033), f.cage, f.error) else
      xage <- (calibt[2] + ((f.cage-exp(-yage/8033))^2) / (2*(f.error^2))) ^ -(calibt[1]+0.5)
  xage.plot <- xrange[1]-((xrange[1]-xrange[2])*calheight)*xage/max(xage)
  pol <- cbind(c(xage.plot, rep(xrange[1], length(xage))), c(yage, rev(yage)))
  polygon(pol, col=dist.col, border="black")

  # draw the highest posterior density (hpd) range of the 14C date
  xage[which(cumsum(xage)/sum(xage) > 1 - (1-prob)/2)] <- 0
  xage[which(cumsum(xage)/sum(xage) < (1-prob)/2)] <- 0
  xage <- xrange[1]-((xrange[1]-xrange[2])*calheight)*(xage/max(xage))
  pol <- cbind(c(xage, rep(xrange[1], length=length(xage))), c(yage, rev(yage)))
  polygon(pol, col=sd.col, border=FALSE)

  # plot the mid and error of the 14C date
  points(xrange[1]-.01*(xrange[1]-xrange[2]), cage, pch=19, col=date.col)
  lines(rep(xrange[1]-.01*(xrange[1]-xrange[2]), 2), c(cage-error, cage+error), lwd=2, col=date.col)

  # now draw the calibration curve
  pol <- cbind(c(theta, rev(theta)),
    c(calcurve[,2]-qnorm(1-(1-prob)/2)*calcurve[,3], rev(calcurve[,2]+qnorm(1-(1-prob)/2)*calcurve[,3])))
  polygon(pol, border=cc.col, col=cc.col)
  }
}


# list the available cores
clam_runs <- list.files("clam_runs/")

.validateDirectoryName <- function(dir) {
  if(!dir.exists(dir))
    dir.create(dir, showWarnings=FALSE, recursive=TRUE)
  dir <- suppressWarnings(normalizePath(dir))
  lastchar <- substr(dir, nchar(dir), nchar(dir))
  if(lastchar != "/" & lastchar != "\\" & lastchar != "" & lastchar != "." )
    dir <- paste(dir, "/", sep="") # does this work in Windows?
  return(dir)
}

# function to load results into global environment
# parameter position defaults to 1, which equals an assignment to the global environment
.assign_to_global <- function(key, val, pos=1) {
  assign(key, val, envir=as.environment(pos) )
}


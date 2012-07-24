



################################################################################
################################################################################
#' @description This function provides a GUI for the function \link{rrisk.BayesZIP}.
#'
#' @name ZIPGUI
#' @aliases ZIPGUI
#' @title Bayes estimation of a zero inflated Poisson (ZIP) model
#' @usage ZIPGUI(data=NULL, prior.lambda=c(1,10), prior.pi=c(0.8,1),
#'  chains=3, burn=1000, update=10000, thin=1)
#' @param data a vector of numeric data, possibly containing zeros, and of minimal length 10.
#' @param prior.lambda numeric vector containing minimum and maximum of a uniform
#' distribution used as prior for the Poisson parameter \code{lambda}, e.g. \cr \code{lambda} \eqn{\sim} \code{prior.lambda(*,*)=unif(*,*)}
#' @param prior.pi numeric vector containing parameters of a beta distribution
#' describing prior knowledge about prevalence (proportion of contaminated samples), e.g. \cr \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}
#' @param chains positive single numeric value, number of independent MCMC chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period (default 1000)
#' @param update positive single numeric value, length of update iterations for estimation (default 10000)
#' @param thin positive single numeric value (default 1). The samples from every kth iteration will be used for 
#'        inference, where k is the value of thin. Setting \code{thin > 1} can help to reduce the autocorrelation
#'        in the sample.
#' @return The function \code{ZIPGUI} returns an instance of the \code{\linkS4class{bayesmodelClass}}
#' class containing following informations
#' \item{\code{convergence}}{logical, whether the model has converged (assessed by the user)}
#' \item{\code{results}}{data frame containing statitsics of the posterior distribution}
#' \item{\code{jointpost}}{data frame giving the joint posterior probability distribution}
#' \item{\code{nodes}}{names of the parameters jointly estimated by the Bayes model}
#' \item{\code{model}}{model in BRugs/Winbugs syntax as a character string}
#' \item{\code{chains}}{number of independent MCMC chains}
#' \item{\code{burn}}{length of burn-in period}
#' \item{\code{update}}{length of update iterations for estimation}
#' @note The convergence of the model will be entered by the user after the simulation process.
#' @seealso \code{\link{rrisk.BayesZIP}}
#' @keywords manip
#' @export
#' @references Bohning, D., E. Dietz, P. Schlattman, L. Mendonca, and U. Kirchner (1999). 
#' The zero-inflated Poisson model and the decayed, missing and filled teeth index in 
#' dental epidemiology. Journal of the Royal Statistical Society, Series A 162, 195-209.
#' @examples
#' \donttest{
#' data <- rpois(30, 4)
#' prior.lambda <- c(1, 10)
#' prior.pi <- c(0.8, 1)
#' ZIPGUI(data, prior.lambda, prior.pi)}


ZIPGUI <- function(data=NULL, prior.lambda=c(1, 10), prior.pi=c(0.8, 1),
  chains=3, burn=1000, update=10000,thin=1){
  #-----------------------------------------------------------------------------
  # define buttons and their actions
  #-----------------------------------------------------------------------------
  # GUIDiag: function for plotting the two diagnosis plots
  GUIDiag<-function(nodes,plotnumber=1){
    len<-length(nodes)
    if(plotnumber==1){
      par(mfrow=c(2,len))
      for(i in 1:len){ 
        if(len<=3){
          maintext<-paste("Gelman-Rubin conv. \nstatistic for",nodes[i])
        } else {
          maintext<-paste("G.-R. statistic for",nodes[i])
        } # end if
        rrisk.samplesBgr(nodes[i],main=maintext,cex.main=1)
        abline(h=1,lty="dashed")  
      } # end first for-loop
      for(i in 1:len){
        plotDensity(nodes[i],main=paste("Density plot of ",nodes[i]),cex.main=1)
      } # end second for-loop
    } else if (plotnumber==2){
      if(len<=3){
        par(mfrow=c(len,1))
      } else {
        par(mfrow=c(2,3))
      } #enf if-else
      for(i in 1:len){
        plotHistory(nodes[i],main=paste("Convergence monitored for node",nodes[i]),cex.main=1)
      } # emd for-loop  
    } else if(plotnumber==3){
      if(len<=3){
        par(mfrow=c(len,1))
      } else {
        par(mfrow=c(2,3))
      } #enf if-else
      for(i in 1:len){
        plotAutoC(nodes[i],main=paste("Autocorrelation plot for node",nodes[i]),cex.main=0.5)
      } # emd for-loop  
    }# end overall if
  } # end function GUIDiag()

  #-----------------------------------------------------------------------------
  # what happens by pressing RUN button
  #-----------------------------------------------------------------------------
  onRun <- function(...){
    prior.pi.l. <- as.numeric(tclvalue(tkget(prior.pi.l)))
    prior.pi.r. <- as.numeric(tclvalue(tkget(prior.pi.r)))
    prior.lambda.l. <- as.numeric(tclvalue(tkget(prior.lambda.l)))
    prior.lambda.r. <- as.numeric(tclvalue(tkget(prior.lambda.r)))
    chains. <- as.numeric(tclvalue(tkget(chainsEntry)))
    burn. <- as.numeric(tclvalue(tkget(burnEntry)))
    update. <- as.numeric(tclvalue(tkget(updateEntry)))
    thin. <- as.numeric(tclvalue(tkget(thinEntry)))
 
    ##try.result <- try(
    mod <- rrisk.BayesZIP(data=data, simulation=TRUE, prior.pi=c(prior.pi.l., prior.pi.r.), 
                          prior.lambda=c(prior.lambda.l., prior.lambda.r.),
                          chains=chains., burn=burn., update=update.,thin=thin.)
    
    if(!is(mod,"bayesmodelClass")){
      assign("mod", value=mod, envir=envirZIP)
      tkdestroy(bayesZIPWindow)
      stop("Any error occured during fitting bayesian ZIP model",call.=FALSE)
    } else {
      ##,silent=TRUE)
      ##if (inherits(try.result, "try-error")) { cat("Error: The simulation could not be completed!\n")}
      assign("mod", value=mod, envir=envirZIP)
  
      if(exists("imgPlot1", where=envirZIP)){
        imgPlot1 <- get("imgPlot1", envir=envirZIP)
        tkpack.forget(imgPlot1)
      }
      if(exists("imgPlot2", where=envirZIP)){
        imgPlot2 <- get("imgPlot2", envir=envirZIP)  
        tkpack.forget(imgPlot2)
      } 
      if(exists("imgPlot3", where=envirZIP)){
        imgPlot2 <- get("imgPlot3", envir=envirZIP)  
        tkpack.forget(imgPlot3)
      }  
    
      nodes <- c("pi", "lambda")
      imgPlot1 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=1),hscale=1.4,vscale=1.4)
      imgPlot2 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=2),hscale=1.4,vscale=1.4)
      imgPlot3 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=3),hscale=1.4,vscale=1.4)
   
      assign("imgPlot1", value=imgPlot1, envir=envirZIP)
      assign("imgPlot2", value=imgPlot2, envir=envirZIP)
      assign("imgPlot3", value=imgPlot3, envir=envirZIP)
  
      tkpack(imgPlot1, side="top")
      assign("plotNumber",value=1,envir=envirZIP)
   
      tkraise(bayesZIPWindow)
    }
  } # end onRun() fucntion
  
  #-----------------------------------------------------------------------------
  # what happens by pressing RESET button
  #-----------------------------------------------------------------------------
  onReset <- function(...){
    tkconfigure(prior.pi.l, text=tclVar(prior.pi[1]))
    tkconfigure(prior.pi.r, text=tclVar(prior.pi[2]))
    tkconfigure(prior.lambda.l, text=tclVar(prior.lambda[1]))
    tkconfigure(prior.lambda.r, text=tclVar(prior.lambda[2]))
    tkconfigure(chainsEntry, text=tclVar(chains))
    tkconfigure(burnEntry, text=tclVar(burn))
    tkconfigure(updateEntry, text=tclVar(update))
    tkconfigure(thinEntry, text=tclVar(thin))
  } # end onReset() function

  #-----------------------------------------------------------------------------
  # what happens by pressing CANCEL button
  #-----------------------------------------------------------------------------
  onCancel <- function(...){
    tkdestroy(bayesZIPWindow) 
  } # end onCancel() function
  
  #-----------------------------------------------------------------------------
  # what happens by pressing Nexplot button
  #-----------------------------------------------------------------------------
  onNextplot <- function(...){
  # get the current plot number (in 1,2) and switch to the other plot
    if(!exists("plotNumber", envirZIP)) return(FALSE) # if there was no model run yet
    
    plotNumber <- get("plotNumber", envir=envirZIP)
    if(plotNumber==1){
      imgPlot1 <- get("imgPlot1", envir=envirZIP)
      tkpack.forget(imgPlot1)
      assign("plotNumber",value=2,envir=envirZIP)
      imgPlot2 <- get("imgPlot2", envir=envirZIP)
      tkpack(imgPlot2, side="top")
    }  else if(plotNumber==2){       
      imgPlot2 <- get("imgPlot2", envir=envirZIP)
      tkpack.forget(imgPlot2)
      assign("plotNumber",value=3,envir=envirZIP)
      imgPlot3 <- get("imgPlot3", envir=envirZIP)
      tkpack(imgPlot3, side="top")
    } else if(plotNumber==3){       
      imgPlot3 <- get("imgPlot3", envir=envirZIP)
      tkpack.forget(imgPlot3)
      assign("plotNumber",value=1,envir=envirZIP)
      imgPlot1 <- get("imgPlot1", envir=envirZIP)
      tkpack(imgPlot1, side="top")
    }
    tkraise(bayesZIPWindow)  
  } # end onNextplot() function

  #-----------------------------------------------------------------------------
  # what happens by pressing CONVERGENCE button
  #-----------------------------------------------------------------------------
  onConv <- function(...){
    mod <- get("mod", envir=envirZIP)
    mod@convergence <- TRUE
    assign("mod", value=mod, envir=envirZIP)
    tkdestroy(bayesZIPWindow)
  } # end fucntion onConv()

  #-----------------------------------------------------------------------------
  # what happens by pressing CONVERGENCE button
  #-----------------------------------------------------------------------------
  onNotconv <- function(...){
   mod <- get("mod", envir=envirZIP)
    mod@convergence <- FALSE
    assign("mod", value=mod, envir=envirZIP)
    tkdestroy(bayesZIPWindow)
  } # end fucntion onNotconv()
  
  #-----------------------------------------------------------------------------
  # define help varriable(s)
  #-----------------------------------------------------------------------------
  assign("envirZIP", value=new.env(), envir=.GlobalEnv)

  #-----------------------------------------------------------------------------
  # define GUI window and frames
  #-----------------------------------------------------------------------------
  bayesZIPWindow <- tktoplevel()
  tkwm.title(bayesZIPWindow, "GUI for the function rrisk.bayesZIP")
  tkwm.resizable(bayesZIPWindow, 0, 0) # fixed size
  tkwm.maxsize(bayesZIPWindow,800,500)
  tkwm.minsize(bayesZIPWindow,800,500 )

  mod <- new("bayesmodelClass")

  leftFrame <- tkframe(bayesZIPWindow)
  rightFrame <- tkframe(bayesZIPWindow)
  imgFrame <- tkframe(rightFrame,height=500,width=200)
  inputFrame <- tkframe(leftFrame)
  lButtonFrame <- tkframe(leftFrame)
  rButtonFrame <- tkframe(leftFrame)
  #rButtonFrame <- tkframe(rightFrame)

  #-----------------------------------------------------------------------------
  # define input fields
  #-----------------------------------------------------------------------------
  prior.lambda.l <- tkentry(inputFrame, text=tclVar(prior.lambda[1]), width=6)
  prior.lambda.lLabel <- tklabel(inputFrame, text="prior.lambda (unif)")
  prior.lambda.r <- tkentry(inputFrame, text=tclVar(prior.lambda[2]), width=6)

  prior.pi.l <- tkentry(inputFrame, text=tclVar(prior.pi[1]), width=6)
  prior.pi.lLabel <- tklabel(inputFrame, text="prior.pi (beta)")
  prior.pi.r <- tkentry(inputFrame, text=tclVar(prior.pi[2]), width=6)

  chainsEntry <- tkentry(inputFrame, text=tclVar(chains))
  chainsLabel <- tklabel(inputFrame, text="chains")

  burnEntry <- tkentry(inputFrame, text=tclVar(burn))
  burnLabel <- tklabel(inputFrame, text="burn")

  updateEntry <- tkentry(inputFrame, text=tclVar(update))
  updateLabel <- tklabel(inputFrame, text="update")
  
  thinEntry <- tkentry(inputFrame, text=tclVar(thin))
  thinLabel <- tklabel(inputFrame, text="thin")
  
  #-----------------------------------------------------------------------------
  # define buttons
  #-----------------------------------------------------------------------------
  runButton <- ttkbutton(lButtonFrame, width=12, text="Run", command=onRun)
  resetButton <- ttkbutton(lButtonFrame, width=12, text="Reset", command=onReset)
  cancelButton <- ttkbutton(lButtonFrame, width=12, text="Cancel", command=onCancel)

  nextplotButton <- ttkbutton(rButtonFrame, width=12, text="Next Plot", command=onNextplot)
  convButton <- ttkbutton(rButtonFrame, width=12, text="Converge", command=onConv)
  notconvButton <- ttkbutton(rButtonFrame, width=12, text="Not Converge", command=onNotconv)

  #-----------------------------------------------------------------------------
  ## tkgrid() and tkpack() the inputs and frames together
  #-----------------------------------------------------------------------------
  tkgrid(prior.lambda.lLabel, prior.lambda.l, prior.lambda.r, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=2)
  tkgrid(prior.pi.lLabel, prior.pi.l, prior.pi.r, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=2)
  tkgrid(chainsLabel, chainsEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(burnLabel, burnEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(updateLabel, updateEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(thinLabel, thinEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(runButton, resetButton, cancelButton, sticky="we", padx=c(5,5))
  tkgrid(nextplotButton, convButton, notconvButton, sticky="swe", padx=c(5,5))

  tkpack(inputFrame, side="top")
  tkpack(rButtonFrame, pady=c(5,25), side="bottom") 
  tkpack(lButtonFrame, side="bottom", pady=c(0,25))                
  tkpack(leftFrame, side="left") 
  
  tkpack(imgFrame, side="top", padx=c(15,0), pady=c(0,10))
  tkpack(rightFrame, side="right")
                   
  tkraise(bayesZIPWindow)

  tkwait.window(bayesZIPWindow) # otherwise, mod@convergence won't be saved

  #-----------------------------------------------------------------------------
  # prepare output
  #-----------------------------------------------------------------------------
  if(exists("mod", where=envirZIP)){
    mod <- get("mod", envir=envirZIP) 
  } # end if
  return(mod)
} # end of function ZIPGUI



################################################################################
################################################################################
#' @description This function provides a GUI for the function \link{rrisk.BayesPEM}.
#' 
#' @name PEMGUI
#' @aliases PEMGUI
#' @title GUI for Bayesian Prevalence estimation under misclassification (PEM)
#' @usage PEMGUI(x=20, n=20, k=10, prior.pi=c(1,19), prior.se=c(1,1),
#'  prior.sp=c(1,1), chains=3, burn=1000, update=10000, thin=1)
#' @param x scalar value for number of pools (\code{k>1}) or single outcomes (\code{k=1}) with positive test result
#' @param n scalar value for number of pools tested (\code{k>1}) or the sample size in application study (\code{k=1})
#' @param k scalar value for number of individual samples physically combined into one pool
#' @param prior.pi numeric vector containing parameters of a beta distribution as prior for prevalence \code{pi}, e.g. \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}
#' @param prior.se numeric vector containing parameters of a beta distribution as prior for sensitivity \code{se}, e.g. \code{se} \eqn{\sim} \code{prior.se(*,*)=beta(*,*)}
#' @param prior.sp numeric vector containing parameters of a beta distribution as prior for specificity \code{sp}, e.g. \code{sp} \eqn{\sim} \code{prior.sp(*,*)=beta(*,*)}
#' @param chains positive single numeric value, number of independent MCMC chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period (default 1000)
#' @param update positive single numeric value, length of update iterations for estimation (default 10000)
#' @param thin positive single numeric value (default 1). The samples from every kth iteration will be used for 
#'        inference, where k is the value of thin. Setting \code{thin > 1} can help to reduce the autocorrelation
#'        in the sample.
#' @return The function \code{PEMGUI} returns an instance of the \code{\linkS4class{bayesmodelClass}}
#' class containing following informations
#' \item{\code{convergence}}{logical, whether the model has converged (assessed by the user)}
#' \item{\code{results}}{data frame containing statitsics of the posterior distribution}
#' \item{\code{jointpost}}{data frame giving the joint posterior probability distribution}
#' \item{\code{nodes}}{names of the parameters jointly estimated by the Bayes model}
#' \item{\code{model}}{model in BRugs/Winbugs syntax as a character string}
#' \item{\code{chains}}{number of independent MCMC chains}
#' \item{\code{burn}}{length of burn-in period}
#' \item{\code{update}}{length of update iterations for estimation}
#' @note The convergence of the model is assessed by the user using diagnostic plots
#' provided by the \pkg{BRugs} package.
#' @seealso \code{\link{rrisk.BayesPEM}}
#' @keywords manip
#' @export
# @references Greiner, M., Belgoroski, N., NN (2011). Estimating prevalence using diagnostic data 
# with pooled samples: R function \code{rrisk.BayesPEM}. J.Stat.Software (in preparation). 
#' @examples
#' #------------------------------------------
#' # Example of PEM model. Without parameters,
#' # the input fields will show default values
#' #------------------------------------------
#' \donttest{mod <- PEMGUI()}

PEMGUI <- function(x=20, n=20, k=10, prior.pi=c(1,19), prior.se=c(1,1),
                    prior.sp=c(1,1), chains=3,burn=1000, update=10000, thin=1){
  #-----------------------------------------------------------------------------
  # GUIDiag: function for plotting the two diagnosis plots
  #-----------------------------------------------------------------------------
  GUIDiag<-function(nodes,plotnumber=1)
    { len<-length(nodes)
      if(plotnumber==1){ par(mfrow=c(2,len))
        for(i in 1:len){ 
          if(len<=3){
            maintext<-paste("Gelman-Rubin conv. \nstatistic for",nodes[i])
          } else maintext<-paste("G.-R. statistic for",nodes[i])
          rrisk.samplesBgr(nodes[i],main=maintext,cex.main=1)
          abline(h=1,lty="dashed")  
        }
        for(i in 1:len){
          plotDensity(nodes[i],main=paste("Density plot of ",nodes[i]),cex.main=1)
        }
      } else if (plotnumber==2){
        if(len<=3){
          par(mfrow=c(len,1))
        } else  par(mfrow=c(2,3))
        for(i in 1:len){
          plotHistory(nodes[i],main=paste("Convergence monitored for node",nodes[i]),cex.main=1)
        }
      } else if(plotnumber==3){
        if(len<=3){
          par(mfrow=c(len,1))
        } else  par(mfrow=c(2,3))
        for(i in 1:len){
          plotAutoC(nodes[i],main=paste("Autocorrelation plot for node",nodes[i]),cex.main=1)
        }
      }
    }
  #-----------------------------------------------------------------------------
  # what happens by pressing RUN button
  #-----------------------------------------------------------------------------
  onRun <- function(...){
    x. <- as.numeric(tclvalue(tkget(xEntry)))
    n. <- as.numeric(tclvalue(tkget(nEntry)))
    k. <- as.numeric(tclvalue(tkget(kEntry)))
    prior.pi.l. <- as.numeric(tclvalue(tkget(prior.pi.l)))
    prior.pi.r. <- as.numeric(tclvalue(tkget(prior.pi.r)))
    prior.se.l. <- as.numeric(tclvalue(tkget(prior.se.l)))
    prior.se.r. <- as.numeric(tclvalue(tkget(prior.se.r)))
    prior.sp.l. <- as.numeric(tclvalue(tkget(prior.sp.l)))
    prior.sp.r. <- as.numeric(tclvalue(tkget(prior.sp.r)))
    chains. <- as.numeric(tclvalue(tkget(chainsEntry)))
    burn. <- as.numeric(tclvalue(tkget(burnEntry)))
    update. <-  as.numeric(tclvalue(tkget(updateEntry)))
    thin. <-  as.numeric(tclvalue(tkget(thinEntry)))
    #misclass. <- misclasses[as.numeric(tclvalue(tcl(misclassEntry, "getvalue")))+1]
  
    # testmodel:
    # rrisk.BayesPEM(20, 20, 10, simulation=FALSE, prior.pi=c(0.1,0.2), prior.se=c(0.1,0.2), prior.sp=c(0.1,0.2), chains=10, misclass="pool")
  
    mod <- rrisk.BayesPEM(x=x., n=n., k=k., simulation=TRUE, prior.pi=c(prior.pi.l., prior.pi.r.), 
                          prior.se=c(prior.se.l., prior.se.r.), prior.sp=c(prior.sp.l., prior.sp.r.),
                          chains=chains., misclass="pool",update=update., burn=burn.,thin=thin.)
  
    assign("mod", value=mod, envir=envirPEM)
  
    if(exists("imgPlot1", where=envirPEM)){
      imgPlot1 <- get("imgPlot1", envir=envirPEM)
      tkpack.forget(imgPlot1)
    }
    if(exists("imgPlot2", where=envirPEM)){
      imgPlot2 <- get("imgPlot2", envir=envirPEM)
      tkpack.forget(imgPlot2)
    } 
    if(exists("imgPlot3", where=envirPEM)){
      imgPlot3 <- get("imgPlot3", envir=envirPEM)
      tkpack.forget(imgPlot3)
    }   
  
    nodes <- c("pi", "se", "sp")
    imgPlot1 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=1),hscale=1.65,vscale=1.5)
    imgPlot2 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=2),hscale=1.55,vscale=1.7)
    imgPlot3 <- tkrplot(imgFrame,fun=function()GUIDiag(nodes,plotnumber=3),hscale=1.55,vscale=1.7)
   
    assign("imgPlot1", value=imgPlot1, envir=envirPEM)
    assign("imgPlot2", value=imgPlot2, envir=envirPEM)
    assign("imgPlot3", value=imgPlot3, envir=envirPEM)
  
    tkpack(imgPlot1, side="top")
    assign("plotNumber",value=1,envir=envirPEM)
   
    tkraise(bayesPEMWindow)
  }
  
  #-----------------------------------------------------------------------------
  # what happens by pressing RESET button
  #-----------------------------------------------------------------------------
  onReset <- function(...){
    tkconfigure(xEntry, text=tclVar(x))
    tkconfigure(nEntry, text=tclVar(n))
    tkconfigure(kEntry, text=tclVar(k))
    tkconfigure(prior.pi.l, text=tclVar(prior.pi[1]))
    tkconfigure(prior.pi.r, text=tclVar(prior.pi[2]))
    tkconfigure(prior.se.l, text=tclVar(prior.se[1]))
    tkconfigure(prior.se.r, text=tclVar(prior.se[2]))
    tkconfigure(prior.sp.l, text=tclVar(prior.sp[1]))
    tkconfigure(prior.sp.r, text=tclVar(prior.sp[2]))
    #tkconfigure(misclassEntry, text=misclass)
    tkconfigure(chainsEntry, text=tclVar(chains))
    tkconfigure(burnEntry, text=tclVar(burn))
    tkconfigure(updateEntry, text=tclVar(update))
    tkconfigure(thinEntry, text=tclVar(thin))
  }
  
  #-----------------------------------------------------------------------------
  # what happens by pressing CANCEL button
  #-----------------------------------------------------------------------------
  onCancel <- function(...){
    tkdestroy(bayesPEMWindow)
  }
  
  #-----------------------------------------------------------------------------
  # what happens by pressing NEXT button
  #-----------------------------------------------------------------------------
  onNextplot <- function(...){
  # get the current plot number (in 1,2) and switch to the other plot
    if(!exists("plotNumber", envirPEM)) return(FALSE) # if there was no model run yet
    
    plotNumber <- get("plotNumber", envir=envirPEM)
    if(plotNumber==1){
      imgPlot1 <- get("imgPlot1", envir=envirPEM)
      tkpack.forget(imgPlot1)
      assign("plotNumber",value=2,envir=envirPEM)
      imgPlot2 <- get("imgPlot2", envir=envirPEM)
      tkpack(imgPlot2, side="top")
    } else if(plotNumber==2){       
      imgPlot2 <- get("imgPlot2", envir=envirPEM)
      tkpack.forget(imgPlot2)
      assign("plotNumber",value=3,envir=envirPEM)
      imgPlot3 <- get("imgPlot3", envir=envirPEM)
      tkpack(imgPlot3, side="top")
    } else if(plotNumber==3){       
      imgPlot3 <- get("imgPlot3", envir=envirPEM)
      tkpack.forget(imgPlot3)
      assign("plotNumber",value=1,envir=envirPEM)
      imgPlot1 <- get("imgPlot1", envir=envirPEM)
      tkpack(imgPlot1, side="top")
    }
    tkraise(bayesPEMWindow)  
  }
  
  #-----------------------------------------------------------------------------
  # what happens by pressing CONVERGE button
  #-----------------------------------------------------------------------------
  onConv <- function(...){
    mod <- get("mod", envir=envirPEM)
    mod@convergence <- TRUE
    assign("mod", value=mod, envir=envirPEM)
    tkdestroy(bayesPEMWindow)
  }
  
  #-----------------------------------------------------------------------------
  # what happens by pressing NOT CONVERGE button
  #-----------------------------------------------------------------------------
  onNotconv <- function(...){
    mod <- get("mod", envir=envirPEM)
    mod@convergence <- FALSE
    assign("mod", value=mod, envir=envirPEM)
    tkdestroy(bayesPEMWindow)
  }
  
  #-----------------------------------------------------------------------------
  # define Dialof window
  #----------------------------------------------------------------------------- 
  tclRequire("BWidget")
  
  assign("envirPEM",value=new.env(),envir=.GlobalEnv)
  
  mod <- new("bayesmodelClass")
  
  #-----------------------------------------------------------------------------
  ## define GUI window and frames
  #-----------------------------------------------------------------------------
  bayesPEMWindow <- tktoplevel()
  tkwm.title(bayesPEMWindow, "GUI for the function rrisk.bayesPEM")
  tkwm.resizable(bayesPEMWindow, 0, 0) # fixed size
  tkwm.maxsize(bayesPEMWindow,900,600)
  tkwm.minsize(bayesPEMWindow,900,600)
  
  leftFrame <- tkframe(bayesPEMWindow)
  rightFrame <- tkframe(bayesPEMWindow)
  imgFrame <- tkframe(rightFrame, height=600,width=200)
  inputFrame <- tkframe(leftFrame)
  lButtonFrame <- tkframe(leftFrame)
  rButtonFrame <- tkframe(leftFrame)
  
  #-----------------------------------------------------------------------------
  ## define input fields
  #-----------------------------------------------------------------------------
  xEntry <- tkentry(inputFrame, text=tclVar(x))
  xLabel <- tklabel(inputFrame, text="x")
  
  nEntry <- tkentry(inputFrame, text=tclVar(n))
  nLabel <- tklabel(inputFrame, text="n")
  
  kEntry <- tkentry(inputFrame, text=tclVar(k))
  kLabel <- tklabel(inputFrame, text="k")
  
  prior.pi.l <- tkentry(inputFrame, text=tclVar(prior.pi[1]), width=6)
  prior.pi.lLabel <- tklabel(inputFrame, text="prior.pi (beta)")
  prior.pi.r <- tkentry(inputFrame, text=tclVar(prior.pi[2]), width=6)
  
  prior.se.l <- tkentry(inputFrame, text=tclVar(prior.se[1]), width=6)
  prior.se.lLabel <- tklabel(inputFrame, text="prior.se (beta)")
  prior.se.r <- tkentry(inputFrame, text=tclVar(prior.se[2]), width=6)
  
  prior.sp.l <- tkentry(inputFrame, text=tclVar(prior.sp[1]), width=6)
  prior.sp.lLabel <- tklabel(inputFrame, text="prior.sp (beta)")
  prior.sp.r <- tkentry(inputFrame, text=tclVar(prior.sp[2]), width=6)
  
  #misclassLabel <- tklabel(inputFrame, text="misclass")
  #misclasses <- c("pool", "individual", "compare")
  #misclassEntry <- tkwidget(inputFrame, "ComboBox", editable=FALSE, values=misclasses, text=misclass)
  
  chainsEntry <- tkentry(inputFrame, text=tclVar(chains))
  chainsLabel <- tklabel(inputFrame, text="chains")
  
  burnEntry <- tkentry(inputFrame, text=tclVar(burn))
  burnLabel <- tklabel(inputFrame, text="burn")
                                                   
  updateEntry <- tkentry(inputFrame, text=tclVar(update))
  updateLabel <- tklabel(inputFrame, text="update")  
  
  thinEntry <- tkentry(inputFrame, text=tclVar(thin))
  thinLabel <- tklabel(inputFrame, text="thin")  
  
  runButton <- ttkbutton(lButtonFrame, width=12, text="Run", command=onRun)
  resetButton <- ttkbutton(lButtonFrame, width=12, text="Reset", command=onReset)
  cancelButton <- ttkbutton(lButtonFrame, width=12, text="Cancel", command=onCancel)
  
  nextplotButton <- ttkbutton(rButtonFrame, width=12, text="Next Plot", command=onNextplot)
  convButton <- ttkbutton(rButtonFrame, width=12, text="Converge", command=onConv)
  notconvButton <- ttkbutton(rButtonFrame, width=12, text="Not Converge", command=onNotconv)
  
  
  ## tkgrid() and tkpack() the inputs and frames together
  ## must be run "inside-out"
  ################################################################
  
  tkgrid(xLabel, xEntry, sticky="nw", padx=c(10,10), pady=c(10,15), columnspan=3)
  tkgrid(nLabel, nEntry, sticky="nw", padx=c(10,10), pady=c(0,15), columnspan=3)
  tkgrid(kLabel, kEntry, sticky="nw", padx=c(10,10), pady=c(0,15), columnspan=3)
  tkgrid(prior.pi.lLabel, prior.pi.l, prior.pi.r, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=2)
  tkgrid(prior.se.lLabel, prior.se.l, prior.se.r, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=2)
  tkgrid(prior.sp.lLabel, prior.sp.l, prior.sp.r, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=2)
  #tkgrid(misclassLabel, misclassEntry, columnspan=3, padx=c(10,0), pady=c(0,15))
  tkgrid(chainsLabel, chainsEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(burnLabel, burnEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(updateLabel, updateEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
  tkgrid(thinLabel, thinEntry, sticky="nw", padx=c(10,0), pady=c(0,15), columnspan=3)
                   
  tkgrid(runButton, resetButton, cancelButton, sticky="we", padx=c(5,5))
  tkgrid(nextplotButton, convButton, notconvButton, sticky="swe", padx=c(5,5))
  
  tkpack(inputFrame, side="top")
  tkpack(rButtonFrame, pady=c(5,25), side="bottom") 
  tkpack(lButtonFrame, side="bottom", pady=c(0,25))                
  tkpack(leftFrame, side="left") 
  
  tkpack(imgFrame, side="top", padx=c(15,0), pady=c(0,10))
  tkpack(rightFrame, side="right") # zum zweiten mal
  
  tkraise(bayesPEMWindow)
  
  tkwait.window(bayesPEMWindow) # otherwise, mod@convergence won't be saved
  
  if(exists("mod", where=envirPEM)){
    mod <- get("mod", envir=envirPEM)
  } # end if statement 
 
  return(mod)
} # end of function PEMGUI




################################################################################  
################################################################################
#' @description S4 class for displaying the output of the functions \code{rrisk.BayesRGE} and 
#' \code{rrisk.BayesZIP}.
#'
#' @name bayesmodelClass-class
#' @aliases bayesmodelClass
#' @docType class
#' @title S4 class for displaying output of Bayesian models
#' @slot convergence
#' @slot results
#' @slot jointposterior
#' @slot nodes
#' @slot model
#' @slot chains
#' @slot burn
#' @slot updat
#' @rdname bayesmodelClass-class
#' @exportClass bayesmodelClass
#' @examples
#' new("bayesmodelClass")
#' new("bayesmodelClass",model="Some model...", nodes="Nodes info...")

setClass(Class="bayesmodelClass",
  representation=representation(
    convergence="logical",
    results="ANY",
    jointpost="ANY",
    nodes="ANY",
    model="ANY",
    chains="ANY",
    burn="ANY",
    update="ANY"),
  prototype=prototype(
    convergence=FALSE,
    results=NULL,
    jointpost=NULL,
    nodes=NULL,
    model=NULL,
    chains=NULL,
    burn=NULL,
    update=NULL))
  
################################################################################
#' @description Show method for \code{\linkS4class{bayesmodelClass}}
#'
#' @name show-methods
#' @aliases show,bayesmodelClass-method
#' @docType methods
#' @title Show method for bayesmodelClass
#' @param object a \code{bayesmodelClass} object
#' @exportMethod show
#' @importFrom methods show
#' @examples
#' show(new("bayesmodelClass"))
 
 setMethod(f="show",
  signature=signature(object="bayesmodelClass"),
  definition=function(object)
  { cat("\n")
    if(is.null(object@convergence))
    { cat("convergence: \n NULL\n\n")
    } else cat("convergence: \n",object@convergence, "\n\n")
    if(is.null(object@results))
    { cat("results: \n NULL\n\n")
    } else {cat("results: \n");print(object@results);cat("\n")}
    if(is.null(object@jointpost))
    { cat("jointpost: \n NULL\n\n")
    } else {cat("jointpost: \n");print(head(object@jointpost));cat( "........\n\n")}
    if(is.null(object@nodes))
    { cat("nodes: \n NULL\n\n")
    } else cat("nodes: \n",object@nodes, "\n\n")
    if(is.null(object@model))
    { cat("model: \n NULL\n\n")
    } else cat("model: \n",object@model, "\n\n")
    if(is.null(object@chains))
    { cat("chains: \n NULL\n\n")
    } else cat("chans: \n",object@chains, "\n\n")
    if(is.null(object@burn))
    { cat("burn: \n NULL\n\n")
    } else cat("burn: \n",object@burn, "\n\n")
    if(is.null(object@update))
    { cat("update: \n NULL\n\n")
    } else cat("update: \n",object@update, "\n\n")})
    
    
    
################################################################################
################################################################################
#' @description This function calculates and plots the Gelman-Rubin convergence statistic, 
#' as modified by Brooks and Gelman (1998).
#'
#' @details This function is an alias of the function \code{\link{samplesBgr}} from the 
#' package \pkg{BRugs}. The original function was modified to create a diagnostic 
#' plots in other format as it was implemented in the original version. For more 
#' details see the function \code{\link{samplesBgr}} from the package \pkg{BRugs}.
#' \cr \cr
#' This function is not intended to be called directly but is internally called
#' by \code{\link{diagnostics}} within the functions                                          
#' \code{\link{rrisk.BayesZIP}} and \code{\link{rrisk.BayesPEM}}.
#'
#' @name rrisk.samplesBgr
#' @aliases rrisk.samplesBgr
#' @title Plot the Gelman-Rubin convergence statistic
#' @usage rrisk.samplesBgr(node, beg=samplesGetBeg(), end=samplesGetEnd(),
#'  firstChain=samplesGetFirstChain(), lastChain=samplesGetLastChain(),
#'  thin=samplesGetThin(), bins=50, plot=TRUE, ...)
#' @param node character vector of length 1, name of a variable in the model
#' @param beg argument to select a slice of monitored values corresponding to iterations \code{beg:end}
#' @param end argument to select a slice of monitored values corresponding to iterations \code{beg:end}
#' @param firstChain argument to select a sub group of chains to calculate the Gelman-Rubin convergence statistics for. Number of chains must be larger than one
#' @param lastChain argument to select a sub group of chains to calculate the Gelman-Rubin convergence statistics for. Number of chains must be larger than one
#' @param thin only use every thin-th value of the stored sample for statistics
#' @param bins number of blocks
#' @param plot logical, whether to plot the BGR statistics or only return the values. If \code{TRUE}, values are returned invisibly
#' @param ... further graphical parameters as in \code{par}
# @return Returns ...
# @note Some notes...
# @seealso Nothing...
#' @keywords manip
#' @export
#' @examples
#' \dontrun{rrisk.samplesBgr("se")}                                              )

rrisk.samplesBgr<-function(node, beg=samplesGetBeg(), end=samplesGetEnd(),
  firstChain=samplesGetFirstChain(), lastChain=samplesGetLastChain(),
  thin=samplesGetThin(), bins=50, plot=TRUE, ...)
{
  #if (plot && is.null(ask))
  if (plot)
  {
        if (is.R())
            ask <- !((dev.cur() > 1) && !dev.interactive())
        else ask <- !((dev.cur() > 1) && !interactive())
  }
  oldBeg <- samplesGetBeg()
  oldEnd <- samplesGetEnd()
  oldFirstChain <- samplesGetFirstChain()
  oldLastChain <- samplesGetLastChain()                  
  oldThin <- samplesGetThin()
  on.exit({
    samplesSetBeg(oldBeg)
    samplesSetEnd(oldEnd)
    samplesSetFirstChain(oldFirstChain)
    samplesSetLastChain(oldLastChain)
    samplesSetThin(oldThin)})
  beg <- max(beg, modelAdaptivePhase())
  samplesSetBeg(beg)
  samplesSetEnd(end)
  samplesSetFirstChain(firstChain)
  samplesSetLastChain(lastChain)
  thin <- max(c(thin, 1))
  samplesSetThin(thin)
  mons <- samplesMonitors(node)
  result <- lapply(mons, plotBgr, bins = bins, plot = plot,...)
  names(result) <- mons
  if (plot)
    invisible(result)
  else return(result)
} # end of function rrisk.samplesBgr()


################################################################################
################################################################################
#' @description This function provides a GUI for diagnostic plots to check convergence in Markov
#' chain Monte-Carlo (MCMC) models provided by 
#' \code{\link{rrisk.BayesZIP}} and \code{\link{rrisk.BayesPEM}}.
#'
#' @details The argument \code{nodes} denotes the node(s) to be used for diagnostic plots 
#' of the MCMC chains. The user is interactively requested to confirm whether 
#' the convergence has been reached. In this case the function returns 
#' \code{TRUE} otherwise \code{FALSE}.
#' \cr \cr
#' This function is not intended to be called directly but is internally called
#' within \code{\link{rrisk.BayesZIP}} or \code{\link{rrisk.BayesPEM}}.
#'
#' @name diagnostics
#' @aliases diagnostics
#' @title Diagnostic plots for MCMC models provided by rrisk.BayesPEM and rrisk.BayesZIP functions
#' @usage diagnostics(nodes, plots=FALSE)
#' @param nodes character string, the name of parameters(s)
#' @param plots logical, if \code{TRUE} the diagnostic plots will be displayed in separate windows
#' @return Returns \code{TRUE} if the user confirm convergence. Otherwise the 
#' function returns \code{FALSE}.
# @note Some notes...
#' @seealso For more details see documentation to the functions \code{\link{samplesBgr}}, 
#' \code{\link{plotHistory}} and \code{\link{plotDensity}} from the package \pkg{BRugs}.
#' @keywords manip
#' @export
# @references Literatur to be added...
#' @examples
#' \dontrun{diagnostics(nodes)}

diagnostics<-function(nodes, plots=FALSE)
{ 
  #-----------------------------------------------------------------------------
  # what happends by pressing "not converged" button
  #-----------------------------------------------------------------------------
  onNO<-function(...)
  {
     assign("output",value=FALSE,envir=tempEnvir)
     tkdestroy(diagPlotWindow)
  } # end of onNO()
  #-----------------------------------------------------------------------------
  # what happends by pressing "converged" button
  #-----------------------------------------------------------------------------
  onYES<-function(...)
  {
     assign("output",value=TRUE,envir=tempEnvir)
     tkdestroy(diagPlotWindow)
  } # end of onNO()
  #-----------------------------------------------------------------------------
  # what happends by pressing "next plot" button
  #-----------------------------------------------------------------------------
  onNext<-function(...)
  { plotNumber<-get("plotNumber",envir=tempEnvir)
    if(plotNumber==1){
      tkpack.forget(imgPlot1)
      assign("plotNumber",value=2,envir=tempEnvir)
      tkpack(imgPlot2,side="top")
    } else if (plotNumber==2){
      tkpack.forget(imgPlot2)
      assign("plotNumber",value=3,envir=tempEnvir)
      tkpack(imgPlot3,side="top")
    } else if (plotNumber==3){
      tkpack.forget(imgPlot3)
      assign("plotNumber",value=1,envir=tempEnvir)
      tkpack(imgPlot1,side="top")
    }
    tkraise(diagPlotWindow)
  }  # end of onNext()

  #-----------------------------------------------------------------------------
  # display diagnostic plots in separate windows
  #-----------------------------------------------------------------------------
  if(plots)
  {
    X11();plotDiag(nodes,plotnumber=1)
    X11();plotDiag(nodes,plotnumber=2)
    X11();plotDiag(nodes,plotnumber=3)
  }

  #-----------------------------------------------------------------------------
  # define help variables for output
  #-----------------------------------------------------------------------------
  assign("tempEnvir",value=new.env(),envir=.GlobalEnv)
  assign("output",value=FALSE,envir=tempEnvir)
  assign("plotNumber",value=1,envir=tempEnvir)

  diagPlotWindow<-tktoplevel(width=100,height=100)
  tkwm.title(diagPlotWindow,"Graphical convergence diagnostics for Bayesian models")
  tkwm.resizable(diagPlotWindow,0,0)  # fixed size, not resizeable
  tkwm.maxsize(diagPlotWindow,800,600)
  tkwm.minsize(diagPlotWindow,800,600)
  allFrame<-tkframe(diagPlotWindow)
  headingFont2<-tkfont.create(weight="bold",size=10)
  imageFrame<-tkframe(allFrame)

  imgPlot1<-tkrplot(imageFrame,fun=function()plotDiag(nodes,plotnumber=1),hscale=2,vscale=1.4);
  imgPlot2<-tkrplot(imageFrame,fun=function()plotDiag(nodes,plotnumber=2),hscale=2.0,vscale=1.4)
  imgPlot3<-tkrplot(imageFrame,fun=function()plotDiag(nodes,plotnumber=3),hscale=1.0,vscale=1.4)

  tkpack(imageFrame,fill="both",expand="yes")
  tkpack(imgPlot1)
  buttonsFrame<-tkframe(allFrame)
  nextButton<-ttkbutton(buttonsFrame,width=nchar("Show next plot")+2, text="Show next plot",command=onNext)
  convButton<-ttkbutton(buttonsFrame,width=nchar("Converged")+2, text="Converged",comman=onYES)
  notconvButton<-ttkbutton(buttonsFrame,width=nchar("Not converged")+2, text="Not converged",command=onNO)
  tkpack(nextButton,convButton,notconvButton,padx=c(20,20),side="left")
  tkpack(buttonsFrame,pady=c(0,0))
  tkpack(allFrame,fill="both",expand="yes",padx=c(15,15),pady=c(15,15))

  tkfocus(diagPlotWindow)
  tkwait.window(diagPlotWindow)
  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  return(get("output",envir=tempEnvir))
} #  end of function diagnostics()


################################################################################
################################################################################
#' @name plotDiag
#' @aliases plotDiag
#' @title Auxiliary function
#' @usage plotDiag(nodes, plotnumber)
#' @description Auxiliary function function for plotting graphical convergence
#'  diagnostics for Bayes models
#' @usage diagnostics(nodes, plots=FALSE)
#' @param nodes character string, the name of parameters(s)
#' @param plotnumber single numerical value (1 for plotting Gelman-Rubin convergence
#'  statistic (\code{samplesBgr()}), 2 for plotting history (\code{plotDensity}) and
#'  3 for plotting autocorrelation plots (\code{plotAutoC}))
#' @keywords manip
#' @export   

  plotDiag<-function(nodes,plotnumber=1)
  { len<-length(nodes)
    if(plotnumber==1){ 
      par(mfrow=c(2,len))
      for(i in 1:len){
        if(len<=3){
          maintext<-paste("Gelman-Rubin conv. \nstatistic for",nodes[i])
        } else {
          maintext<-paste("Gelman-Rubin conv. \nstatistic for",nodes[i])
        }
        rrisk.samplesBgr(nodes[i],main=maintext,cex.main=1)
        abline(h=1,lty="dashed")
      }
      for(i in 1:len){
        plotDensity(nodes[i],main=paste("Density plot of ",nodes[i]),cex.main=1)
      }
    } else if (plotnumber==2){ 
      if(len<=3){
       par(mfrow=c(len,1))
      } else  par(mfrow=c(2,3))
      #par(mfrow=c(len,1))
      for(i in 1:len){ 
        plotHistory(nodes[i],main=paste("Convergence monitored for node",nodes[i]),cex.main=1)
      }
    } else if(plotnumber==3){
      if(len<=3){
       par(mfrow=c(len,1))
      } else  par(mfrow=c(2,3))
      for(i in 1:len){ 
      plotAutoC(nodes[i],main=paste("Autocorrelation plot for node",nodes[i]),cex.main=1)
      }
    } # end if statement
  } # end of function plotDiag()


################################################################################
################################################################################
#' @description Bayesian PEM models provide the posterior distribution for the true prevalence (\code{pi}),
#' diagnostic sensitivity (\code{se}) and specificity (\code{sp}) for a given empirical prevalence 
#' estimate using physically pooled samples (if \code{k>1}) and priors for the model parameters. 
#' The misclassification parameters (\code{se} and \code{sp}) can be specified
#' at the level of the pool or individual level of testing. On the other side,
#' the function estimates the true prevalence based on the results
#' (\code{x/n}) of an application study with individual samples (if \code{k=1}) using a diagnostic test, for
#' which some prior information on sensitivity and specificity is available.
#' 
#' @details The Bayesian model for estimation prevalence, sensitivity and specificity has
#' in BRugs/Winbugs syntax following form for misclassification at the pool-level
#' (\code{k>1} and \code{misclass="pool"})
#' \preformatted{model{
#'
#'        pi ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'        se ~ dbeta(prior.se[1],prior.se[2])
#'
#'        sp ~ dbeta(prior.sp[1],prior.sp[2])
#'
#'        p.neg <- pow(1-pi,k)
#'
#'        p <- (1-p.neg)*se + p.neg*(1-sp)
#'
#'        x ~ dbin(p,n)
#'
#'      }}
#' for misclassifications at the individual level (\code{k>1} and \code{misclass="individual"})
#' \preformatted{model{
#'
#'        pi ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'        se ~ dbeta(prior.se[1],prior.se[2])
#'
#'        sp ~ dbeta(prior.sp[1],prior.sp[2])
#'
#'        ap <- pi*se + (1-pi)*(1-sp)
#'
#'        p <- 1- pow(1-ap,k)
#'
#'        x ~ dbin(p,n)
#'
#'      }}
#' and for comparison (\code{k>1})
#' \preformatted{model{
#'
#'        pi1 ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'        pi2 ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'        se ~ dbeta(prior.se[1],prior.se[2])
#'
#'        sp ~ dbeta(prior.sp[1],prior.sp[2])
#'
#'        x1 <- x 
#'
#'        x2 <- x  
#'  
#'        p.neg <- pow(1-pi1,k)
#'
#'        p.pos <- (1-p.neg)*se + p.neg*(1-sp)
#'
#'        x1 ~ dbin(p.pos,n)     
#'
#'        ap <- pi2*se + (1-pi2)*(1-sp)
#'
#'        p <- 1- pow(1-ap,k)
#'
#'        x2 ~ dbin(p,n) 
#'     
#'      }}
#' The application data (\code{k=1}) has one degree of freedom while the underlying model
#' has three unknown parameters. Thus, the model is not identifiable and informative
#' priors on at least two model parameters is required. The Bayesian model for estimation
#' prevalence, sensitivity and specificity takes a form
#' \preformatted{model{
#'
#'      x ~ dbin(p,n)
#'
#'      p <- pi * se + (1-pi) * (1-sp)
#'
#'      se ~ dbeta(prior.se[1],prior.se[2])
#'
#'      sp ~ dbeta(prior.sp[1],prior.sp[2])
#'
#'      pi  ~ dbeta(prior.pi[1],prior.pi[2])
#'
#'    }}
#'
#' @name rrisk.BayesPEM
#' @aliases rrisk.BayesPEM
#' @title Bayesian Prevalence estimation under misclassification (PEM)
#' @usage rrisk.BayesPEM(x, n, k, simulation=FALSE, prior.pi, prior.se, prior.sp,
#'  misclass="pool",chains=3, burn=1000, thin=1, update=10000,
#'  workdir=getwd(), plots=FALSE)
#' @param x scalar value for number of pools (\code{k>1}) or individual outcomes (\code{k=1}) with positive test result
#' @param n scalar value for number of pools tested (\code{k>1}) or the sample size in application study (\code{k=1})
#' @param k scalar value for number of individual samples physically combined into one pool;
#' set \code{k>1} for pooled sampling and \code{k=1} for individual sampling 
#' @param simulation logical, value \code{TRUE} means the function will be called within any simulation routine,
#' in this case the graphical diagnostic interface will not be invoked (default \code{FALSE})
#' @param prior.pi numeric vector containing parameters of a beta distribution as prior for prevalence \code{pi}, e.g. \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}
#' @param prior.se numeric vector containing parameters of a beta distribution as prior for sensitivity \code{se}, e.g. \code{se} \eqn{\sim} \code{prior.se(*,*)=beta(*,*)}
#' @param prior.sp numeric vector containing parameters of a beta distribution as prior for specificity \code{sp}, e.g. \code{sp} \eqn{\sim} \code{prior.sp(*,*)=beta(*,*)}
#' @param misclass character with legal character entries \code{pool}, \code{individual} or \code{compare}; ignored if k=1
#' @param chains positive single numeric value, number of independent MCMC chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period (default 1000)
#' @param thin positive single numeric value (default 1). The samples from every kth iteration will be used for 
#'        inference, where k is the value of thin. Setting \code{thin > 1} can help to reduce the autocorrelation
#'        in the sample.
#' @param update positive single numeric value, length of update iterations for estimation (default 10000)
#' @param workdir character string giving working directory to store temporary data (default \code{getwd()})
#' @param plots logical, if \code{TRUE} the diagnostic plots will be displayed in separate windows 
#' @return The function \code{rrisk.BayesPEM} returns an instance of the \code{\linkS4class{bayesmodelClass}}
#' class containing following informations
#' \item{\code{convergence}}{logical, whether the model has converged (assessed by the user)}
#' \item{\code{results}}{data frame containing statistics of the posterior distribution}
#' \item{\code{jointpost}}{data frame giving the joint posterior probability distribution}
#' \item{\code{nodes}}{names of the parameters jointly estimated by the Bayes model}
#' \item{\code{model}}{model in BRugs/Winbugs syntax as a character string}
#' \item{\code{chains}}{number of independent MCMC chains}
#' \item{\code{burn}}{length of burn-in period}
#' \item{\code{update}}{length of update iterations for estimation}
#' @note The convergence of the model is assessed by the user using diagnostic plots
#' provided by the \pkg{BRugs} package.
# @seealso nothing...
#' @keywords manip
#' @export
# @references Greiner, M., Belgoroski, N., NN (2011). Estimating prevalence using diagnostic data 
# with pooled samples: R function \code{rrisk.BayesPEM}. J.Stat.Software (in preparation).
#' @references Cowling, D.W., I.A. Gardner and W.O. Johnson (1999). Comparison of methods for estimation 
#'   of individual-level prevalence based on pooled samples, Prev.Vet.Med. 39: 211-225.
#' \cr
#' \cr
#' Rogan, W.J. and B. Gladen (1978). Estimating prevalence from the results of a screening test. Am. J. Epidemiol. 107: 71-76.
#' @examples
#' \donttest{
#' #------------------------------------------
#' # Example of PEM model (k>1)
#' #------------------------------------------
#' pi <- 0.01
#' se <- 0.96
#' se.n <- 1000
#' sp <- 0.99
#' sp.n <- 1000
#' n <- sample(10:1000,1,replace=TRUE)  # stochatsic sample size
#' k <- sample(5:50,1,replace=FALSE)    # stochastic pool size
#' 
#' # Parameters for beta priors 
#' se.a <- se.n*se+1
#' se.b <- se.n*(1-se)+1
#' sp.a <- sp.n*sp+1
#' sp.b <- sp.n*(1-sp)+1
#'            
#' # Random number of positive pools (x) considering uncertainty of se and sp
#' ap <- pi*se + (1-pi)*(1-sp)
#' p.pos <- 1-(1-ap)^k                             
#' x <- rbinom(1,prob=p.pos,size=n)
#' 
#' # Estimate using Bayes model at individual level
#' resPEM1 <- rrisk.BayesPEM(x=x, n=n,k=k, 
#'      prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
#'      misclass="individual")
#' resPEM1@@results
#' 
#' # Estimate using Bayes model at pool level
#' resPEM2 <- rrisk.BayesPEM(x=x, n=n,k=k, 
#'      prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
#'      misclass="pool")
#' resPEM2@@results
#' 
#' # Estimate using Bayes model compared
#' resPEM3 <- rrisk.BayesPEM(x=x, n=n,k=k, 
#'      prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
#'      misclass="compare")
#' resPEM3@@results
#'
#' #------------------------------------------
#' # Example of PEM model (k=1)
#' #------------------------------------------
#' # informative priors -> convergence is o.k.
#' resPEM4<-rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(12,22),
#'  prior.sp=c(22,55),prior.pi=c(1,1))
#' resPEM4@@results
#'
#' # non-informative priors -> convergence is not o.k.
#' resPEM5<-rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(1,1),
#'  prior.sp=c(1,1),prior.pi=c(1,1))
#' resPEM5@@results
#'
#' # informative priors -> convergence is o.k., without invoking
#' # graphical diagnostic interface
#' resPEM6<-rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(12,22),
#'  prior.sp=c(22,55),prior.pi=c(1,1))
#' resPEM6@@results
#' }

rrisk.BayesPEM <- function(x, n, k, simulation=FALSE,
    prior.pi, prior.se, prior.sp, misclass="pool",
    chains=3, burn=1000, thin=1, update=10000, workdir=getwd(), plots=FALSE)
{
  #-----------------------------------------------------------------------------
  # checking input
  #-----------------------------------------------------------------------------
  if(missing(x) | missing(n) | missing(k))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the following arguments: 'x', 'n', 'k'!",call.=FALSE)
  }
  if(missing(prior.pi) | missing(prior.se) | missing(prior.sp))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the following prior arguments: 'prior.pi', 'prior.se', 'prior.sp'!",call.=FALSE)
  }
  if(!is.numeric(chains) | !is.numeric(burn) | !is.numeric(update) | !is.numeric(x) | !is.numeric(n) | !is.numeric(k))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not numeric: 'chains', 'burn', 'update'!",call.=FALSE)
  }
  max.l <- max(length(x),length(n),length(k)) # variable pool size disabled
  if((length(x) == max.l)*(length(n) == max.l)*(length(k) == max.l) != 1)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, equal length of arguments 'x', 'n' and 'k' required.",call.=FALSE)
  }
  if(chains<=0 | burn<=0 | update<=0 | thin<1)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not positive: 'chains', 'burn', 'update', 'thin'!",call.=FALSE)
  }
  if(!is.element(misclass,c("pool","individual","compare")))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, misclass should be 'pool','individual' or 'compare'!",call.=FALSE)
  }
  try.setwd<-try(setwd(workdir),silent=TRUE)
  if(inherits(try.setwd, "try-error"))
  { on.exit(return(invisible(NA)))
    error.mess<-paste("INVALID INPUT, the working directory could not be found! Your input is \n",workdir)
    stop(error.mess,call.=FALSE)
  }
  if(!is.logical(plots))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'plots' should be of type logical!",call.=FALSE)
  }  
  #-----------------------------------------------------------------------------
  # plausibility check on priors
  #----------------------------------------------------------------------------- 
  if(min(prior.pi)<=0 | min(prior.se)<=0 | min(prior.sp)<=0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the beta prior distribution for prevalence, sensitivity and specificity should be strictly positive!",call.=FALSE)
  }
  if(length(prior.pi)!=2 | length(prior.se)!=2 | length(prior.sp)!=2)
  { on.exit(return(invisible(NA)))
    stop("Two parameters for the beta prior distribution for prevalence, sensitivity and specificity are required",call.=FALSE)
  }
      
  #-----------------------------------------------------------------------------
  # create a temporary directory
  #-----------------------------------------------------------------------------
  setwd(workdir)
  #if(!file.exists("tempBayes"))
  #{
  #  suppressWarnings(dir.create("tempBayes"))
  #} else if(file.exists("tempBayes"))
  #{
  #  # system("rmdir tempBayes")
  #  unlink("tempBayes", recursive=TRUE)
  #  suppressWarnings(dir.create("tempBayes"))
  #}
  # setwd(file.path(workdir,"tempBayes")) # no more `cd`

  #-----------------------------------------------------------------------------
  # replace the out list below with S4 return objects ...
  #-----------------------------------------------------------------------------
  out<-new("bayesmodelClass")
  #out <- list()

  #-----------------------------------------------------------------------------
  # create nodes (parameters to be estimated) for the model and table for parameters of prior distribution for each node
  #-----------------------------------------------------------------------------
  
  nodes <- c("pi","se","sp")
  nodes.prior <- matrix(1,nrow=3, ncol=2)
  colnames(nodes.prior) <- c("Parameter1","Parameter2")
  rownames(nodes.prior) <- nodes
  nodes.prior["pi",] <- prior.pi
  nodes.prior["se",] <- prior.se
  nodes.prior["sp",] <- prior.sp
  
  #-----------------------------------------------------------------------------
  # write model
  #-----------------------------------------------------------------------------
  if(k==1)
  { cat(
    "model {
      se ~ dbeta(prior.se[1],prior.se[2])
      sp ~ dbeta(prior.sp[1],prior.sp[2])
      pi  ~ dbeta(prior.pi[1],prior.pi[2])
      x ~ dbin(p,n)
      p <- pi*se+(1-pi)*(1-sp)
    }",file="model.txt")
  } else if (k>1)
  { if(misclass == "pool")
    {  cat(
      "model{
        pi ~ dbeta(prior.pi[1],prior.pi[2])
        se ~ dbeta(prior.se[1],prior.se[2])
        sp ~ dbeta(prior.sp[1],prior.sp[2])
        p.neg <- pow(1-pi,k)
        p <- (1-p.neg)*se + p.neg*(1-sp)
        x ~ dbin(p,n)
      }",file="model.txt")
    } else if(misclass == "individual")
    { cat(
      "model{
        pi ~ dbeta(prior.pi[1],prior.pi[2])
        se ~ dbeta(prior.se[1],prior.se[2])
        sp ~ dbeta(prior.sp[1],prior.sp[2])
        ap <- pi*se + (1-pi)*(1-sp)
        p <- 1- pow(1-ap,k)
        x ~ dbin(p,n)
      }",file="model.txt")
    } else if(misclass == "compare")
    { cat(
      "model{
        pi1 ~ dbeta(prior.pi[1],prior.pi[2])
        pi2 ~ dbeta(prior.pi[1],prior.pi[2])
        se ~ dbeta(prior.se[1],prior.se[2])
        sp ~ dbeta(prior.sp[1],prior.sp[2])
        x1 <- x 
        x2 <- x 
       
        p.neg <- pow(1-pi1,k)
        p.pos <- (1-p.neg)*se + p.neg*(1-sp)
        x1 ~ dbin(p.pos,n)
        
        ap <- pi2*se + (1-pi2)*(1-sp)
        p <- 1- pow(1-ap,k)
        x2 ~ dbin(p,n)
        
        d <- pi1 - pi2
        
      }",file="model.txt")
     nodes <- c("pi1","pi2","se","sp","d")
    }
  }

  #-----------------------------------------------------------------------------
  # write data
  #-----------------------------------------------------------------------------
  if(k>1)
  {data.set <- list(x=x,n=n,k=k,prior.se=prior.se,prior.sp=prior.sp,prior.pi=prior.pi)
  } else if(k==1)
  { data.set <- list(x=x,n=n,prior.se=prior.se,prior.sp=prior.sp,prior.pi=prior.pi)
  }

  dput(data.set, file = "data.txt",control="keepNA")

  #-----------------------------------------------------------------------------
  # check model, load data and compile three chains
  #-----------------------------------------------------------------------------
  cat("-------------------------------------------------------------------\n")
  cat("Begin model fitting...\n")
  modelCheck("model.txt")
  modelData("data.txt")
  modelCompile(numChains = chains)

  #-----------------------------------------------------------------------------
  # generate inits
  #-----------------------------------------------------------------------------
  modelGenInits()

  #-----------------------------------------------------------------------------
  # run model with burn-in
  #-----------------------------------------------------------------------------
  modelUpdate(burn,thin=thin)
  samplesSet(nodes)
  # run model with iterations to reach stationary joint distribution
  modelUpdate(update,thin=thin)
  
  #-----------------------------------------------------------------------------
  # überprüfe, ob model erfolgreich gefittet wurde
  #-----------------------------------------------------------------------------
  plotDiag.check<-function(nodes){
    X11()
    plot1<-plotDiag(nodes,plotnumber=1)
    X11()
    plot2<-plotDiag(nodes,plotnumber=2)
    plot3<-plotDiag(nodes,plotnumber=3)
    output<-list(plot1=plot1,plot2=plot2,plot3=plot3)
  }
  try.result<-try(plotDiag.check(nodes),silent=FALSE)
  if (inherits(try.result, "try-error")) {
    if(!is.null(dev.list())) {graphics.off()}
    on.exit(return("ERROR"))
    stop("Error occured during model fitting...",call.=FALSE)
  }

  #-----------------------------------------------------------------------------
  # check diagnostic plots for convergence
  #-----------------------------------------------------------------------------
  if(!simulation)
  { out@convergence <- FALSE
    convergence <- diagnostics(nodes,plots)
    if(!is.logical(convergence)){ # falls "ERROR" zurückgeliefert
      on.exit(return(convergence))
      stop("Any error occured during model fitting...ABORD function",call.=FALSE)
    } else { # falls TRUE oder FALSE zurückgeliefert
      if(!convergence)
      {
        on.exit(return(invisible(NA)))
        # tidy-up
        # setwd(workdir)
        #unlink("tempBayes",recursive=TRUE)
        cat("End model fitting...\n")
        cat("-------------------------------------------------------------------\n")
        stop("Process has been cancelled by the user due to non-convergence",call.=FALSE)
      }
      out@convergence <- convergence
    }
  }
  
  #-----------------------------------------------------------------------------
  # estimation results from BRugs package
  #-----------------------------------------------------------------------------

  results <- samplesStats(nodes,beg = samplesGetBeg(), end = samplesGetEnd(),
    firstChain = samplesGetFirstChain(),
    lastChain = samplesGetLastChain())
  out@results <- results


  #-----------------------------------------------------------------------------
  # collect posterior joint distribution
  #-----------------------------------------------------------------------------
  if(misclass != "compare") pi.out <- samplesSample("pi")
  if(misclass == "compare") {
    pi.out <- samplesSample("pi1")
    pi.out <- samplesSample("pi2")
    pi.out <- samplesSample("d")
    }
  se.out <- samplesSample("se")
  sp.out <- samplesSample("sp")
  if(misclass == "compare") {
    out@jointpost <- data.frame(pi=NA,se=NA,sp=NA)
    } else {
    out@jointpost <- data.frame(pi=pi.out,se=se.out,sp=sp.out)
    }

    
  #-----------------------------------------------------------------------------
  # full model description (which can be run in Winbugs if necessary)
  #-----------------------------------------------------------------------------

  out@model <- paste(" 
  # pi=prevalence on element level 
  # se=sensitivity of diagnostic method
  # sp=specificity of diagnostic method
  
  # model
  ",paste(suppressWarnings(read.table(file="model.txt",sep="\n",stringsAsFactors=FALSE))[,1],collapse="\n"),"
  # data
  ",paste(suppressWarnings(read.table(file="data.txt",sep="\n",stringsAsFactors=FALSE))[,1],collapse="\n"),"\n",collapse=" ")

  #-----------------------------------------------------------------------------
  # nodes names
  #-----------------------------------------------------------------------------
  out@chains <- chains
  out@nodes <- nodes
  out@burn<-burn
  out@update<-update
  cat("End model fitting...\n")
  cat("-------------------------------------------------------------------\n")

  #-----------------------------------------------------------------------------
  # tidy-up
  #-----------------------------------------------------------------------------
  file.remove("data.txt")
  file.remove("model.txt")
  # setwd(workdir)
  # unlink("tempBayes",recursive=TRUE)

  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  #write.table(out$model,file="doc.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  return(out)
} # end of function 




################################################################################
################################################################################
#' @description Zero-inflated Poisson data are count data with an excess number of zeros. The
#' ZIP model involves the Poisson parameter \code{lambda} and the prevalence
#' parameter \code{pi}.
#'
#' @details The ZIP model applies to count data and can be interpreted as a mixture
#' distribution with one component comprising the 'true' zeros and another component
#' of Poisson distributed values with density parameter \code{lambda}. The prevalence
#' parameter \code{pi} refers to the proportion of the second, true non-zero
#' component.
#' \cr
#' The Bayesian model for estimation prevalence and lambda parameter has
#' in BRugs/Winbugs syntax following form 
#' \preformatted{model{
#'
#'    lambda ~ dunif(prior.lambda[1],prior.lambda[2])
#'
#'    pi ~ dbeta(prior.pi[1],prior.pi[2]) 
#'
#'    for (i in 1:n) {  
#'                    
#'                 y[i]  ~ dpois(mu[i])
#' 
#'                 mu[i] <- I[i] * lambda  
#'
#'                I[i] ~ dbern(pi)  
#'        
#'    }   
#'                        
#'  }}
#'
#' @name rrisk.BayesZIP
#' @aliases rrisk.BayesZIP
#' @title Bayes estimation of a zero inflated Poisson (ZIP) model
#' @usage rrisk.BayesZIP(data, prior.lambda=c(1,10), prior.pi=c(0.8,1), simulation=FALSE,
#'  chains=3, burn=1000, thin=1, update=10000, workdir=getwd(), plots=FALSE)
#' @param data matrix or data frame, data set with positive integers, including zeros and of the minimal length 10.
#' @param prior.lambda numeric vector containing minimum and maximum of a uniform
#' distribution used as prior for the Poisson parameter \code{lambda}, e.g. \cr \code{lambda} \eqn{\sim} \code{prior.lambda(*,*)=unif(*,*)}
#' @param prior.pi numeric vector containing parameters of a beta distribution
#' describing prior knowledge about prevalence (proportion of contaminated samples), e.g. \cr \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}
#' @param simulation logical, value \code{TRUE} means the function will be called within any simulation routine,
#' in this case the graphical diagnostic interface will not be invoked (default \code{FALSE})
#' @param chains positive single numeric value, number of independent MCMC chains (default 3)
#' @param burn positive single numeric value, length of the burn-in period (default 1000)
#' @param thin positive single numeric value (default 1). The samples from every kth iteration will be used for 
#'        inference, where k is the value of thin. Setting \code{thin > 1} can help to reduce the autocorrelation
#'        in the sample.
#' @param update positive single numeric value, length of update iterations for estimation (default 10000)
#' @param workdir character string giving working directory to store temporary data (default \code{getwd()})
#' @param plots logical, if \code{TRUE} the diagnostic plots will be displayed in separate windows 
#' @return The function \code{rrisk.BayesZIP} returns an instance of the \code{\linkS4class{bayesmodelClass}}
#' class containing following informations
#' \item{\code{convergence}}{logical, whether the model has converged (assessed by the user)}
#' \item{\code{results}}{data frame containing statitsics of the posterior distribution}
#' \item{\code{jointpost}}{data frame giving the joint posterior probability distribution}
#' \item{\code{nodes}}{names of the parameters jointly estimated by the Bayes model}
#' \item{\code{model}}{model in BRugs/Winbugs syntax as a character string}
#' \item{\code{chains}}{number of independent MCMC chains}
#' \item{\code{burn}}{length of burn-in period}
#' \item{\code{update}}{length of update iterations for estimation}
#' @note The convergence of the model should be checked using the diagnostic plots
#' see the package \pkg{BRugs}, see also \pkg{zicounts}.
# @seealso nothing...
#' @keywords manip
#' @export
#' @references Bohning, D., E. Dietz, P. Schlattman, L. Mendonca, and U. Kirchner (1999). 
#' The zero-inflated Poisson model and the decayed, missing and filled teeth index in 
#' dental epidemiology. Journal of the Royal Statistical Society, Series A 162, 195-209.
#' @examples
#' \donttest{
#' #------------------------------------------
#' # Example of ZIP model
#' #------------------------------------------
#' # generate ZIP data
#' pi<-0.01
#' n<-200
#' lambda<-3.5
#' zip.data<-rep(0,n)
#' zip.data[sample(1:n,n*pi,replace=FALSE)]<-rpois(n*pi,lambda=lambda)
#'
#' # estimate using Bayes model for zero inflated data
#' resZIP<-rrisk.BayesZIP(data=zip.data, prior.lambda=c(0,100),prior.pi=c(1,1),
#'  burn=100,update=1000)
#' resZIP@@results
#'
#' # estimate using Bayes model for zero inflated data without invoking
#' # graphical diagnostic interface
#' rrisk.BayesZIP(data=zip.data, prior.lambda=c(0,100),prior.pi=c(1,1),
#'  burn=100,update=1000,simulation=TRUE)
#'
#' # compare with naive results ignoring ZIP model
#' pi.crude <- sum(zip.data>0)/n
#' lambda.crude <- mean(zip.data)
#' print(pi.crude)
#' print(lambda.crude)
#' resZIP@@results
#' }

rrisk.BayesZIP <- function(data, prior.lambda=c(1,10), prior.pi=c(0.8,1), simulation=FALSE,
  chains=3, burn=1000, thin=1, update=10000, workdir=getwd(), plots=FALSE)
{
  #-----------------------------------------------------------------------------
  # checking input
  #-----------------------------------------------------------------------------
  if(missing(data))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, missing one or more of the function arguments: 'data', 'prior.lambda', 'prior.pi'!",call.=FALSE)
  }
  if(!is.numeric(data)| !is.numeric(prior.lambda) | !is.numeric(prior.pi) | !is.numeric(chains) | !is.numeric(burn) | !is.numeric(update))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not numeric: 'data', 'prior.lambda', 'prior.pi', 'chains', 'burn', 'update'!",call.=FALSE)
  }
  if(!is.logical(plots))
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'plots' should be of type logical!",call.=FALSE)
  }
  if(chains<=0 | burn<=0 | update<=0 | thin<1)
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, one or more of the following arguments is not positive: 'chains', 'burn', 'update', 'thin'!",call.=FALSE)
  }
  if(length(prior.lambda)!=2 |length(prior.pi)!=2 )
  { on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the arguments 'prior.lambda' and 'prior.pi' should be of length 2!",call.=FALSE)
  }
  try.setwd<-try(setwd(workdir),silent=TRUE)
  if(inherits(try.setwd, "try-error"))
  { on.exit(return(invisible(NA)))
    error.mess<-paste("INVALID INPUT, the working directory could not be found! Your input is \n",workdir)
    stop(error.mess,call.=FALSE)
  }  
  #-----------------------------------------------------------------------------
  # plausibility check on data
  #-----------------------------------------------------------------------------
  if(length(data)<10)
  { on.exit(return(invisible(NA)))
    stop("Data set too small for this purpose!",call.=FALSE)
  }
  if(min(data)<0)
  { on.exit(return(invisible(NA)))
    stop("Negative counts in data set are not allowed!",call.=FALSE)
  }
  if(any(abs(round(data)-data)!=0))
  { on.exit(return(invisible(NA)))
    stop("Data set contains non-integer values!",call.=FALSE)
  }
  #-----------------------------------------------------------------------------
  # plausibility check on priors
  #----------------------------------------------------------------------------- 
  if(min(prior.pi)<=0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the beta prior distribution for prevalence should not be negative!",call.=FALSE)
  }
  if(min(prior.lambda)<0)
  { on.exit(return(invisible(NA)))
    stop("Parameters of the uniform prior distribution used as prior for the Poisson parameter should be strictly positive!",call.=FALSE)
  }

  #-----------------------------------------------------------------------------
  # create a temporary directory
  #-----------------------------------------------------------------------------
  setwd(workdir)
#  if(!file.exists("tempBayes"))
#  {
#    suppressWarnings(dir.create("tempBayes"))
#  } else if(file.exists("tempBayes"))
#  {
#    # system("rmdir tempBayes")
#    unlink("tempBayes", recursive=TRUE)
#    suppressWarnings(dir.create("tempBayes"))
#  }
  # setwd(file.path(workdir,"tempBayes")) # no more `cd`

  #-----------------------------------------------------------------------------
  # replace the out list below with S4 return objects ...
  #-----------------------------------------------------------------------------
  out<-new("bayesmodelClass")
  #out <- list()

  #-----------------------------------------------------------------------------
  # write model
  #-----------------------------------------------------------------------------
  cat(
  "model{
    lambda ~ dunif(prior.lambda[1],prior.lambda[2])
    pi ~ dbeta(prior.pi[1],prior.pi[2])
    for (i in 1:n) {
    y[i]  ~ dpois(mu[i])
    mu[i] <- I[i] * lambda
    I[i] ~ dbern(pi)
    }
  }",file="model.txt")

  #-----------------------------------------------------------------------------
  # write data
  #-----------------------------------------------------------------------------
  data.set <- list(y=data, n=length(data), prior.lambda=prior.lambda,prior.pi=prior.pi)
  dput(data.set, file = "data.txt",control="keepNA")

  #-----------------------------------------------------------------------------
  # check model, load data and compile three chains
  #-----------------------------------------------------------------------------
  cat("-------------------------------------------------------------------\n")
  cat("Begin model fitting...\n")
  modelCheck("model.txt")
  modelData("data.txt")
  modelCompile(numChains = chains)

  #-----------------------------------------------------------------------------
  # generate inits
  #-----------------------------------------------------------------------------
  modelGenInits()

  #-----------------------------------------------------------------------------
  # run model with burn-in
  #-----------------------------------------------------------------------------
  modelUpdate(burn,thin=thin)
  samplesSet(c("pi","lambda"))
  # run model with iterations to reach stationary joint distribution
  modelUpdate(update,thin=thin)

  #-----------------------------------------------------------------------------
  # überprüfe, ob model erfolgreich gefittet wurde
  #----------------------------------------------------------------------------- 
  plotDiag.check<-function(nodes){
    X11()
    plot1<-plotDiag(nodes,plotnumber=1)
    X11()
    plot2<-plotDiag(nodes,plotnumber=2)
    plot3<-plotDiag(nodes,plotnumber=3)
    output<-list(plot1=plot1,plot2=plot2,plot3=plot3)
  }
  try.result<-try(plotDiag.check(nodes=c("pi","lambda")),silent=FALSE)
  if (inherits(try.result, "try-error")) {
    if(!is.null(dev.list())) {graphics.off()}
    on.exit(return("ERROR"))
    stop("Error occured during model fitting...",call.=FALSE)
  }

  #-----------------------------------------------------------------------------
  # check diagnostic plots for convergence
  #-----------------------------------------------------------------------------
  if(!simulation)
  { out@convergence <- FALSE
    convergence <- diagnostics(nodes=c("pi","lambda"),plots)
    if(!is.logical(convergence)){ # falls "ERROR" zurückgeliefert
      on.exit(return(convergence))
      stop("Any error occured during model fitting...ABORD function",call.=FALSE)
    } else { # falls TRUE oder FALSE zurückgeliefert
      if(!convergence)
      {
        on.exit(return(invisible(NA)))
        # tidy-up
        # setwd(workdir)
        #unlink("tempBayes",recursive=TRUE)
        cat("End model fitting...\n")
        cat("-------------------------------------------------------------------\n")
        stop("Process has been cancelled by the user due to non-convergence",call.=FALSE)
      }
      out@convergence <- convergence
    }
  }

  #-----------------------------------------------------------------------------
  # estimation results from BRugs package
  #-----------------------------------------------------------------------------
  results <- samplesStats(c("pi","lambda"),beg = samplesGetBeg(), end = samplesGetEnd(),
    firstChain = samplesGetFirstChain(),
    lastChain = samplesGetLastChain())
  out@results <- results

  #-----------------------------------------------------------------------------
  # collect posterior joint distribution
  #-----------------------------------------------------------------------------
  pi.out <- samplesSample("pi")
  lambda.out <- samplesSample("lambda")
  out@jointpost <- data.frame(pi=pi.out,lambda=lambda.out)

  #-----------------------------------------------------------------------------
  # full model description (which can be run in Winbugs if necessary)
  #-----------------------------------------------------------------------------
  out@model <- paste(" # lambda=Poisson density parameter
  # pi=prevalence of contaminated sampes
  #
  # model
  ",paste(suppressWarnings(read.table(file="model.txt",sep="\n",stringsAsFactors=FALSE))[,1],collapse="\n"),"
  # data
  ",paste(suppressWarnings(read.table(file="data.txt",sep="\n",stringsAsFactors=FALSE))[,1],collapse="\n"),"\n",collapse=" ")

  #-----------------------------------------------------------------------------
  # nodes names
  #-----------------------------------------------------------------------------
  out@nodes <- c("lambda","pi")
  out@burn<-burn
  out@update<-update
  out@chains<-chains
  cat("End model fitting...\n")
  cat("-------------------------------------------------------------------\n")

  #-----------------------------------------------------------------------------
  # tidy-up
  #-----------------------------------------------------------------------------
  file.remove("data.txt")
  file.remove("model.txt")
  #  setwd(workdir)
  #  unlink("tempBayes",recursive=TRUE)

  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  #write.table(out$model,file="doc.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  return(out)
} # end of function rrisk.BayesZIP()

    
    
    

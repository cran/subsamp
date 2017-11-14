#' Subsample Winner Algorithm for Variable Selection in Linear Regression with a Large Number of Variables
#'
#' This subsample winner algorithm (SWA) for regression with a large-p data (X, Y) selects the important variables (or features) among the p features X in explaining the response Y.  The SWA first uses a base procedure, here a linear regression, on each of subsamples randomly drawn from the p variables, and then computes the scores of all features, i.e., the p variables, according to the performance of these features collected in each of the subsample analyses. It then obtains the 'semifinalist' of the features based on the resulting scores and determines the 'finalists', i.e., the important features, from the 'semifinalist'.
#'
#' @param x: The n by p design matrix where n is the sample size and p is the demension.
#' @param y: The vector with length n.
#' @param s: The subsample size which must be greater than or equal to 2.
#' @param m: The number of subsample repetition. The default value is 1000.
#' @param qnum: The number of semi-finalists. It should be at least as big as s.
#' @param wplot: Boolean input showing whether the multi-panel weights plot is drawn to guide a selection of s, the subsample size. The default value is F (False), without the plot. If it is plotted, a user should look at the points above the elbow points in all panels, i.e. upper arm sets, where the current value of s is in the red color. Find the smallest s0 such that the upper arm sets, corresponding to next few s >= s0,  are similar.
#'
#' @return Output a list of the finalists, i.e. the important features with their p-values, estimated coefficients and other corresponding statistics.
#'
#' @keywords internal
#'
#' @export
#' @import stats graphics
#'
#' @author Yiying Fan, Xingye Qiao, and Jiayang Sun
#'
#' @references http:sr2c.case.edu/swa-reg/
#'
#' @examples
#' n <- 80; p <- 100
#' set.seed(2017)
#' x <- matrix(rnorm(n*p),n);
#' coefs <- c(.1, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5)
#' y <- x[, 1:length(coefs)]%*%coefs + rnorm(n)
#' d <- data.frame(y,x)
#'
#' ## The number of true features in d is 10. The first 3 features,
#' ## X1-X3, are not distinguishable from noises, while the next two
#' ## X4 and X5 have a weak signal, X6 has a moderate signal,
#' ## and X7-X10 have a strong signal in comparison with the noises
#' ## in the data. A good feature selection procedure should capture
#' ## X7-X10, at least. The ideal s here should be 10 or slightly
#' ## bigger.
#'
#' ## Using the default values for m, wplot and a smaller value of s
#' ## than desired:
#' subsamp(x, y, s=8, qnum=10)  #1st run
#' ## This run captured X7-X10 and X5, with an adjusted
#' ## R^2=0.9015.
#'
#' ## Next try a bigger s and include a diagnostic plot:
#' #  subsamp(x, y, s=10, qnum=10, wplot=TRUE) #2nd run
#' ## It captured X5-X10, as expected from the truth.
#' ## It is also good enough by looking at either the
#' ## adjusted R^2, 0.955.  The diagnostic weights plot indicated that
#' ## s=10 is a good choice.
#'
#' ## However, if a conservative user decided to try for an even
#' ## bigger m and q:
#' #  subsamp(x, y, s=10, qnum=12, m=1200, wplot=TRUE) #3rd run
#' ## It now definitely suggests selecting s=10, but this run only
#' ## captured X5,X6,X8,X9,X10 without X7, and added a spurious X62,
#' ## with a **smaller** adjusted R2 = 0.8827. Thus a user should
#' ## now stop and conclude with using the outcome from the 2nd run.
#'
#' ## Regardless, if the user kept increasing both s& m, we would have
#' #  subsamp(x, y, s=12, qnum=12, m=1500, wplot=TRUE)  #run4
#' ## Its outcome is same as that from run2.
#'
#' ## Double Assurance Procedure: This is to further assure
#' ## the outcome, by  applying the base procedure to
#' ## the combined features from reasonable SWA runs. Combining run1
#' ## and run2, with and without run3, lead to the same important
#' ## features as those from the run2:
#' ##
#' #  g <- lm(formula = y ~ X8 + X10 + X9 + X6 + X5 + X95 + X61 +X20+X7+X17+X73+X82+X47,d)
#' #  summary(g)
#' #  step(g)
#'
#' ## We did not include the outcome from run4 into the double
#' ## assurance procedure as its outputs is same as that of run2.


subsamp <- function(x, y, s, m=1000, qnum, wplot=FALSE)
{

  sub.weight <- function(x, y, s, m=1000){
    # x is data matrix and y is the vector of response:
    p<-dim(x)[2]

    # intial values for RSS vector and t-value matrix
    rss.full<-NULL; tval.full<-NULL

    for (k in 1:m){
      tval.new<-rep(0,p)

      index.new<- sort(sample(1:p, s,replace = FALSE))
      subsamp.new <- x[,index.new]
      Data.new <- data.frame(y,subsamp.new)
      ssreg.new <-lm(y~., Data.new)
      tval.new[index.new]<-summary(ssreg.new)$coefficients[-1,3]
      rss.new<-summary(ssreg.new)$sigma

      rss.full<-c(rss.full, rss.new)
      tval.full<-cbind(tval.full, tval.new)
    }

    tval.sel<-tval.full[,order(rss.full)[1:s]]
    rss.sel<-rss.full[order(rss.full)[1:s]]

    weight1<-rep(0,p)
    for (j in 1:p){
      if (length(which(tval.sel[j,]!=0)))
        weight1[j]<-sum(abs(tval.sel[j,])/rss.sel)/length(which(tval.sel[j,]!=0))
    }

    ord.final<-order(weight1,decreasing = T)
    weight.final<- cbind(ord.final, weight1[ord.final])
    return(weight.final)
  }

  weight0=sub.weight(x, y, s, m)

  ### Plot feature weigts
  if(wplot==TRUE){
    p<- dim(x)[2]
    num<- floor(p/4)  # max number of features ploted on the x-axis

    ### generate a vector of s
    svec=function(s){
      if (s>=7) 		ss=c((s-4), (s-2), s, (s+2), (s+4))
      else if (s>=5&s<=6)  	ss=c(3, 4, s, (s+2), (s+4))
      else if (s==4)  	ss=c(3, s, (s+2), (s+4), (s+6))
      else if (s>=2&s<=3)  	ss=c(s, (s+1), (s+2), (s+4), (s+6))
      else print("invalid s")
      return(ss)
    }

    ### generate a sequence of s
    ss=svec(s); s.index=which(ss==s)

    ### compute and plot weights
    par(mfrow=c(2,3))

    for(i in 1:length(ss)){
      if (i==s.index) {
        weight1=weight0
        plot(weight1[,2][1:num], type="n", ylab="weights of features"); title(paste("s=",ss[i]), col.main = "red")
        text(1:num, weight1[,2][1:num], as.character(weight1[,1][1:num]))

      }
      else   {
        weight1=sub.weight(x, y, ss[i], m)
        plot(weight1[,2][1:num], type="n", ylab="weights of features"); title(paste("s=",ss[i]))
        text(1:num, weight1[,2][1:num], as.character(weight1[,1][1:num]))
      }
    }
  }

  ### use semi-finalists to constuct a data frame
  temp<-data.frame(x);	ord.final<-weight0[,1][1:qnum]
  x.final<-temp[,ord.final];	Data.final<-data.frame(y,x.final)

  ### find the finalists
  lm.final=lm(y~.,Data.final)
  reg.final=summary(step(lm.final,trace=F))
  return(reg.final)
}



mbrglm.control <-
function (mbr.epsilon = 1e-06, mbr.maxit = 500,mbr.trace = FALSE
             , ...) 
{
   if (!is.numeric(mbr.epsilon) || mbr.epsilon <= 0) 
      stop("value of 'epsilon' must be > 0")
   if (!is.numeric(mbr.maxit) || mbr.maxit <= 0) 
      stop("maximum number of iterations must be > 0")
   list(mbr.epsilon = mbr.epsilon, mbr.maxit = mbr.maxit,mbr.trace = mbr.trace)
}

#' Model list2 function
#'
#' @param yorder.dex
#' @param ydf.dex.0
#' @param xdf.dex
#'
#' @return
#' @export
#'
#' @examples
model_list <- function(yorder.dex,ydf.dex.0,xdf.dex){

  mod.dex      <- list()
  mod.dex[[1]] <- yorder.dex
  mod.dex[[2]] <- list()
  mod.dex[[3]] <- list()
  for(i in 1:length(yorder.dex)){ mod.dex[[3]][[i]] <- list() }

  list.mod <- NULL
  for(i in 1:length(yorder.dex)){

    if(yorder.dex[i]==0){  ydf.dex <- 0 }
    if(yorder.dex[i]!=0){  ydf.dex <- ydf.dex.0 }

    mod.dex[[2]][[i]] <- ydf.dex

    for(j in 1:length(ydf.dex)){

      mod.dex[[3]][[i]][[j]] <- xdf.dex

      for(k in 1:length(xdf.dex)){

        list.mod <- rbind(list.mod,c(i,j,k))

      }

    }

  }

  save(list.mod, file = paste("model_list_",dataset,".RData", sep="") )

  return(list.mod)
}

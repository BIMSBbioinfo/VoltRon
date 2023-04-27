slotApply <- function(x,FUN,...){
  cl <- class(x)
  result <- list()
  for(i in slotNames(cl)){
    result[[i]] <- FUN(slot(x,i),...)
  }
  result
}

slotToList <- function(x){
  returnlist <- list()
  namesslot <- slotNames(x)
  for(cur_slot in namesslot)
    returnlist[[cur_slot]] <- slot(x, name = cur_slot)
  returnlist
}

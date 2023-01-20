#' transform_magick_image
#' 
#' apply given transformations to a magick image
#'
#' @param image magick image
#' @param extension name extension for the shiny input parameter
#' @param input input 
#' @param session session
#'
#' @return magick image
#' @export
#'
transform_magick_image <- function(image, extension, input, session){
  
  # rotate 
  input_rotate <- isolate(as.numeric(input[[paste0("rotate_", extension)]]))
  image <- image_rotate(image, degrees = input_rotate)
  
  # flip flop
  input_flipflop <- isolate(input[[paste0("flipflop_", extension)]])
  if(input_flipflop == "Flip"){
    image <- image_flip(image)
  } else if(input_flipflop == "Flop"){
    image <- image_flop(image)
  }
  
  # ggplot output
  image <- image_ggplot(image)
  
}
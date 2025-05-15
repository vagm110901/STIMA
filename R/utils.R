#' selectCoord(image)
#'
#' Allows the user to select multiple points on the provided image by clicking on it.
#'
#' @param image 
#' @return list containing the coordinates of the selected points with two elements: x and y.
#' @import imager
#' @export
selectCoord <- function(image) {
  plot(image)
  coordinates <- locator(type = "p")
  return(coordinates)
}
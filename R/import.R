####
# Xenium ####
####

#' ImportXenium
#'
#' Importing a Xenium data
#'
#' @param dir.path path to Xenium
#' @param selected_assay selected assay
#' @param assay_name the assay name of the SR object
#' @param ... additional parameters passed to \code{CreateSpaceRover}
#'
#' @import magick
#'
#' @export
#'
ImportXenium <- function (dir.path, selected_assay = "Gene Expression", assay_name = "Xenium", ...)
{
  # raw counts
  datafile <- paste0(dir.path, "/cell_feature_matrix.h5")
  if(file.exists(datafile)){
    rawdata <- Read10X_h5(filename = datafile)
    rawdata <- as.matrix(rawdata[[selected_assay]])
  } else {
    stop("There are no files named 'filtered_feature_bc_matrix.h5' in the path")
  }

  # image
  image_file <- paste0(dir.path, "/morphology_lowres.tif")
  if(file.exists(image_file)){
    image <-  magick::image_read(image_file)
  } else {
    stop("There are no spatial image files in the path")
  }

  # cell boundaries
  bound_file <- paste0(dir.path, "/cell_boundaries.csv.gz")
  if(file.exists(bound_file)){
    Xenium_boundaries <- read.csv(bound_file)
    Xenium_box <- apply(Xenium_boundaries[,-1], 2, range)
  } else {
    stop("There are no files named 'cell_boundaries.csv.gz' in the path")
  }

  # coordinates
  coord_file <- paste0(dir.path, "/cells.csv.gz")
  if(file.exists(coord_file)){
    Xenium_coords <- read.csv(file = coord_file)
    coords <- as.matrix(Xenium_coords[,c("x_centroid", "y_centroid")])
    colnames(coords) <- c("x","y")
    coords[,2] <- max(coords[,2]) - coords[,2] + min(coords[,2])
    coords <- rescaleXeniumCells(coords, Xenium_box, image)
  } else {
    stop("There are no files named 'cells.csv.gz' in the path")
  }

  # create SpaceRover
  CreateSpaceRover(rawdata, metadata = NULL, image, coords, main.assay = assay_name, assay.type = "cell", ...)
}

#' rescaleXeniumCells
#'
#' rescale Xenium cells coordinates for image registration
#'
#' @param cells coordinates of the cells from the Xenium assays
#' @param bbox the surrounding box of the Xenium cell coordinates
#' @param image reference image
#'
rescaleXeniumCells <- function(cells, bbox, image){

  # get image scales
  scales <- unlist(image_info(image)[c("width","height")])

  # rescale cell locations
  cells[,1] <- (cells[,1] - bbox[1,1])/(bbox[2,1] - bbox[1,1])
  cells[,1] <- cells[,1] * scales[1]
  cells[,2] <- (cells[,2] - bbox[1,2])/(bbox[2,2] - bbox[1,2])
  cells[,2] <- cells[,2] * scales[2]

  # return
  return(cells)
}

####
# Visium ####
####

#' ImportVisium
#'
#' @param dir.path path to Xenium
#' @param assay_name the assay name
#' @param inTissue if TRUE, only barcodes that are in the tissue will be kept (default: TRUE)
#' @param ... additional parameters passed to \code{CreateSpaceRover}
#'
#' @import hdf5r
#' @import magick
#' @import jsonlite
#'
#' @export
#'
ImportVisium <- function(dir.path, assay_name = "Visium", InTissue = TRUE, ...)
{
  # raw counts
  listoffiles <- list.files(dir.path)
  datafile <- listoffiles[grepl("filtered_feature_bc_matrix.h5", listoffiles)][1]
  datafile <- paste0(dir.path, "/", datafile)
  if(file.exists(datafile)){
    rawdata <- Read10X_h5(filename = datafile)
    rawdata <- as.matrix(rawdata)
  } else {
    stop("There are no files named 'filtered_feature_bc_matrix.h5' in the path")
  }

  # image
  image_file <- paste0(dir.path, "/spatial/tissue_lowres_image.png")
  if(file.exists(image_file)){
    image <-  magick::image_read(image_file)
  } else {
    stop("There are no spatial image files in the path")
  }

  # coordinates
  coords_file <- paste0(dir.path, "/spatial/tissue_positions.csv")
  if(file.exists(coords_file)){
    coords <- read.csv(file = coords_file)
  } else {
    stop("There are no files named 'tissue_positions.csv' in the path")
  }
  coords$pxl_row_in_fullres <- max(coords$pxl_row_in_fullres) - coords$pxl_row_in_fullres + min(coords$pxl_row_in_fullres)
  if(InTissue){
    coords <- coords[coords$in_tissue==1,]
    rawdata <- rawdata[,colnames(rawdata) %in% coords$barcode]
  }
  spotID <- coords$barcode
  coords <- as.matrix(coords[,c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  colnames(coords) <- c("x", "y")
  rownames(coords) <- spotID

  # scale coordinates
  scale_file <- paste0(dir.path, "/spatial/scalefactors_json.json")
  if(file.exists(scale_file)){
    scalefactors <- jsonlite::read_json(path = scale_file)
    scales <- scalefactors$tissue_lowres_scalef
    coords <- coords*scales
  } else {
    stop("There are no files named 'scalefactors_json.json' in the path")
  }

  # create SpaceRover
  CreateSpaceRover(rawdata, metadata = NULL, image, coords, main.assay = assay_name, assay.type = "spot", ...)
}

####
# GeoMx ####
####

#' ImportGeoMx
#'
#' @param dir.path path to GeoMx run
#' @param assay_name the assay name
#' @param pkc_file path to the pkc file
#' @param ... additional parameters passed to \code{CreateSpaceRover}
#'
#' @import dplyr
#' @import jsonlite
#' @import GeomxTools
#' @import xlsx
#'
#' @export
#'
ImportGeoMx <- function(dir.path, pkc_file, summarySegment, summarySegmentSheetName, assay_name = "GeoMx", ...)
{
  # Get dcc file
  dcc_files <- dir(dir.path, pattern = ".dcc$", full.names = TRUE)
  dcc_files <- dcc_files[!grepl("A01.dcc$", dcc_files)]
  dcc_filenames <- dir(dir.path, pattern = ".dcc$", full.names = FALSE)
  dcc_filenames <- dcc_filenames[!grepl("A01.dcc$", dcc_filenames)]
  dcc_filenames <- gsub(".dcc$", "", dcc_filenames)
  dcc_filenames <- gsub("-", "_", dcc_filenames)
  dccData <- sapply(dcc_files, GeomxTools::readDccFile, simplify = FALSE, USE.NAMES = FALSE)
  names(dccData) <- dcc_filenames

  # merge dcc files
  rawdata <- NULL
  for(i in 1:length(dccData)){
    cur_data <- dccData[[i]]$Code_Summary
    colnames(cur_data) <- c("RTS_ID", dcc_filenames[i])
    if(i == 1){
      rawdata <- cur_data
    } else {
      suppressMessages(rawdata <- rawdata %>% full_join(cur_data))
    }
  }
  rawdata[is.na(rawdata)] <- 0
  rownames(rawdata) <- rawdata$RTS_ID
  rawdata <- rawdata[,!colnames(rawdata) %in% "RTS_ID"]

  # get pkc file
  pkcdata <- readPKCFile(pkc_file)

  # get genes
  NegProbes <- pkcdata$RTS_ID[pkcdata$Target == "NegProbe-WTX"]
  rawdata <- rawdata[!rownames(rawdata) %in% NegProbes, ]
  rownames(rawdata) <- pkcdata$Target[match(rownames(rawdata), pkcdata$RTS_ID)]
  rawdata <- as.matrix(rawdata)

  # get segment summary
  segmentsummary <- xlsx::read.xlsx(summarySegment, sheetName = summarySegmentSheetName)

  # get image
  image <- image_read(paste0(dir.path, "/geomx_lowres.tif"))
  geomx_image_info <- image_info(image)

  # get coordinates
  coords <- segmentsummary[,c("X","Y")]
  colnames(coords) <- c("x", "y")
  rownames(coords) <- segmentsummary$ROI.name

  # convert coordinates
  xRatio = geomx_image_info$width/segmentsummary$Scan.Width[1]
  yRatio = geomx_image_info$height/segmentsummary$Scan.Height[1]
  coords$x = (coords$x - segmentsummary$Scan.Offset.X[1]) * xRatio
  coords$y = (coords$y - segmentsummary$Scan.Offset.Y[1]) * yRatio
  coords$y = geomx_image_info$height - coords$y
  coords <- as.matrix(coords)

  # create SpaceRover
  CreateSpaceRover(rawdata, metadata = NULL, image, coords, main.assay = assay_name, assay.type = "ROI", ...)
}

#' DemuxGeoMx
#'
#' @param object a SpaceRover object
#'
#' @export
#'
DemuxGeoMx <- function(object, scale_param = 800)
{
  # get images
  images <- Image(object)

  # check if there are only one assay in the object
  if(length(images) > 1)
    stop("You can only subset a spaceRover assay with one image")

  # scale
  imageinfo <- image_info(images[[1]])
  scale_factor <- imageinfo$width/800
  scale_param <- paste0(scale_param, "x")
  images <- image_scale(images[[1]], scale_param)

  # get the ui and server
  if (interactive()){
    # ui <- tagList(
    ui <- fluidPage(
      # use javascript extensions for Shiny
      useShinyjs(),

      sidebarLayout(position = "left",

                    # Side bar
                    sidebarPanel(

                      h4("Selected Sections"),
                      fluidRow(
                        htmlOutput("summary"),
                        br(),
                        column(12,shiny::actionButton("resetpoints", "Reset Points")),
                        br(),
                        column(12,shiny::actionButton("addbox", "Add Box")),
                        br(),
                        column(12,shiny::actionButton("done", "Done")),
                      ),
                      br(),
                      # panel options
                      width = 3,
                    ),

                    mainPanel(

                      # main image
                      br(),
                      br(),
                      fluidRow(
                        imageOutput("cropped_image", click = "choose_corner", height = "1000px")
                      ),

                      # panel options
                      width = 9
                    )
      )
    )

    server <- function(input, output, session) {

      # selected corners
      selected_corners <- reactiveVal(tibble(x = numeric(), y = numeric()))

      # selected corner list
      selected_corners_list <- reactiveVal(tibble(box = character()))

      ### Main observable ####
      observe({

        # update summary
        output[["summary"]] <- renderUI({
          if(nrow(selected_corners_list()) > 0){
            corners <- selected_corners_list()$box
            print_selected <- paste0("Subset ", 1:length(corners), ": ", corners)
            HTML(paste(print_selected, collapse = '<br/>'))
          }
        })


        # output image
        output[["cropped_image"]] <- renderPlot({
          if(nrow(selected_corners()) < 2){
            corners <- apply(as.matrix(selected_corners()),2,as.numeric)
            # corners <- as.data.frame(corners)
            # colnames(corners) <- c("x", "y")
            # print(corners)
            # image_ggplot(images) +
            #   geom_point(mapping = aes(x = x, y = y), corners, size = 5, shape = 21, fill = "white")
            image_ggplot(images)
          } else {
            corners <- apply(as.matrix(selected_corners()),2,as.numeric)
            image_ggplot(images) +
              geom_rect(aes(xmin = corners[1,1], xmax = corners[2,1], ymin = corners[1,2], ymax = corners[2,2]),
                        fill = "green", alpha = 0.3, color = "black")
          }
        })
      })

      ## reset points ####
      observeEvent(input$resetpoints, {
        selected_corners() %>%
          filter(FALSE) %>% selected_corners()
      })

      ## add box ####
      observeEvent(input$addbox, {
        next_ind <- length(selected_corners_list()) + 1
        corners <- selected_corners()
        corners <- corners*scale_factor
        corners <- apply(corners,2,ceiling)
        print(corners)
        corners <- paste0(abs(corners[2,1]-corners[1,1]), "x",
                       abs(corners[2,2]-corners[1,2]), "+",
                       min(corners[,1]), "+", imageinfo$height - max(corners[,2]))
        selected_corners_list() %>%
          add_row(box = corners) %>%
          selected_corners_list()
        selected_corners() %>%
          filter(FALSE) %>% selected_corners()
      })

      ## select points on the image ####
      observeEvent(input$choose_corner, {
        if(nrow(selected_corners()) < 2) {
          keypoint <- data.frame(x = input[["choose_corner"]]$x,
                                 y = input[["choose_corner"]]$y)
          selected_corners() %>%
            add_row(x = keypoint$x, y = keypoint$y) %>%
            selected_corners()
        } else {
          selected_corners() %>%
            filter(FALSE) %>% selected_corners()
        }
      })

      ## select points on the image ####
      observeEvent(input$done, {

        RegisteredSpatialDataList <- list()
        corners_list <- selected_corners_list()
        print(corners_list)
        for(i in 1:length(corners_list$box)){
          print(corners_list$box[i])
          RegisteredSpatialDataList[[i]] <- subset(object, image = corners_list$box[i])
        }
        stopApp(list(RegisteredSpatialDataList = RegisteredSpatialDataList))
      })

    }

    shiny::runApp(shiny::shinyApp(ui, server))
  }
}

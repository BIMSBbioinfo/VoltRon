####
# Main Shiny App ####
####

#' registerSpatialData
#'
#' A mini shiny app to for registering images and spatial coordinates of multiple consequtive spatial datasets
#'
#' @param object_list a list of VoltRon (or Seurat) objects
#' @param reference_spatdata a reference spatial data set, used only if \code{object_list} is \code{NULL}
#' @param query_spatdata a query spatial data set, used only if \code{object_list} is \code{NULL}
#' @param keypoints a list of tables, each points to matching keypoints from registered images.
#' @param shiny.options a list of shiny options (launch.browser, host, port etc.) passed \code{options} arguement of \link{shinyApp}. For more information, see \link{runApp}
#'
#' @import shiny
#' @importFrom shinyjs useShinyjs show hide
#' @importFrom stats median
#' @importFrom magick image_read
#'
#' @export
registerSpatialData <- function(object_list = NULL, reference_spatdata = NULL, query_spatdata = NULL, keypoints = NULL, 
                                shiny.options = list(launch.browser = getOption("shiny.launch.browser", interactive()))) {

  ## Importing images ####

  # check object classes
  # if the input is not a list, switch to reference vs query mode
  if(!is.null(object_list)){
    # reference and query indices
    spatdata_list <- object_list
    centre <- floor(stats::median(1:length(spatdata_list)))
    register_ind <- setdiff(1:length(spatdata_list), centre)

    # reference vs query mode
  } else {

    spatdata_list <- list(reference_spatdata, query_spatdata)

    # reference and query indices
    centre <- 1
    register_ind <- 2
  }

  # get images from the list of objects
  orig_image_query_list_full <- lapply(spatdata_list, function(spat){
    assayname <- vrAssayNames(spat)
    channel_names <- vrImageChannelNames(spat[[assayname]])
    sapply(channel_names, function(chan){
      # vrImages(spat, assay = assayname, channel = chan)
      img <- vrImages(spat[[assayname]], channel = chan, as.raster = TRUE)
      if(!inherits(img, "Image_Array")){
        img <- magick::image_read(img)
      }
      img
    }, USE.NAMES = TRUE)
  })
  orig_image_query_list <- lapply(orig_image_query_list_full, function(spat_img) {
    return(spat_img[[1]])
  })
  orig_image_channelname_list <- lapply(spatdata_list, function(spat){
    assayname <- vrAssayNames(spat)
    vrImageChannelNames(spat[[assayname]])
  })
  
  ## UI and Server ####
  ui <- fluidPage(
    # use javascript extensions for Shiny
    waiter::useWaiter(),
    shinyjs::useShinyjs(),
    
    sidebarLayout(position = "left",
                  
                  # Side bar
                  sidebarPanel(
                    tags$style(make_css(list('.well', 'margin', '7%'))),
                    
                    # # specific settings for dealing with simultaneous click and brush events
                    # # https://jokergoo.github.io/2021/02/20/differentiate-brush-and-click-event-in-shiny/
                    tags$script(HTML("
                        $('#plot').mousedown(function(e) {
                            var parentOffset = $(this).offset();
                            var relX = e.pageX - parentOffset.left;
                            var relY = e.pageY - parentOffset.top;
                            Shiny.setInputValue('x1', relX);
                            Shiny.setInputValue('y1', relY);
                        }).mouseup(function(e) {
                            var parentOffset = $(this).offset();
                            var relX = e.pageX - parentOffset.left;
                            var relY = e.pageY - parentOffset.top;
                            Shiny.setInputValue('x2', relX);
                            Shiny.setInputValue('y2', relY);
                            Shiny.setInputValue('action', Math.random());
                        });
                    ")),
                    
                    # side bar for configuration
                    getSideBar(),
                    
                    # panel options
                    width = 3,
                  ),
                  
                  mainPanel(
                    
                    # Interface for the reference image
                    br(),
                    column(6,
                           
                           # Reference Images
                           getImageTabPanels(length(orig_image_query_list), orig_image_channelname_list, type = "ref"),
                           
                           br(),
                           
                           # Matching Alignment
                           getAlignmentTabPanel(length(orig_image_query_list), centre, register_ind),
                    ),
                    
                    # Interface for the query images
                    column(6,
                           
                           # Query Images
                           getImageTabPanels(length(orig_image_query_list), orig_image_channelname_list, type = "query"),
                           
                           br(),
                           
                           # Registered Query Images
                           getRegisteredImageTabPanels(length(orig_image_query_list), centre, register_ind)
                    ),
                    
                    # panel options
                    width = 9
                  )
    )
  )
  
  server <- function(input, output, session) {
    
    ## Manage interface ####
    updateParameterPanels(length(orig_image_query_list), input, output, session)
    updateSequentialTabPanels(input, output, session, centre, register_ind)
    
    ## Transform images ####
    trans_image_query_list <- transformImageQueryList(orig_image_query_list, input)
    
    ## get image and zoom info ####
    orig_image_query_info_list <- getImageInfoList(orig_image_query_list)
    zoom_list <- initiateZoomOptions(orig_image_query_info_list)
    manageImageZoomOptions(centre, register_ind, zoom_list, orig_image_query_list, orig_image_query_info_list, input, output, session)
    
    ## Manage reference and query keypoints ####
    xyTable_list <- initateKeypoints(length(orig_image_query_list), keypoints)
    manageKeypoints(centre, register_ind, xyTable_list, orig_image_query_list, orig_image_query_info_list, zoom_list, input, output, session)
    
    ## Image registration ####
    registration_mapping_list <- initiateMappings(length(spatdata_list))
    getManualRegisteration(registration_mapping_list, spatdata_list, orig_image_query_list, xyTable_list,
                           centre, register_ind, input, output, session)
    getAutomatedRegisteration(registration_mapping_list, spatdata_list, orig_image_query_list_full, orig_image_channelname_list,
                              centre, register_ind, input, output, session)
    
    ## Main observable ####
    observe({
      
      # output the list of query images
      getImageOutput(orig_image_query_list_full, orig_image_query_info_list, xyTable_list, zoom_list, centre, input, output, session)
      
    })
    
    ## Return values for the shiny app ####
    observeEvent(input$done, {
      
      # keypoints
      keypoints <- reactiveValuesToList(xyTable_list)
      
      # get keypoints and registered spatial datasets
      stopApp(
        list(keypoints = keypoints,
             registered_spat = getRegisteredObject(spatdata_list,
                                                   registration_mapping_list,
                                                   register_ind,
                                                   centre,
                                                   input,
                                                   reg_mode = ifelse(input$automatictag, "auto", "manual"),
                                                   image_list = orig_image_query_list))
      )
    })
  }
  
  # configure options
  shiny.options <- configure_shiny_options(shiny.options)
  
  # run app
  shiny::runApp(
    shiny::shinyApp(ui, server, options = list(host = shiny.options[["host"]], port = shiny.options[["port"]], launch.browser = shiny.options[["launch.browser"]]),
                    onStart = function() {
                      cat("Doing application setup\n")
                      onStop(function() {
                        cat("Doing application cleanup\n")
                      })
                    })
  )
}

####
# User Interface ####
####

#' getSideBar
#'
#' The UI for the app side bar
#'
#' @import shiny
#' 
#' @noRd
getSideBar <- function(len_images, channel_names, type){
  list(
    h4("Spatial Data Alignment"),
    fluidRow(
      column(12,shiny::checkboxInput("automatictag", "Automated", value = FALSE)),
      br(),
      column(12,selectInput("Method", "Method", 
                            choices = c("Homography", "Non-Rigid", "Homography + Non-Rigid"), 
                            selected = "Homography")),
      br(),
      column(12,selectInput("Matcher", "Matcher", 
                            choices = c("FLANN", "BRUTE-FORCE"), 
                            selected = "FLANN")),
      br(),
      column(12,textInput("GOOD_MATCH_PERCENT", "Match %", value = "0.20", width = "80%", placeholder = NULL)),
      column(12,textInput("MAX_FEATURES", "# of Features", value = "1000", width = "80%", placeholder = NULL)),
      br(),
      column(12,shiny::actionButton("register", "Register!")),
      br(),
    ),
    br(),
    fluidRow(
      column(12,shiny::htmlOutput("summary"))
    ),
    br(),
    fluidRow(
      column(12,shiny::actionButton("done", "Done")),
      br()
    ),
    br(),
    h4("How to use"),
    p(style="font-size: 12px;", strong("Single-L-click:"), "Select point"),
    p(style="font-size: 12px;", strong("Single-L-hold-drag:"), "Select area"),
    p(style="font-size: 12px;", strong("Double-L-click (selected area):"), "Zoom in"),
    p(style="font-size: 12px;", strong("Double-L-click (no area):"), "Zoom out")
  )
}

#' getImageTabPanels
#'
#' The UI for a set of reference/query spatial slides
#'
#' @param len_images the number of query images
#' @param channel_names the list of channel names for each image
#' @param type Either reference (ref) or query (query) image
#'
#' @noRd
getImageTabPanels <- function(len_images, channel_names, type){

  # get panel label
  label <- ifelse(type == "ref", "Ref. ", "Query ")

  # call panels
  do.call(tabsetPanel, c(id=paste0('image_tab_panel_',type), lapply(1:len_images, function(i) {
    tabPanel(paste0(label,i),
             br(),
             fluidRow(
               column(4, selectInput(paste0("rotate_", type, "_image",i), "Rotate (ClockWise):", choices = c(0,90,180,270), selected = 0)),
               column(4, selectInput(paste0("flipflop_", type, "_image",i), "Transform:", choices = c("None", "Flip", "Flop"), selected = "None")),
               column(4, selectInput(paste0("negate_", type, "_image",i), "Negate Image:", choices = c("No", "Yes"), selected = "No"))
             ),
             fluidRow(
               column(4, selectInput(paste0("channel_", type, "_image",i), "Channel:", choices = channel_names[[i]])),
               column(4, sliderInput(paste0("scale_", type, "_image",i), "Scale Parameter", min = 0, max = 1,  value = 1)),
               textOutput(paste0("scaleinfo_", type, "_image",i))
             ),
             fluidRow(imageOutput(paste0("plot_", type, i), 
                                  click = paste0("click_plot_", type, i),
                                  dblclick = paste0("dblclick_plot_", type, i),
                                  brush = brushOpts(paste0("brush_plot_", type, i), fill = "green", resetOnNew = TRUE))),
             br(),
             fluidRow(
               shiny::actionButton(paste0("remove_", type, i), "Remove Point")
             ),
    )
  })))
}

#' getRegisteredImageTabPanels
#'
#' The UI for a set of query spatial slides
#'
#' @param len_images the number of query images
#' @param centre center image index
#' @param register_ind query image indices
#'
#' @noRd
getAlignmentTabPanel <- function(len_images, centre, register_ind){

  # tab panels
  do.call(tabsetPanel, c(id='image_tab_panel_alignment',lapply(register_ind, function(i) {
    tabPanel(paste0("Ali. ",i, "->", centre),
             br(),
             fluidRow(imageOutput(paste0("plot_alignment",i)))
    )
  })))
}

#' getRegisteredImageTabPanels
#'
#' The UI for a set of query spatial slides
#'
#' @param len_images the number of query images
#' @param centre center image index
#' @param register_ind query image indices
#'
#' @return tabsetpanel
#'
#' @noRd
getRegisteredImageTabPanels <- function(len_images, centre, register_ind){

  # tab panels
  do.call(tabsetPanel, c(id='image_tab_panel_reg_query',lapply(register_ind, function(i) {
    tabPanel(paste0("Reg. ",i, "->", centre),
             br(),
             # column(6, sliderInput(paste0("plot_query_reg_alpha",i), label = "Alpha Level", min = 0, max = 1, value = 0.2)),
             fluidRow(
               column(12, align="center",
                      imageOutput(paste0("plot_query_reg",i))
               )
             )
    )
  })))
}

#' updateSequentialTabPanels
#'
#' A function for automatized selection of reference/query tab panels
#'
#' @param input input
#' @param output output
#' @param session session
#' @param centre center image index
#' @param register_ind query image indices
#'
#' @noRd
updateSequentialTabPanels <- function(input, output, session, centre, register_ind){

  # number of panels
  npanels <- length(register_ind) + 1

  # observe changes in the reference tab panel
  observeEvent(input$image_tab_panel_ref,{
    selected_panel <- input$image_tab_panel_ref
    selected_panel_ind <- as.numeric(strsplit(selected_panel, split = " ")[[1]][2])

    query_panel_ind <- (selected_panel_ind + 1)
    if(query_panel_ind == 1) query_panel_ind <- npanels
    updateTabsetPanel(session, "image_tab_panel_query", paste0("Query ", query_panel_ind))
    updateTabsetPanel(session, "image_tab_panel_reg_query", paste0("Reg. ",selected_panel_ind, "->", centre))

    if(selected_panel_ind == npanels)
      updateTabsetPanel(session, "image_tab_panel_ref", paste0("Ref. ", selected_panel_ind-1))
  })

  # observe changes in the query tab panel
  observeEvent(input$image_tab_panel_query,{
    selected_panel <- input$image_tab_panel_query
    selected_panel_ind <- as.numeric(strsplit(selected_panel, split = " ")[[1]][2])

    query_panel_ind <- (selected_panel_ind - 1)
    if(query_panel_ind == 0) query_panel_ind <- 1
    updateTabsetPanel(session, "image_tab_panel_ref", paste0("Ref. ", query_panel_ind))

    if(selected_panel_ind == 1){
      updateTabsetPanel(session, "image_tab_panel_query", paste0("Query ", selected_panel_ind+1))
      updateTabsetPanel(session, "image_tab_panel_reg_query", paste0("Reg. ",selected_panel_ind+1, "->", centre))
    } else {
      query_panel_ind <- selected_panel_ind
      updateTabsetPanel(session, "image_tab_panel_reg_query", paste0("Reg. ",query_panel_ind, "->", centre))
    }
  })

  # observe changes in the registered query tab panel
  observeEvent(input$image_tab_panel_reg_query,{
    selected_panel <- input$image_tab_panel_reg_query
    selected_panel_ind <- strsplit(selected_panel, split = " ")[[1]][2]
    selected_panel_ind <- as.numeric(strsplit(selected_panel_ind, split = "->")[[1]][1])
    updateTabsetPanel(session, "image_tab_panel_query", paste0("Query ", selected_panel_ind))
    selected_panel_ali <- gsub("Reg.", "Ali.", selected_panel)
    updateTabsetPanel(session, "image_tab_panel_alignment", selected_panel_ali)
  })

  # observe changes in the registered query tab panel
  observeEvent(input$image_tab_panel_alignment,{
    selected_panel <- input$image_tab_panel_alignment
    selected_panel_reg <- gsub("Ali.", "Reg.", selected_panel)
    updateTabsetPanel(session, "image_tab_panel_reg_query", selected_panel_reg)
  })
}

#' updateParameterPanels
#'
#' A function for managing which parameter panels or input boxes to appear on UI
#'
#' @param len_images the length of images
#' @param input input
#' @param output output
#' @param session session
#'
#' @noRd
updateParameterPanels <- function(len_images, input, output, session){

  # done event
  shinyjs::hide(id = "done")
  observeEvent(input$register, {
    shinyjs::show(id = "done")
  })

  # registration panels/buttons
  shinyjs::hide(id = "GOOD_MATCH_PERCENT")
  shinyjs::hide(id = "MAX_FEATURES")

  # hide scale parameters
  for(i in 1:len_images){
    shinyjs::hide(id = paste0("scale_ref_image",i))
    shinyjs::hide(id = paste0("scale_query_image",i))
    shinyjs::hide(id = paste0("scaleinfo_ref_image",i))
    shinyjs::hide(id = paste0("scaleinfo_query_image",i))
  }

  observeEvent(input$automatictag, {
    if(input$automatictag){
      
      # Method and Matcher
      updateSelectInput(session, "Method", choices = c("Homography", "Homography + Non-Rigid"), selected = "Homography")
      shinyjs::show(id = "Matcher")

      # show automatic registration parameters of BRUTE-FORCE
      if(input$Matcher == "BRUTE-FORCE"){
        shinyjs::show(id = "GOOD_MATCH_PERCENT")
        shinyjs::show(id = "MAX_FEATURES")
      } 
      if(input$Matcher == "FLANN"){
        shinyjs::hide(id = "GOOD_MATCH_PERCENT")
        shinyjs::hide(id = "MAX_FEATURES")
      }

      # show scale parameters
      for(i in 1:len_images){
        shinyjs::show(id = paste0("scale_ref_image",i))
        shinyjs::show(id = paste0("scale_query_image",i))
        shinyjs::show(id = paste0("scaleinfo_ref_image",i))
        shinyjs::show(id = paste0("scaleinfo_query_image",i))
      }

    } else {
      
      # Method and Matcher
      updateSelectInput(session, "Method", choices = c("Non-Rigid", "Homography + Non-Rigid"), selected = "Non-Rigid")
      shinyjs::hide(id = "Matcher")

      # hide automatic registration parameters of BRUTE-FORCE
      if(input$Matcher == "FLANN"){
        shinyjs::hide(id = "GOOD_MATCH_PERCENT")
        shinyjs::hide(id = "MAX_FEATURES")
      }

      # hide scale parameters
      for(i in 1:len_images){
        shinyjs::hide(id = paste0("scale_ref_image",i))
        shinyjs::hide(id = paste0("scale_query_image",i))
        shinyjs::hide(id = paste0("scaleinfo_ref_image",i))
        shinyjs::hide(id = paste0("scaleinfo_query_image",i))
      }
    }
  })

  observeEvent(input$Method, {
    if(grepl("FLANN", input$Matcher)){
      shinyjs::hide(id = "GOOD_MATCH_PERCENT")
      shinyjs::hide(id = "MAX_FEATURES")
    } else {
      shinyjs::show(id = "GOOD_MATCH_PERCENT")
      shinyjs::show(id = "MAX_FEATURES")
      if(grepl("Non-Rigid", input$Method)){
        updateSelectInput(session, "Method", selected = "Homography") 
        showNotification("Brute-Force Matching can't be used with Non-Rigid Registration\n")
      }
    }
  })
  
  observeEvent(input$Matcher, {
    if(grepl("FLANN", input$Matcher)){
      shinyjs::hide(id = "GOOD_MATCH_PERCENT")
      shinyjs::hide(id = "MAX_FEATURES")
    } else {
      shinyjs::show(id = "GOOD_MATCH_PERCENT")
      shinyjs::show(id = "MAX_FEATURES")
      if(grepl("Non-Rigid", input$Method)){
        updateSelectInput(session, "Method", selected = "Homography") 
        showNotification("Brute-Force Matching can't be used with Non-Rigid Registration\n")
      }
    }
  })
}

####
# Registering Objects ####
####

#' getRegisteredObject
#'
#' Get registered list of VoltRon objects
#'
#' @param obj_list a list of VoltRon objects
#' @param mapping_list a list of transformation matrices
#' @param register_ind the indices of query images/spatialdatasets
#' @param centre the index of the central reference image/spatialdata
#' @param input input
#' @param reg_mode the registration mode, either "auto" or "manual"
#' @param image_list the list of query/ref images
#' @param aligned_image_list the list of aligned query/ref images
#'
#' @noRd
getRegisteredObject <- function(obj_list, mapping_list, register_ind, centre, input, reg_mode = "manual", image_list = NULL, aligned_image_list = NULL){

  # initiate registered VoltRon objects
  ref_ind <- centre
  registered_sr <- list()

  # the original reference object
  registered_sr[[ref_ind]] <- obj_list[[ref_ind]]

  # waiter start
  withProgress(message = 'Register Coordinates (and Segments)', value = 0, {
  # waiter::waiter_show(html = waiter::spin_ring(), color = paste0("rgba(128,128,128,", 0.15, ")"))

  # register all assays
  for(i in register_ind){

    # choose image query and ref order
    if(i > ref_ind){
      ref_extension = paste0("ref_image",ref_ind)
      query_extension = paste0("query_image",i)
    } else {
      ref_extension = paste0("query_image",ref_ind)
      query_extension = paste0("ref_image",i)
    }

    # register the VoltRon object
    for(assy in vrAssayNames(obj_list[[i]], assay = "all")){

      # Increment the progress bar, and update the detail text.
      incProgress(1/length(register_ind), detail = paste("Register", assy, "of Layer", i, sep = " "))

      # register assay
      obj_list[[i]] <- applyPerspectiveTransform(obj_list[[i]],
                                                 assay = assy,
                                                 mapping = mapping_list[[paste0(i)]],
                                                 reference_image = image_list[[ref_ind]],
                                                 input = input,
                                                 reg_mode = reg_mode,
                                                 ref_extension = ref_extension,
                                                 query_extension = query_extension)

    }
    registered_sr[[i]] <- obj_list[[i]]

  }

  # waiter end
  # waiter::waiter_hide()

  })
  return(registered_sr)
}

#' applyPerspectiveTransform
#'
#' Applying a perspective transformation to the VoltRon object
#'
#' @param object a VoltRon objects
#' @param mapping a list of transformation matrices
#' @param reference_image the reference image
#' @param input input
#' @param reg_mode the registration mode, either "auto" or "manual"
#' @param ref_extension the shiny extension of reference image
#' @param query_extension the shiny extension of query image
#'
#' @importFrom magick image_info
#'
#' @noRd
applyPerspectiveTransform <- function(object,
                                      assay = NULL,
                                      mapping,
                                      reference_image,
                                      input,
                                      reg_mode,
                                      ref_extension,
                                      query_extension){

  # check assay
  if(is.null(assay))
    assay <- vrAssayNames(object)

  # get coordinates, segments and spatial points
  coords <- vrCoordinates(object, assay = assay)
  segments <- vrSegments(object, assay = assay)

  if(reg_mode == "manual"){

    # get the multiplication of all homography matrices
    # cur_mapping <- Reduce("%*%", mapping)
    mapping <- manageMapping(mapping)
    
    # get registered coordinates
    coords_reg <- coords
    coords_reg[,c("x", "y")] <- applyTransform(coords[,c("x", "y")], mapping)
    rownames(coords_reg) <- rownames(coords)
    colnames(coords_reg) <- colnames(coords)

    # get registered segments
    if(length(segments) > 0){
      segments_reg <- do.call(rbind, segments)
      segments_reg <- as.matrix(segments_reg)
      segments_reg[,colnames(segments_reg) %in% c("x", "y")] <- applyTransform(segments_reg[,colnames(segments_reg) %in% c("x", "y")], mapping)
      segments_reg <- as.data.frame(segments_reg)
      segments_reg <- split(segments_reg, segments_reg[,1])
      names(segments_reg) <- names(segments)
    } else {
      segments_reg <- segments
    }

    # get registered image (including all channels)
    image_reg_list <- sapply(vrImageChannelNames(object[[assay]]), function(x) NULL, USE.NAMES = TRUE)
    for(channel_ind in names(image_reg_list)){
      query_image <- vrImages(object[[assay]], channel = channel_ind, as.raster = TRUE)
      if(!inherits(query_image, "Image_Array")){
        query_image <- magick::image_read(query_image)
      }
      warped_image <- getRcppWarpImage(ref_image = reference_image,
                                       # query_image = vrImages(object, assay = assay, channel = channel_ind),
                                       query_image = query_image,
                                       mapping = mapping)
      image_reg_list[[channel_ind]] <- warped_image
    }

  } else if(reg_mode == "auto"){

    # get the multiplication of all homography matrices
    # cur_mapping <- Reduce("%*%", mapping)
    mapping <- manageMapping(mapping)
    
    # images
    ref_image <- transformImage(reference_image, ref_extension, input)
    # query_image <- transformImage(vrImages(object, assay = assay), query_extension, input)
    query_image <- vrImages(object[[assay]], as.raster = TRUE)
    if(!inherits(query_image, "Image_Array")){
      query_image <- magick::image_read(query_image)
    }
    query_image <- transformImage(query_image, query_extension, input)

    # image info
    # query_info <- magick::image_info(query_image)
    # ref_info <- magick::image_info(ref_image)
    query_info <- getImageInfo(query_image)
    ref_info <- getImageInfo(ref_image)

    # get registered coordinates
    coords_reg <- as.data.frame(coords)
    coords_reg <- transformImageKeypoints(query_image, coords_reg[,c("x","y")], query_extension, input)$keypoints

    coords_reg[,2] <- query_info$height - coords_reg[,2]
    coords_reg <- as.matrix(coords_reg)
    coords_reg <- applyTransform(coords_reg, mapping)
    coords_reg <- as.data.frame(coords_reg)
    coords_reg[,2] <- ref_info$height - coords_reg[,2]

    colnames(coords_reg) <- c("x", "y")
    coords_reg <- transformKeypoints(ref_image, coords_reg, ref_extension, input)
    coords_reg <- as.matrix(coords_reg)
    rownames(coords_reg) <- rownames(coords)
    
    # fix 3rd dimension 
    coords[,c("x", "y")] <- coords_reg[,c("x", "y")]
    coords_reg <- coords

    # get registered segments
    if(length(segments) > 0){
      segments_reg <- do.call(rbind, segments)
      segments_reg <- as.data.frame(segments_reg)
      segments_reg <- transformImageKeypoints(query_image, segments_reg, query_extension, input)$keypoints

      segments_reg[,colnames(segments_reg) %in% c("y")] <- query_info$height - segments_reg[,colnames(segments_reg) %in% c("y")]
      segments_reg <- as.matrix(segments_reg)
      segments_reg[,colnames(segments_reg) %in% c("x", "y")] <- applyTransform(segments_reg[,colnames(segments_reg) %in% c("x", "y")], mapping)
      segments_reg <- as.data.frame(segments_reg)
      segments_reg[,colnames(segments_reg) %in% c("y")]  <- ref_info$height - segments_reg[,colnames(segments_reg) %in% c("y")]

      segments_reg <- transformKeypoints(ref_image, segments_reg, ref_extension, input)

      segments_reg <- split(segments_reg, segments_reg[,1])
      names(segments_reg) <- names(segments)
    } else {
      segments_reg <- segments
    }

    # get registered image (including all channels)
    image_reg_list <- sapply(vrImageChannelNames(object[[assay]]), function(x) NULL, USE.NAMES = TRUE)
    for(channel_ind in names(image_reg_list)){

      # rotate, flip and flop before warping in C++
      ref_image <- transformImage(reference_image, ref_extension, input)
      # query_image <- transformImage(vrImages(object, assay = assay, channel = channel_ind),
      #                               query_extension, input)
      query_image <- vrImages(object[[assay]], channel = channel_ind, as.raster = TRUE)
      if(!inherits(query_image, "Image_Array")){
        query_image <- magick::image_read(query_image)
      }
      query_image <- transformImage(query_image, query_extension, input)
      query_image <- getRcppWarpImage(ref_image, query_image, mapping = mapping)
      query_image <- transformImageReverse(query_image, ref_extension, input)

      image_reg_list[[channel_ind]] <- query_image
    }
  }

  # make new image object
  vrImages(object[[assay]], reg = TRUE) <- formImage(coords = coords_reg, segments = segments_reg, image = image_reg_list)
  
  # set up the spatial coordinate name 
  vrMainSpatial(object[[assay]]) <- paste0(vrMainSpatial(object[[assay]]), "_reg")

  # return object
  return(object)
}

####
# Managing Mappings ####
####

manageMapping <- function(mappings){
  
  # check if all transformations are homography
  allHomography <- suppressWarnings(all(lapply(mappings, function(map){
    nrow(map[[1]] > 0) && is.null(map[[2]])
  })))
  
  # change the mapping
  new_mappings <- list()
  if(allHomography){
    mappings <- lapply(mappings, function(map) map[[1]])
    new_mappings <- list(
      list(Reduce("%*%", mappings),
           NULL)
    )
  } else {
    new_mappings <- mappings
  }
  
    
  # return
  return(new_mappings)
}
  
####
# Managing Keypoints ####
####

#' initateKeypoints
#'
#' Initiate shiny reactive values for keypoint dataframes for pairwise reference and query images
#'
#' @param len_images the length of images
#' @param keypoints_list the list of keypoint pairs
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @importFrom dplyr tibble
#'
#' @noRd
initateKeypoints <- function(len_images, keypoints_list, input, output, session){

  # initiate keypoints
  if(is.null(keypoints_list)){
    keypoints_list <- lapply(1:(len_images-1), function(i) {
      list(ref = dplyr::tibble(KeyPoint = numeric(), x = numeric(), y = numeric()),
           query = dplyr::tibble(KeyPoint = numeric(), x = numeric(), y = numeric()))

    })

    # set names for keypoints
    names(keypoints_list) <- paste0(1:(len_images-1),"-",2:len_images)
  }

  # return keypoints as reactive values
  do.call("reactiveValues", keypoints_list)
}

#' manageKeypoints
#'
#' A list of shiny observe events for tables and auxiliary operations for pairwise reference and query image
#'
#' @param centre center image index
#' @param register_ind query image indices
#' @param xyTable_list a list of keypoints x,y coordinates for each magick image
#' @param image_list a lost of magick image
#' @param info_list a list of magick image info on width and height
#' @param zoom_list a list of x,y ranges of query and ref images
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @noRd
manageKeypoints <- function(centre, register_ind, xyTable_list, image_list, info_list, zoom_list, input, output, session){

  # get image types
  image_types <- c("ref","query")

  # get the length of tables
  len_tables <- length(xyTable_list)

  # set click operations for reference and query points
  lapply(1:len_tables, function(i){
    lapply(image_types, function(type){

      # listen to click operations for reference/query plots
      observeEvent(input[[paste0("click_plot_", type ,i)]], {

        # get brush information
        brush <- input[[paste0("brush_plot_", type ,i)]]
        limits <- cbind(zoom_list[[paste0(i)]][[type]]$x, zoom_list[[paste0(i)]][[type]]$y)
        if (is.null(brush)) {

          # get image
          image <- image_list[[i]]

          # get and transform keypoints
          keypoint <- data.frame(x = input[[paste0("click_plot_",type,i)]]$x,
                                 y = input[[paste0("click_plot_",type,i)]]$y)
          
          # get the transformed zoom info first and calculate width, then record transformed image
          limits_trans <- data.frame(x = limits[,1], y = limits[,2])
          limits_trans <- transformImageKeypoints(image, limits_trans, paste0(type, "_image",i), input)
          image_trans <- limits_trans$image
          limits_trans <- data.frame(x = range(limits_trans$keypoints[,1]), y = range(limits_trans$keypoints[,2]))
          
          # correct for scaling, scale factor = 1000
          width <- limits_trans[2,1]-limits_trans[1,1]
          height <- limits_trans[2,2]-limits_trans[1,2]
          if(max(height,width) > 1000){
            if(inherits(image_trans, "Image_Array")){
              n.series <- ImageArray::len(image_trans)
              cur_width <- width
              cur_height <- height
              for(ii in 2:n.series){
                cur_width <- width/(2^(ii-1))
                cur_height <- height/(2^(ii-1))
                if(max(cur_height, cur_width) <= 1000){
                  break
                }
              }
              keypoint <- keypoint*width/ceiling(cur_width)
            } else {
              keypoint <- keypoint*width/1000 
            }
          }
          
          # correct for zoom information
          keypoint <- keypoint + limits_trans[1,]
          
          # correct for flipflop and rotate
          keypoint <- transformKeypoints(image_trans, keypoint, paste0(type, "_image",i), input)
          
          # insert keypoint to associated table
          ref_ind <- ifelse(type == "ref", i, i-1) # select reference image
          
          # insert keypoint to associated table
          temp <- xyTable_list[[paste0(ref_ind, "-", ref_ind+1)]][[type]]
          temp <- temp %>%
            add_row(KeyPoint = nrow(temp)+1, x = keypoint$x, y = keypoint$y)
          xyTable_list[[paste0(ref_ind, "-", ref_ind+1)]][[type]] <- temp

        }
      })
    })
  })

  # remove keypoints from images
  lapply(1:len_tables, function(i){
    lapply(image_types, function(type){
      observeEvent(input[[paste0("remove_", type, i)]], {
        ref_ind <- ifelse(type == "ref", i, i-1) # select reference image
        temp <- xyTable_list[[paste0(ref_ind, "-", ref_ind+1)]][[type]]
        if(nrow(temp) > 0){
          temp <- temp %>% filter(KeyPoint != nrow(temp))
          xyTable_list[[paste0(ref_ind, "-", ref_ind+1)]][[type]] <- temp
        }
      })
    })
  })
}

#' transformImageKeypoints
#'
#' Apply given transformations to a magick image and keypoints for plotting
#'
#' @param image magick image
#' @param extension name extension for the shiny input parameter
#' @param keypoints a set of keypoints
#' @param input shiny input
#' @param session shiny session
#'
#' @importFrom magick image_negate image_rotate image_flip image_flop image_info
#'
#' @noRd
transformImageKeypoints <- function(image, keypoints, extension, input, session){

  if(is.null(keypoints))
    return(list(image = image, keypoints = keypoints))

  # negate image
  input_negate <- input[[paste0("negate_", extension)]]
  if(input_negate == "Yes"){
    # image <- magick::image_negate(image)
    image <- negateImage(image)
  }

  # get unrotated image info
  # image_limits <- unlist(magick::image_info(image)[1,c("width", "height")])
  image_limits <- unlist(getImageInfo(image)[1,c("width", "height")])
  image_origin <- image_limits/2

  # rotate image and keypoints
  input_rotate <- as.numeric(input[[paste0("rotate_", extension)]])
  # image <- magick::image_rotate(image, input_rotate)
  image <- rotateImage(image, input_rotate)

  # get rotated image info
  # rotated_image_limits <- unlist(magick::image_info(image)[1,c("width", "height")])
  rotated_image_limits <- unlist(getImageInfo(image)[1,c("width", "height")])
  rotated_image_origin <- rotated_image_limits/2

  # rotate keypoints
  keypoints <- rotateKeypoint(keypoints, input_rotate, image_origin, image_limits, rotated_image_origin, rotated_image_limits)

  # flip flop image and keypoints
  input_flipflop <- input[[paste0("flipflop_", extension)]]
  if(input_flipflop == "Flip"){
    # image <- magick::image_flip(image)
    image <- flipImage(image)
  } else if(input_flipflop == "Flop"){
    # image <- magick::image_flop(image)
    image <- flopImage(image)
  }

  # flipflop keypoints
  keypoints <- flipflopKeypoint(keypoints, rotated_image_limits, input_flipflop)

  # return both the image and the keypoints
  return(list(image = image, keypoints = keypoints))
}

#' transformKeypoints
#'
#' Apply transformations to keypoints given transformed images to find the keypoints locations in the original image
#'
#' @param image magick image
#' @param keypoints keypoints visualized on image
#' @param extension name extension for the shiny input parameter
#' @param input shiny input
#'
#' @importFrom magick image_flip image_flop image_rotate
#'
#' @noRd
transformKeypoints <- function(image, keypoints, extension, input){

  # get unrotated image info
  # image_limits <- unlist(image_info(image)[1,c("width", "height")])
  image_limits <- unlist(getImageInfo(image)[1,c("width", "height")])
  image_origin <- image_limits/2

  # flip flop image and keypoints
  input_flipflop <- input[[paste0("flipflop_", extension)]]
  if(input_flipflop == "Flip"){
    # image <- magick::image_flip(image)
    image <- flipImage(image)
  } else if(input_flipflop == "Flop"){
    # image <- magick::image_flop(image)
    image <- flopImage(image)
  }
  keypoints <- flipflopKeypoint(keypoints, image_limits, input_flipflop)

  # rotate image (reverse) and keypoints
  input_rotate <- 360 - as.numeric(input[[paste0("rotate_", extension)]])
  # image <- magick::image_rotate(image, input_rotate)
  image <- rotateImage(image, input_rotate)

  # get rotated image info
  # rotated_image_limits <- unlist(image_info(image)[1,c("width", "height")])
  rotated_image_limits <- unlist(getImageInfo(image)[1,c("width", "height")])
  rotated_image_origin <- rotated_image_limits/2

  # rotate keypoints
  keypoints <- rotateKeypoint(keypoints, input_rotate, image_origin, image_limits, rotated_image_origin, rotated_image_limits)

  return(keypoints)
}

#' rotateKeypoint
#'
#' Find transformations of keypoints under clockwise rotations of the image
#'
#' @param keypoints dataset of keypoints
#' @param angle angle of rotation [0,360]
#' @param origin center of the image
#' @param limits limits of the image
#' @param rotated_origin center of the rotated image
#' @param rotated_limits limits of the rotated image
#'
#' @noRd
rotateKeypoint <- function(keypoints, angle, origin, limits, rotated_origin, rotated_limits){

  # if there are no points, return
  if(nrow(keypoints) == 0)
    return(keypoints)

  # get coordinates from the keypoints dataset
  points <- keypoints[,c("x","y")]

  # set rotation matrix for angles
  radii <- ((360-angle)*pi/180)
  s = sin(radii);
  c = cos(radii);
  rotation_mat <- matrix(c(c, s, -s, c), nrow = 2, byrow = F)

  # rotate point
  # points <- t(apply(points, 1, function(x) return(x - origin)))
  # points <- t(apply(points, 1, function(x) return(x/limits)))
  # rotated_points <- t(rotation_mat %*% t(points))
  # rotated_points <- t(apply(rotated_points, 1, function(x) return(x*rotated_limits)))
  # rotated_points <- t(apply(rotated_points, 1, function(x) return(x + rotated_origin)))
  points <- points - matrix(rep(origin, nrow(points)), nrow = nrow(points), byrow = T)
  points <- points * matrix(rep(1/limits, nrow(points)), nrow = nrow(points), byrow = T)
  rotated_points <- t(rotation_mat %*% t(points))
  rotated_points <- rotated_points * matrix(rep(rotated_limits, nrow(points)), nrow = nrow(rotated_points), byrow = T)
  rotated_points <- rotated_points + matrix(rep(rotated_origin, nrow(points)), nrow = nrow(rotated_points), byrow = T)

  # put rotated points back to keypoints
  keypoints[,c("x","y")] <- rotated_points

  return(keypoints)
}

#' flipflopKeypoint
#'
#' Find transformed keypoints on image given any flip or flop action by magick
#'
#' @param keypoints dataset of keypoints
#' @param image_limits limits of the images
#' @param flipflop a flip or flop action as string
#'
#' @noRd
flipflopKeypoint <- function(keypoints, image_limits, flipflop){

  if(nrow(keypoints) == 0)
    return(keypoints)

  if(grepl("Flop", flipflop))
    keypoints$x = image_limits[1] - keypoints$x

  if(grepl("Flip", flipflop))
    keypoints$y = image_limits[2] - keypoints$y

  return(keypoints)
}

#' imageKeypoint
#'
#' add keypoints as points on ggplot object
#'
#' @param image magick image
#' @param keypoints keypoints to draw on image
#'
#' @noRd
imageKeypoint <- function(image, keypoints){

  if(is.null(keypoints))
    return(image)

  # select keypoints and texts on image
  image <- image +
    geom_point(mapping = aes(x = x, y = y), keypoints, size = 8, shape = 21, fill = "white") +
    geom_text(mapping = aes(x = x, y = y, label = KeyPoint), keypoints, size = 5)
}

####
# Managing Zoom Options ####
####

#' imageZoom
#'
#' zoom image
#'
#' @param image magick image
#' @param zoom_info zoom info to draw on image
#'
#' @importFrom magick image_info
#' 
#' @noRd
imageZoom <- function(image, zoom_info = NULL){
  
  if(is.null(zoom_info))
    return(image)
  
  # get image info
  # imageinfo <- magick::image_info(image)
  imageinfo <- getImageInfo(image)
  
  # get info of zoom
  zoom_info <- FromBoxToCrop(as.data.frame(zoom_info), imageinfo)
  
  # return
  return(zoom_info)
}

#' initiateZoomOptions
#'
#' Initiate shiny reactive values for capturing zoom/brush limits
#'
#' @param info_list the list of image information
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @noRd
initiateZoomOptions <- function(info_list, input, output, session){
  
  # length of images
  len_images <- length(info_list)
  
  # initiate zoom options list
  zoom_list <- lapply(1:len_images, function(i) {
    list(ref = list(x = c(0, info_list[[i]][1]), y = c(0, info_list[[i]][2])),
         query = list(x = c(0, info_list[[i]][1]), y = c(0, info_list[[i]][2])))
  })
  
  # set names for keypoints
  names(zoom_list) <- paste0(1:len_images)
  
  # return keypoints as reactive values
  do.call("reactiveValues", zoom_list)
}

#' manageImageZoomOptions
#'
#' A list of shiny observe events for handling zoom options of image outputs
#'
#' @param centre center image index
#' @param register_ind query image indices
#' @param zoom_list a list of x,y ranges of query and ref images
#' @param image_list a list of transformed magick image
#' @param info_list the list of image information
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @noRd
manageImageZoomOptions <- function(centre, register_ind, zoom_list, image_list, info_list, input, output, session){

  # get image types
  image_types <- c("ref","query")
  
  # get the length of tables
  len_tables <- length(zoom_list)
  
  # set click operations for reference and query points
  lapply(1:len_tables, function(i){
    lapply(image_types, function(type){
      
      # listen to click operations for reference/query plots
      observeEvent(input[[paste0("dblclick_plot_", type ,i)]], {
        
        # get brush information
        brush <- input[[paste0("brush_plot_", type ,i)]]
        limits <- cbind(zoom_list[[paste0(i)]][[type]]$x, zoom_list[[paste0(i)]][[type]]$y)
        if (!is.null(brush)) {
          
          # get brush variables
          brush_mat <- data.frame(x = c(brush$xmin, brush$xmax), 
                                  y = c(brush$ymin, brush$ymax))

          # get image
          image <- image_list[[i]]
          
          # get the transformed limits first and calculate width, then record transformed image
          limits_trans <- data.frame(x = limits[,1], y = limits[,2])
          limits_trans <- transformImageKeypoints(image, limits_trans, paste0(type, "_image",i), input)
          image_trans <- limits_trans$image
          limits_trans <- data.frame(x = range(limits_trans$keypoints[,1]), y = range(limits_trans$keypoints[,2]))

          # if width is large, then correct the brush event for the downsize effect
          width <- limits_trans[2,1]-limits_trans[1,1]
          height <- limits_trans[2,2]-limits_trans[1,2]
          if(max(height,width) > 1000){
            if(inherits(image_trans, "Image_Array")){
              n.series <- ImageArray::len(image_trans)
              cur_width <- width
              cur_height <- height
              for(ii in 2:n.series){
                cur_width <- width/(2^(ii-1))
                cur_height <- height/(2^(ii-1))
                if(max(cur_height, cur_width) <= 1000){
                  break
                }
              }
              brush_mat <- brush_mat*width/ceiling(cur_width)
            } else {
              brush_mat <- brush_mat*width/1000
            }
          }

          # correct brush for the zoom effect
          brush_mat[,1] <- brush_mat[,1] + limits_trans[1,1]
          brush_mat[,2] <- brush_mat[,2] + limits_trans[1,2]

          # correct for flipflop and rotate using the transformed image from above
          brush_mat <- transformKeypoints(image_trans, as.data.frame(brush_mat), paste0(type, "_image",i), input)
          brush_mat <- data.frame(x = range(brush_mat[,1]), y = range(brush_mat[,2]))
          brush_mat <- as.matrix(brush_mat)

          # make new zoom information
          zoom_list[[paste0(i)]][[type]]$x <- brush_mat[,1]
          zoom_list[[paste0(i)]][[type]]$y <- brush_mat[,2]
        
        } else {
          zoom_list[[paste0(i)]][[type]]$x <- c(0, info_list[[i]][1])
          zoom_list[[paste0(i)]][[type]]$y <- c(0, info_list[[i]][2])
        }
      })
    })
  })
}

####
# Managing Images ####
####

#' getImageOutput
#'
#' Shiny outputs for a set of magick images with keypoints
#'
#' @param image_list a list of magick images
#' @param info_list a list of magick image info on width and height
#' @param keypoints_list a list of data frames, each having a set of keypoints
#' @param zoom_list a list of x,y ranges of query and ref images
#' @param centre the center image index
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @importFrom magick image_ggplot image_resize image_crop geometry_size_percent
#'
#' @noRd
getImageOutput <- function(image_list, info_list, keypoints_list = NULL, zoom_list, centre, input, output, session){

  # get image types
  image_types <- c("ref","query")

  # get the length of images
  len_images <- length(image_list)

  # output query images
  lapply(1:len_images, function(i){
    lapply(image_types, function(type){

      # image output
      output[[paste0("plot_", type, i)]] <- renderPlot({

        # select keypoints
        ref_ind <- ifelse(type == "ref", i, i-1) # select reference image
        keypoints <- keypoints_list[[paste0(ref_ind, "-", ref_ind+1)]][[type]]
        
        # transform image and keypoints
        img <- image_list[[i]][[input[[paste0("channel_", type, "_image", i)]]]]
        img_trans <- transformImageKeypoints(img, keypoints, paste0(type, "_image",i), input, session)
        
        # zoom images and keypoints
        limits <- as.data.frame(zoom_list[[paste0(i)]][[type]])
        img_limits <- transformImageKeypoints(img, limits, paste0(type, "_image",i), input, session)
        img_limits$keypoints <- data.frame(x = range(img_limits$keypoints[,1]), y = range(img_limits$keypoints[,2]))
        imgzoom <- imageZoom(img_trans$image, zoom_info = img_limits$keypoints)
        if(!is.null(img_trans$keypoints)){
          if(nrow(img_trans$keypoints) > 0){
            temp <- as.matrix(img_trans$keypoints[,c("x","y")])
            temp <- temp - matrix(unlist(rep(img_limits$keypoints[1,], nrow(img_trans$keypoints))), nrow = nrow(img_trans$keypoints), byrow = T)
            img_trans$keypoints[,c("x","y")] <- temp
          }
        }
        
        # img_trans$image <- magick::image_crop(img_trans$image, geometry = imgzoom)
        img_trans$image <- cropImage(img_trans$image, geometry = imgzoom)
        
        # lower resolution
        width <- img_limits$keypoints[2,1]-img_limits$keypoints[1,1]
        height <- img_limits$keypoints[2,2]-img_limits$keypoints[1,2]
        if(max(height, width) > 1000){
          
          # scale keypoints
          if(inherits(img_trans$image, "Image_Array")){
            n.series <- ImageArray::len(img_trans$image)
            cur_width <- width
            cur_height <- height
            for(ii in 2:n.series){
              cur_width <- width/(2^(ii-1))
              cur_height <- height/(2^(ii-1))
              if(max(cur_height, cur_width) <= 1000){
                break
              }
            }
            img_trans$keypoints[,c("x","y")] <- img_trans$keypoints[,c("x","y")]*(cur_width/width)
          } else {
            img_trans$keypoints[,c("x","y")] <- img_trans$keypoints[,c("x","y")]*(1000/width)
          }
        } 


        # visualize
        # imgggplot <- magick::image_ggplot(img_trans$image)
        img_ggplot <- plotImage(img_trans$image, max.pixel.size = 1000)
        img_ggplot <- imageKeypoint(img_ggplot, img_trans$keypoints)
        
        # return
        return(img_ggplot)
      })

      # update info
      output[[paste0("scaleinfo_", type, "_image", i)]] <- renderText({
        cur_info <- info_list[[i]] * input[[paste0("scale_", type, "_image", i)]]
        paste(cur_info, collapse = "x")
      })
    })
  })
}

#' plotImage
#'
#' plot image
#'
#' @param image a magick image or DelayedArray object
#' 
#' @importFrom magick image_ggplot
#'
#' @noRd
plotImage <- function(image, max.pixel.size = NULL){
  
  if(inherits(image, "magick-image")){
    imageinfo <- getImageInfo(image)
    if(max(imageinfo$width, imageinfo$height) > max.pixel.size){
      image <- magick::image_resize(image, geometry = as.character(max.pixel.size))
    }
    imgggplot <- magick::image_ggplot(image)
  } else if(inherits(image, "Image_Array")){
    img_array <- as.array(image, max.pixel.size = max.pixel.size)
    img_raster <- as.raster_array(aperm(img_array, perm = c(3,2,1)), max = 255)
    info <- list(width = dim(img_raster)[2], height = dim(img_raster)[1])
    imgggplot <- ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes_string("x", "y")) + 
      ggplot2::geom_blank() + 
      ggplot2::theme_void() + 
      ggplot2::coord_fixed(expand = FALSE, 
                           xlim = c(0, info$width), 
                           ylim = c(0, info$height)) + 
      ggplot2::annotation_raster(img_raster, 0, info$width, info$height, 0, interpolate = FALSE)
  }
  imgggplot
}

#' getImageInfoList
#'
#' get information on list of images
#'
#' @param image_list a list of magick images or DelayedArray objects
#'
#' @noRd
getImageInfoList <- function(image_list){
  lapply(image_list, function(x){
    # imginfo <- magick::image_info(x)
    imginfo <- getImageInfo(x)
    c(imginfo$width, imginfo$height)
  })
}

#' getImageInfo
#'
#' get information on images
#'
#' @param image a magick image or DelayedArray object
#'
#' @importFrom magick image_info
#'
#' @noRd
getImageInfo <- function(image){
  
  if(inherits(image, "magick-image")){
    imginfo <- magick::image_info(image)
  } else if(inherits(image, "Image_Array")){
    # dim_image <- dim(image)
    # imginfo <- list(width = dim_image[2], height = dim_image[3])
    imginfo <- ImageArray::getImageInfo(image)
  }
  as.data.frame(imginfo)
}

#' rotateImage
#'
#' rotate images
#'
#' @param image a magick image or DelayedArray object
#' @param degrees value between 0 and 360 for how many degrees to rotate
#'
#' @importFrom magick image_rotate
#'
#' @noRd
rotateImage <- function(image, degrees){
  
  if(inherits(image, "magick-image")){
    image <- magick::image_rotate(image, degrees = degrees)
  } else if(inherits(image, "Image_Array")){
    image <- ImageArray::rotate(image, degrees)
  }
  image
}

#' negateImage
#'
#' negate images
#'
#' @param image a magick image or DelayedArray object
#' 
#' @importFrom magick image_negate
#'
#' @noRd
negateImage <- function(image){
  
  if(inherits(image, "magick-image")){
    image <- magick::image_negate(image)
  } else if(inherits(image, "Image_Array")){
    # image <- 255 - image
    image <- ImageArray::negate(image)
  }
  image
}

#' flipImage
#'
#' flip images
#'
#' @param image a magick image or DelayedArray object
#' 
#' @importFrom magick image_negate
#'
#' @noRd
flipImage <- function(image){
  
  if(inherits(image, "magick-image")){
    image <- magick::image_flip(image)
  } else if(inherits(image, "Image_Array")){
    # dim_img <- dim(image)
    # image <- image[ , , dim_img[3]:1, drop = FALSE]
    image <- ImageArray::flip(image)
  }
  image
}

#' flopImage
#'
#' flop images
#'
#' @param image a magick image or DelayedArray object
#' 
#' @importFrom magick image_negate
#'
#' @noRd
flopImage <- function(image){
  
  if(inherits(image, "magick-image")){
    image <- magick::image_flop(image)
  } else if(inherits(image, "Image_Array")){
    # dim_img <- dim(image)
    # image <- image[ , dim_img[2]:1, , drop = FALSE]
    image <- ImageArray::flop(image)
  }
  image
}

#' cropImage
#'
#' crop images
#'
#' @param image a magick image or DelayedArray object
#' @param geometry a geometry string specifying area (for cropping) or size (for resizing).
#' 
#' @importFrom magick image_crop
#'
#' @noRd
cropImage <- function(image, geometry){
  
  if(inherits(image, "magick-image")){
    image <- magick::image_crop(image, geometry = geometry)
  } else if(inherits(image, "Image_Array")){
    crop_info_int <- as.integer(strsplit(geometry, split = "[x|+]")[[1]])
    # image <- image[,crop_info_int[3]:(crop_info_int[3]+crop_info_int[1]), crop_info_int[4]:(crop_info_int[4]+crop_info_int[2]), drop = FALSE]
    image <- ImageArray::crop(image, ind = list(crop_info_int[3]:(crop_info_int[3]+crop_info_int[1]), crop_info_int[4]:(crop_info_int[4]+crop_info_int[2])))
  }
  image
}

#' resizeImage
#'
#' resize images
#'
#' @param image a magick image or DelayedArray object
#' @param geometry a geometry string specifying area (for cropping) or size (for resizing).
#' 
#' @importFrom magick image_resize image_info image_read geometry_size_percent
#'
#' @noRd
resize_Image <- function(image, geometry){
  
  # get image info
  image_info_large <- getImageInfo(image)
  
  if(inherits(image, "magick-image")){
    image <- magick::image_resize(image, geometry = geometry)
  } else if(inherits(image, "Image_Array")){
    
    # get scale factor 
    if(grepl("%$", geometry)){
      scale_factor <- as.numeric(gsub("%$", "", geometry))/100
    } else if(grepl("x$", geometry)){
      scale_factor <- (as.numeric(gsub("x$", "", geometry))/image_info_large$width)
    }

    # get scaled array 
    scaled_image_info <- ceiling(image_info_large*scale_factor)
    image <- as.array(image, min.pixel.size = max(scaled_image_info))
    
    # convert to magick image
    image <- magick::image_read(array(as.raw(image), dim = dim(image)))
    image_info <- magick::image_info(image)
    image <- magick::image_resize(image, geometry = geometry_size_percent(100*scaled_image_info[1]/image_info$width))
  }
  image
}

#' transformImage
#'
#' Apply given transformations to a magick image
#'
#' @param image magick image
#' @param extension name extension for the shiny input parameter
#' @param input shiny input
#'
#' @importFrom magick image_flip image_flop image_rotate
#'
#' @noRd
transformImage <- function(image, extension, input){

  # rotate image and keypoints
  input_rotate <- as.numeric(input[[paste0("rotate_", extension)]])
  # image <- magick::image_rotate(image, input_rotate)
  image <- rotateImage(image, input_rotate)
  
  # flip flop image and keypoints
  input_flipflop <- input[[paste0("flipflop_", extension)]]
  if(input_flipflop == "Flip"){
    # image <- magick::image_flip(image)
    image <- flipImage(image)
  } else if(input_flipflop == "Flop"){
    # image <- magick::image_flop(image)
    image <- flopImage(image)
  }

  # return image
  image
}

#' transformImageReverse
#'
#' Apply given transformations to a magick image in reverse fashion
#'
#' @param image magick image
#' @param extension name extension for the shiny input parameter
#' @param input shiny input
#'
#' @importFrom magick image_flip image_flop image_rotate
#'
#' @noRd
transformImageReverse <- function(image, extension, input){

  # flip flop image and keypoints
  input_flipflop <- input[[paste0("flipflop_", extension)]]
  if(input_flipflop == "Flip"){
    # image <- magick::image_flip(image)
    image <- flipImage(image)
  } else if(input_flipflop == "Flop"){
    # image <- magick::image_flop(image)
    image <- flopImage(image)
  }

  # rotate image and keypoints
  input_rotate <- 360 - as.numeric(input[[paste0("rotate_", extension)]])
  # image <- magick::image_rotate(image, input_rotate)
  image <- rotateImage(image, input_rotate)

  # return image
  image
}

#' transformImageQueryList
#'
#' Apply given transformations to a list of magick image and return shiny reactive
#'
#' @param image_list magick image
#' @param input shiny input
#'
#' @noRd
transformImageQueryList <- function(image_list, input){

  # length of images
  len_register <- length(image_list) - 1

  trans_query_list <- lapply(1:len_register, function(i){
    reactive({
      list(ref = transformImage(image_list[[i]], paste0("ref_image",i), input),
           query = transformImage(image_list[[i+1]], paste0("query_image",i+1), input))
    })
  })

  ####
  names(trans_query_list) <- paste0(1:(length(image_list)-1),"-",2:length(image_list)) # REMOVE LATER, or decide not to
  ####

  return(trans_query_list)
}

#' getRcppWarpImage
#'
#' Warping a query image given a homography image
#'
#' @param ref_image reference image
#' @param query_image query image
#' @param mapping a list of the homography matrices and TPS keypoints 
#'
#' @importFrom magick image_read image_data
#' 
#' @export
getRcppWarpImage <- function(ref_image, query_image, mapping){
  
  # ref image
  if(inherits(ref_image, "Image_Array")){
    ref_image <- as.array(ref_image)
    ref_image <- array(as.raw(ref_image), dim = dim(ref_image))
  } else {
    ref_image <- magick::image_data(ref_image, channels = "rgb")
  }
  
  # query image
  if(inherits(query_image, "Image_Array")){
    query_image <- as.array(query_image)
    query_image <- array(as.raw(query_image), dim = dim(query_image))
  } else {
    query_image <- magick::image_data(query_image, channels = "rgb")
  }
  
  # warp image
  query_image <- warpImage(ref_image = ref_image, 
                           query_image = query_image, 
                           mapping = mapping,
                           width1 = dim(ref_image)[2], height1 = dim(ref_image)[3],
                           width2 = dim(query_image)[2], height2 = dim(query_image)[3])
  magick::image_read(query_image)
}

####
# Manual Image Registration ####
####

#' initiateMappings
#'
#' Initiate shiny reactive values for registration matrices
#'
#' @param len_images the number of query images
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @noRd
initiateMappings <- function(len_images, input, output, session){

  # initiate matrices
  matrix_list <- lapply(1:len_images, function(i) return(NULL))
  names(matrix_list) <- 1:len_images

  # return matrices as reactive values
  do.call("reactiveValues", matrix_list)
}

#' getManualRegisteration
#'
#' Manual registration of images using manually entered keypoints
#'
#' @param registration_mapping_list a list of mapping matrices used for registering VoltRon objects
#' @param spatdata_list a list of Spatial data object of the query images
#' @param image_list the list of query images
#' @param keypoints_list a list of keypoints x,y coordinates for query image
#' @param centre center image index
#' @param register_ind query image indices
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @import ggplot2
#' @importFrom magick image_write image_join image_read image_resize
#' @importFrom shiny reactiveValuesToList
#'
#' @noRd
getManualRegisteration <- function(registration_mapping_list, spatdata_list, image_list, keypoints_list,
                                   centre, register_ind, input, output, session){

  # the number of registrations
  len_register <- length(image_list) - 1

  # Registration events
  observeEvent(input$register, {

    # get key points as list
    keypoints_list <- shiny::reactiveValuesToList(keypoints_list)
    
    # Manual Registration
    if(!input$automatictag){

      # waiter start
      withProgress(message = paste0('Manual Registration (TPS)'), value = 0, {
      # waiter::waiter_show(html = waiter::spin_ring(), color = paste0("rgba(128,128,128,", 0.15, ")"))

      # Check keypoints
      keypoints_check_flag <- sapply(keypoints_list, function(key_list){
        nrow(key_list$ref) > 0 | nrow(key_list$query) > 0
      })
      if(!all(unlist(keypoints_check_flag))){
        showNotification("Please select keypoints for all images\n")
        return(NULL)
      }

      keypoints_check_flag <- sapply(keypoints_list, function(key_list){
        nrow(key_list$ref) == nrow(key_list$query)
      })
      if(!all(unlist(keypoints_check_flag))){
        showNotification("The number of reference and query keypoints should be equal! \n")
        return(NULL)
      }

      # Register keypoints
      mapping_list <- list()
      aligned_image_list <- list()
      for(i in register_ind){

        # Increment the progress bar, and update the detail text.
        incProgress(1/length(register_ind), detail = paste("Registering Image", i, sep = " "))

        # get a sequential mapping between a query and reference image
        results <- computeManualPairwiseTransform(image_list, keypoints_list, query_ind = i, ref_ind = centre, input = input)

        # save transformation mapping
        registration_mapping_list[[paste0(i)]] <- results$mapping

        # save matches
        aligned_image_list[[i]] <- results$aligned_image
      }

      # waiter end
      # waiter::waiter_hide()
      })

      # Plot registered images
      lapply(register_ind, function(i){
        output[[paste0("plot_query_reg",i)]] <- renderImage({

          # get image list
          # image_view_list <- list(rep(magick::image_resize(image_list[[centre]], geometry = "400x"),5),
          #                         rep(magick::image_resize(aligned_image_list[[i]], geometry = "400x"),5))
          image_view_list <- list(rep(resize_Image(image_list[[centre]], geometry = "400x"),5),
                                  rep(resize_Image(aligned_image_list[[i]], geometry = "400x"),5))

          # make slide show
          image_view_list <- image_view_list %>%
            magick::image_join() %>%
            magick::image_write(tempfile(fileext='gif'), format = 'gif')
          list(src = image_view_list, contentType = "image/gif")
        }, deleteFile = TRUE)
      })

      # Output summary
      output[["summary"]] <- renderUI({
        str1 <- paste0(" Registration Summary:")
        str2 <- paste0("# of Images: ", length(image_list))
        str3 <- paste0("# of Registrations: ", len_register)
        all_str <- c(str1, str2, str3)
        shiny::HTML(paste(all_str, collapse = '<br/>'))
      })
    }
  })
}

#' computeManualPairwiseTransform
#'
#' Computing transformation matrix of manual registration
#'
#' @param image_list the list of images
#' @param keypoints_list the list of keypoint matrices
#' @param query_ind the index of the query image
#' @param ref_ind the index of the reference image
#' @param input input
#'
#' @noRd
computeManualPairwiseTransform <- function(image_list, keypoints_list, query_ind, ref_ind, input){

  # determine the number of transformation to map from query to the reference
  indices <- query_ind:ref_ind
  mapping_mat <- rep(indices,c(1,rep(2,length(indices)-2),1))
  mapping_mat <- matrix(mapping_mat,ncol=2,byrow=TRUE)

  # reference and target landmarks/keypoints
  mapping <- list()
  aligned_image <- image_list[[query_ind]]
  for(kk in 1:nrow(mapping_mat)){
    cur_map <- mapping_mat[kk,]
    ref_image <- image_list[[cur_map[2]]]
    if(which.min(cur_map) == 1){
      key_ind <- paste0(cur_map[1], "-", cur_map[2])
      keypoints <- keypoints_list[[key_ind]]
      target_landmark <- as.matrix(keypoints[["ref"]][,c("x","y")])
      reference_landmark <- as.matrix(keypoints[["query"]][,c("x","y")])
    } else {
      key_ind <- paste0(cur_map[2], "-", cur_map[1])
      keypoints <- keypoints_list[[key_ind]]
      reference_landmark <- as.matrix(keypoints[["ref"]][,c("x","y")])
      target_landmark <- as.matrix(keypoints[["query"]][,c("x","y")])
    }

    if(which.max(cur_map) == 1){
      ref_label = "ref"
      query_label = "query"
    } else {
      ref_label = "query"
      query_label = "ref"
    }

    # get registered image (including all channels)
    reg <- getRcppManualRegistration(aligned_image, ref_image, target_landmark, reference_landmark, 
                                         method = input$Method)
    
    # return transformation matrix and images
    mapping[[kk]] <- list(reg$transmat[[1]], 
                          list(reference = reg$transmat[[2]][[1]],
                               query = reg$transmat[[2]][[2]]))
    aligned_image <- reg$aligned_image
  }

  return(list(mapping = mapping, 
              aligned_image = aligned_image))
}

#' getRcppManualRegistration
#'
#' Manual registration workflow with Rcpp
#'
#' @param query_image query image
#' @param ref_image reference image
#' @param query_landmark query landmark points
#' @param reference_landmark refernece landmark points
#' @param method the automated registration method, either TPS or Homography+TPS
#'
#' @importFrom magick image_read image_data
#'
#' @export
getRcppManualRegistration <- function(query_image, ref_image, query_landmark, reference_landmark, 
                                      method = "TPS") {
  
  
  # ref_image_rast <- magick::image_data(ref_image, channels = "rgb")
  # query_image_rast <- magick::image_data(query_image, channels = "rgb")
  
  # ref image
  if(inherits(ref_image, "Image_Array")){
    ref_image <- as.array(ref_image)
    ref_image <- array(as.raw(ref_image), dim = dim(ref_image))
  } else {
    ref_image <- magick::image_data(ref_image, channels = "rgb")
  }
  
  # query image
  if(inherits(query_image, "Image_Array")){
    query_image <- as.array(query_image)
    query_image <- array(as.raw(query_image), dim = dim(query_image))
  } else {
    query_image <- magick::image_data(query_image, channels = "rgb")
  }
  
  reference_landmark[,2] <- dim(ref_image)[3] - reference_landmark[,2]
  query_landmark[,2] <- dim(query_image)[3] - query_landmark[,2]
  reg <- manual_registeration_rawvector(ref_image = ref_image, query_image = query_image,
                                        reference_landmark = reference_landmark, query_landmark = query_landmark,
                                        width1 = dim(ref_image)[2], height1 = dim(ref_image)[3],
                                        width2 = dim(query_image)[2], height2 = dim(query_image)[3], 
                                        method = method)
  return(list(transmat = reg[[1]], 
              aligned_image = magick::image_read(reg[[2]])))
}

####
# Automated Image Registration ####
####

#' getManualRegisteration
#'
#' Manual registeration of images using manually entered keypoints
#'
#' @param registration_mapping_list a list of mapping matrices used for registering VoltRon objects
#' @param spatdata_list a list of Spatial data object of the query images
#' @param image_list the list of query images
#' @param channel_names the list of channel names for each image
#' @param centre center image index
#' @param register_ind query image indices
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @importFrom magick image_info image_ggplot image_write image_join image_resize
#' @importFrom grid rasterGrob
#' @importFrom ggplot2 ggplot coord_fixed annotation_raster annotation_custom
#'
#' @noRd
getAutomatedRegisteration <- function(registration_mapping_list, spatdata_list, image_list, channel_names, centre, register_ind,
                                      input, output, session){

  # the number of registrations
  len_register <- length(image_list) - 1

  # Registration events
  observeEvent(input$register, {
    
    # Automated registration
    if(input$automatictag){

      # waiter start
      withProgress(message = paste0('Automated Registration (', input$AutoMethod,')'), value = 0, {
      # waiter::waiter_show(html = waiter::spin_ring(), color = paste0("rgba(128,128,128,", 0.15, ")"))

      # Register keypoints
      mapping_list <- list()
      dest_image_list <- list()
      overlayed_image_list <- list()
      aligned_image_list <- list()
      alignment_image_list <- list()
      for(i in register_ind){

        # Increment the progress bar, and update the detail text.
        incProgress(1/length(register_ind), detail = paste("Registering Image", i, sep = " "))

        # get a sequential mapping between a query and reference image
        results <- computeAutomatedPairwiseTransform(image_list, channel_names, query_ind = i, ref_ind = centre, input)

        # save transformation matrix
        registration_mapping_list[[paste0(i)]] <- results$mapping

        # destination image
        dest_image_list[[i]] <- results$dest_image

        # save aligned images
        aligned_image_list[[i]] <- results$aligned_image

        # save alignment
        overlayed_image_list[[i]] <- results$overlay_image

        # save matches
        alignment_image_list[[i]] <- results$alignment_image
      }

      # waiter end
      # waiter::waiter_hide()
      })

      # Plot registered images
      lapply(register_ind, function(i){
        output[[paste0("plot_query_reg",i)]] <- renderImage({

          # get images
          image_view_list <- list(rep(magick::image_resize(dest_image_list[[i]], geometry = "400x"),5),
                                  rep(magick::image_resize(overlayed_image_list[[i]], geometry = "400x"),5))

          # make slide show
          image_view_list <- image_view_list %>%
            magick::image_join() %>%
            magick::image_write(tempfile(fileext='gif'), format = 'gif')
          list(src = image_view_list, contentType = "image/gif")
        }, deleteFile = TRUE)
      })

      # Plot Alignment
      lapply(register_ind, function(i){
        cur_alignment_image <- alignment_image_list[[i]]
        output[[paste0("plot_alignment",i)]] <- renderPlot({
          magick::image_ggplot(cur_alignment_image)
        })
      })

      # Output summary
      output[["summary"]] <- renderUI({
        str1 <- paste0(" Registration Summary:")
        str2 <- paste0("# of Images: ", length(image_list))
        str3 <- paste0("# of Registrations: ", len_register)
        all_str <- c(str1, str2, str3)
        shiny::HTML(paste(all_str, collapse = '<br/>'))
      })
    }
  })
}

#' computeAutomatedPairwiseTransform
#'
#' Computing the registration matrix necessary for automated registration
#'
#' @param image_list the list of images
#' @param channel_names the list of channel names for each image
#' @param query_ind the index of the query image
#' @param ref_ind the index of the reference image
#' @param input input
#'
#' @noRd
computeAutomatedPairwiseTransform <- function(image_list, channel_names, query_ind, ref_ind, input){

  # determine the number of transformation to map from query to the reference
  indices <- query_ind:ref_ind
  mapping_mat <- rep(indices,c(1,rep(2,length(indices)-2),1))
  mapping_mat <- matrix(mapping_mat,ncol=2,byrow=TRUE)

  # reference and target landmarks/keypoints
  mapping <- list()
  query_image <- image_list[[query_ind]]
  for(kk in 1:nrow(mapping_mat)){
    cur_map <- mapping_mat[kk,]
    ref_image <- image_list[[cur_map[2]]]

    # compute and get transformation matrix
    if(which.max(cur_map) == 1){
      ref_label = "ref"
      query_label = "query"
    } else {
      ref_label = "query"
      query_label = "ref"
    }

    # get channels 
    query_image <- query_image[[input[[paste0("channel_", query_label, "_image", cur_map[1])]]]]
    ref_image <- ref_image[[input[[paste0("channel_", ref_label, "_image", cur_map[2])]]]]
    
    # scale parameters
    query_scale <- input[[paste0("scale_", query_label, "_image", cur_map[1])]]
    ref_scale <- input[[paste0("scale_", ref_label, "_image", cur_map[2])]]

    # scale images
    # query_image <- magick::image_resize(query_image, geometry = magick::geometry_size_percent(100*query_scale))
    # ref_image <- magick::image_resize(ref_image, geometry = magick::geometry_size_percent(100*ref_scale))
    query_image <- resize_Image(query_image, geometry = magick::geometry_size_percent(100*query_scale))
    ref_image <- resize_Image(ref_image, geometry = magick::geometry_size_percent(100*ref_scale))

    # register images with OpenCV
    reg <- getRcppAutomatedRegistration(ref_image = ref_image, query_image = query_image,
                                        GOOD_MATCH_PERCENT = as.numeric(input$GOOD_MATCH_PERCENT), MAX_FEATURES = as.numeric(input$MAX_FEATURES),
                                        invert_query = input[[paste0("negate_", query_label, "_image", cur_map[1])]] == "Yes",
                                        invert_ref = input[[paste0("negate_", ref_label, "_image", cur_map[2])]] == "Yes",
                                        flipflop_query = input[[paste0("flipflop_", query_label, "_image", cur_map[1])]],
                                        flipflop_ref = input[[paste0("flipflop_", ref_label, "_image", cur_map[2])]],
                                        rotate_query = input[[paste0("rotate_", query_label, "_image", cur_map[1])]],
                                        rotate_ref = input[[paste0("rotate_", ref_label, "_image", cur_map[2])]],
                                        matcher = input$Matcher, method = input$Method)

    # update transformation matrix
    reg[[1]][[1]] <- solve(diag(c(ref_scale,ref_scale,1))) %*% reg[[1]][[1]] %*% diag(c(query_scale,query_scale,1))

    # return transformation matrix and images
    mapping[[kk]] <- reg[[1]]
    dest_image <- reg$dest_image
    aligned_image <- reg$aligned_image
    alignment_image <- reg$alignment_image
    overlay_image <- reg$overlay_image
  }

  return(list(mapping = mapping, 
              dest_image = dest_image, 
              aligned_image = aligned_image, 
              alignment_image = alignment_image, 
              overlay_image = overlay_image))
}

#' getRcppAutomatedRegistration
#'
#' Automated registration workflos with Rcpp
#'
#' @param ref_image reference image
#' @param query_image query image
#' @param GOOD_MATCH_PERCENT the percentage of good matching keypoints, used by "Brute force" method
#' @param MAX_FEATURES maximum number of detected features, i.e. keypoints, used by "Brute force" method
#' @param invert_query invert query image?
#' @param invert_ref invert reference image
#' @param flipflop_query flip or flop the query image
#' @param flipflop_ref flip or flop the reference image
#' @param rotate_query rotation of query image
#' @param rotate_ref rotation of reference image
#' @param matcher the matching method for landmarks/keypoints FLANN or BRUTE-FORCE
#' @param method the automated registration method, Homography or Homography+TPS
#'
#' @importFrom magick image_read image_data
#'
#' @export
getRcppAutomatedRegistration <- function(ref_image, query_image,
                                         GOOD_MATCH_PERCENT = 0.15, MAX_FEATURES = 500,
                                         invert_query = FALSE, invert_ref = FALSE,
                                         flipflop_query = "None", flipflop_ref = "None",
                                         rotate_query = "0", rotate_ref = "0", 
                                         matcher = "FLANN", method = "Homography") {
  ref_image_rast <- magick::image_data(ref_image, channels = "rgb")
  query_image_rast <- magick::image_data(query_image, channels = "rgb")
  reg <- automated_registeration_rawvector(ref_image = ref_image_rast, query_image = query_image_rast,
                                           width1 = dim(ref_image_rast)[2], height1 = dim(ref_image_rast)[3],
                                           width2 = dim(query_image_rast)[2], height2 = dim(query_image_rast)[3],
                                           GOOD_MATCH_PERCENT = GOOD_MATCH_PERCENT, MAX_FEATURES = MAX_FEATURES,
                                           invert_query = invert_query, invert_ref = invert_ref,
                                           flipflop_query = flipflop_query, flipflop_ref = flipflop_ref,
                                           rotate_query = rotate_query, rotate_ref = rotate_ref,
                                           matcher = matcher, method = method)
  
  # check for null keypoints
  if(suppressWarnings(all(lapply(reg[[1]][[2]], is.null)))){
    reg[[1]] <- list(reg[[1]][[1]], NULL)
  }
  
  return(list(transmat = reg[[1]],
              dest_image = magick::image_read(reg[[2]]),
              aligned_image = magick::image_read(reg[[3]]),
              alignment_image = magick::image_read(reg[[4]]),
              overlay_image = magick::image_read(reg[[5]])))
}

library(shiny)

ui <- fluidPage(
  tags$style("#myplot {cursor: crosshair;}"), # Changes cursor to crosshair over the plot
  plotOutput("myplot", 
             hover = hoverOpts(id = "plot_hover", delay = 100, delayType = "throttle")),
  tags$script(HTML("
    // Capture the scroll event and cursor position
    document.getElementById('myplot').addEventListener('wheel', function(e) {
      var bounds = this.getBoundingClientRect();
      var x = e.clientX - bounds.left;
      var y = e.clientY - bounds.top;
      
      // Determine whether it's zooming in or out
      var zoomDirection = (e.deltaY < 0) ? 'in' : 'out';
      
      // Send the zoom direction and cursor coordinates to Shiny
      Shiny.setInputValue('zoom_info', {zoom: zoomDirection, x: x, y: y}, {priority: 'event'});
      
      e.preventDefault(); // Prevent the page from scrolling
    });
  "))
)

server <- function(input, output, session) {
  
  # Initial plot ranges
  x_range <- reactiveVal(range(cars$speed))
  y_range <- reactiveVal(range(cars$dist))
  
  # Observe zoom information and adjust plot limits
  observeEvent(input$zoom_info, {
    zoom_factor <- ifelse(input$zoom_info$zoom == 'in', 1.5, 0.5)
    
    # Convert pixel coordinates to plot coordinates
    hover_info <- input$plot_hover
    if (!is.null(hover_info)) {
      cursor_x <- hover_info$x
      cursor_y <- hover_info$y
      
      # Get current ranges
      xlim <- x_range()
      ylim <- y_range()
      
      # Compute new ranges centered around the cursor position
      x_diff <- (xlim[2] - xlim[1]) * (1 - zoom_factor) / 2
      y_diff <- (ylim[2] - ylim[1]) * (1 - zoom_factor) / 2
      
      new_xlim <- c(cursor_x - (cursor_x - xlim[1]) * zoom_factor, cursor_x + (xlim[2] - cursor_x) * zoom_factor)
      new_ylim <- c(cursor_y - (cursor_y - ylim[1]) * zoom_factor, cursor_y + (ylim[2] - cursor_y) * zoom_factor)
      
      print(new_xlim)
      print(new_ylim)
      
      x_range(new_xlim)
      y_range(new_ylim)
    }
  })
  
  # Render the plot
  output$myplot <- renderPlot({
    plot(cars, xlim = x_range(), ylim = y_range(), main = "Scroll to Zoom (Cursor-Centered)")
  })
}

shinyApp(ui, server)
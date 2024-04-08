# Load required libraries
library(shiny)
library(ggplot2)

# UI
ui <- fluidPage(
  titlePanel("Gene Analysis Dashboard"),
  
  # Summary of analyses (mandatory)
  mainPanel(
    h3("Summary of Analyses"),
    textOutput("summary")
  ),
  
  # Gene query interface (optional)
  sidebarPanel(
    h3("Gene Query"),
    textInput("gene_query", label = "Enter Gene Name"),
    verbatimTextOutput("gene_info"),
    
    # Select plot type (optional)
    selectInput("plot_type", label = "Select Plot Type",
                choices = c("Scatter Plot", "Bar Plot"))
  ),
  
  # Plot output
  plotOutput("plot")
)

# Server
server <- function(input, output) {
  
  # Load data
  data <- read.csv("gene_data.csv")
  
  # Summary of analyses (dummy output for demonstration)
  output$summary <- renderText({
    "This is a summary of analyses."
  })
  
  # Gene query functionality (optional)
  output$gene_info <- renderPrint({
    gene <- input$gene_query
    if (!is.null(gene)) {
      gene_info <- data[data$Gene == gene, ]
      if (nrow(gene_info) > 0) {
        gene_info
      } else {
        "Gene not found"
      }
    }
  })
  
  # Render different types of plots based on user input
  output$plot <- renderPlot({
    # Check the selected plot type
    if (input$plot_type == "Scatter Plot") {
      # Scatter plot
      ggplot(data, aes(x = X, y = Y, color = Cell_Type)) +
        geom_point() +
        labs(title = "Scatter Plot of Cell Type Clusters")
    } else if (input$plot_type == "Bar Plot") {
      # Bar plot
      ggplot(data, aes(x = Gene, y = Expression, fill = Cell_Type)) +
        geom_bar(stat = "identity") +
        labs(title = "Bar Plot of Gene Expression by Cell Type")
    }
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
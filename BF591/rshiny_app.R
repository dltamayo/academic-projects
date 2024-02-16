## Dylan Tamayo
## BU BF591
## Final

library(shiny)
library(DT)
library(bslib)
library(colourpicker)
library(tidyverse)
library(bigPint)
library(igraph)
library(matrixStats, include.only = 'colVars')
library(RColorBrewer)
library(glue)

# UI
ui <- fluidPage(
  tags$head(
    # Note the wrapping of the string in HTML()
    tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Holtwood+One+SC&display=swap');
      body {
        background-color: #99F2A5;
        color: #0c0f0a;
      }
      h1 {
        font-family: 'Holtwood One SC';
      }
      .shiny-input-container {
        color: #89fc00;
      }
      h2 {
        font-family: 'Holtwood One SC';
        color: #3333FF;
      }
      .shiny-input-container {
        color: #474747;
      }
      
      h3 {
        font-family: 'Holtwood One SC';
      }
      .shiny-input-container {
        color: #0000FF;
      }"))
  ),
  #color: #0000FF; fbff12; 0c0f0a ff206e
    
  #titlePanel(
    #h2('Test title panel'),
  tabsetPanel(
    tabPanel(h1('Samples'), h1(''),
      sidebarLayout(
        sidebarPanel(h3('User inputs'),
          fileInput('sample_input', 'Input CSV file'),
          
          submitButton('Update settings', icon('refresh'))),
        mainPanel(
          tabsetPanel(
            tabPanel(h3('Summary'),
              DTOutput('sample_summary')),
            tabPanel(h3('Table'),
              DTOutput('sample_table')),
            tabPanel(h3('Plots'),
                     sidebarLayout(
                       sidebarPanel(h3('User inputs'),
                                    selectInput('sample_x', 'X-axis:',
                                                choices = c('Developmental stage', 'Age at collection (days after planting)',
                                                            'Lighting (μmol m-2 s-1)', 'Temperature (C)', 'Relative Humidity (%)'),
                                                selected = 'Age at collection (days after planting)'),
                                    selectInput('sample_y', 'Y-axis:',
                                                choices = c('Developmental stage', 'Age at collection (days after planting)',
                                                            'Lighting (μmol m-2 s-1)', 'Temperature (C)', 'Relative Humidity (%)'),
                                                selected = 'Temperature (C)'),
                                    selectInput('sample_group', 'Group by:',
                                                choices = c('Developmental stage', 'Age at collection (days after planting)',
                                                            'Lighting (μmol m-2 s-1)', 'Temperature (C)', 'Relative Humidity (%)'),
                                                selected = 'Developmental stage'),
                                    selectInput('sample_plot', 'Plot type:', choices = c('Histogram', 'Line', 'Point')),
                                    submitButton('Update settings', icon('refresh'))),
                       mainPanel(
                         plotOutput('sample_plot', width = "500px", height = "500px"))))
          )
        )
      )
    ),
    
    tabPanel(h1('Counts'), h1(''),
      sidebarLayout(
        sidebarPanel(h3('User inputs'),
          fileInput('counts_input', 'Input CSV file'),
          sliderInput('counts_variance', 'Include genes above X percentile variance:', 0, 100, 80),
          sliderInput('counts_expression', 'Include genes with at least X non-zero samples:', 0, 9, 4),
          submitButton('Update settings', icon('refresh'))),
        mainPanel(
          tabsetPanel(
            tabPanel(h3('Filtered Stats'),
              tableOutput('counts_table')),
            tabPanel(h3('Scatter Plot'),
              sidebarLayout(
                sidebarPanel(h3('User inputs'),
                  radioButtons('counts_scatter_x', 'Scale x-axis:', c(TRUE, FALSE), inline = TRUE),
                  radioButtons('counts_scatter_y', 'Scale y-axis:', c(TRUE, FALSE), inline = TRUE),
                  submitButton('Update settings', icon('refresh'))),
                mainPanel(
                  plotOutput('counts_scatter', width = "500px", height = "500px")))),
            tabPanel(h3('Heatmap'),
              plotOutput('counts_heatmap', width = "500px", height = "500px")),
            tabPanel(h3('PCA'),
              sidebarLayout(
                sidebarPanel(h3('User inputs'),
                  sliderInput('counts_pca_x', 'PC for x-axis:', 1, 9, 1),
                  sliderInput('counts_pca_y', 'PC for y-axis:', 1, 9, 2),
                  submitButton('Update settings', icon('refresh'))),
                mainPanel(
                  plotOutput('counts_pca', width = "500px", height = "500px")))),
          )       
        )
      )
    ),
    
    tabPanel(h1('Diff Expr'), h1(''),
      sidebarLayout(
        sidebarPanel(h3('User inputs'),
          fileInput('de_input', 'Input CSV file'),
          sliderInput('de_slider', 'Magnitude for false discovery rate cutoff:',
                      0, 1, 0.05, step = 0.005),
          submitButton('Update settings', icon('refresh'))),
        mainPanel(
          tabsetPanel(
            tabPanel(h3('DE Results'),
              DTOutput('de_table')),
            tabPanel(h3('DE Plot'),
              sidebarLayout(
                sidebarPanel(h3('User inputs'),
                  selectInput('de_plot_x', 'Column for x-axis:',
                              choices = c('logFC', 'logCPM', 'LR', 'PValue', 'FDR'),
                              selected = 'logFC'),
                  selectInput('de_plot_y', 'Column for y-axis:',
                              choices = c('logFC', 'logCPM', 'LR', 'PValue', 'FDR'),
                              selected = 'FDR'),
                  colourInput('de_color1', 'Color above FDR cutoff:', '#780000'),
                  colourInput('de_color2', 'Color below FDR cutoff:', '#0077B6'),
                  submitButton('Update settings', icon('refresh'))),
                mainPanel(
                  plotOutput('de_plot', width = "500px", height = "500px"))))
          )       
        )
      )
    ),
  
    
    tabPanel(h1('Network'), h1(''),
      sidebarLayout(
        sidebarPanel(h3('User inputs'),
          fileInput('network_input', 'Input CSV file'),
          textAreaInput('network_subset', 'Input genes to subset:',
                        height = '215px', value = 'Glyma18g00690.1
Glyma03g29150.1
Glyma02g40610.1
Glyma10g31780.1
Glyma13g26600.1
Glyma16g08810.1
Glyma08g22380.1
Glyma14g05270.1
Glyma01g42820.2
Glyma08g19245.1'),
          sliderInput('network_slider', 'Minimum correlation for drawing edge:',
                      0, 1, 0.995, step = 0.005),
          submitButton('Update settings', icon('refresh'))),
        mainPanel(
          tabsetPanel(
            tabPanel(h3('Heatmap'),
              plotOutput('network_heatmap', width = "600px", height = "500px")),
            tabPanel(h3('Network Plot'),
              sidebarLayout(
                sidebarPanel(h3('User inputs'),
                  colourInput('network_node_color', 'Node color:', '#00BB22'),
                  colourInput('network_edge_color', 'Edge color:', '#7700DD'),
                  colourInput('network_label_color', 'Label color:', '#000000'),
                  colourInput('network_background_color', 'Background color:', '#FFFFFF'),
                  submitButton('Update settings', icon('refresh'))),
                mainPanel(
                  plotOutput('network_plot', width = "500px", height = "500px")))),
            tabPanel(h3('Metrics'),
              tableOutput('network_table')),
          )       
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # Network functions
  load_network_data <- reactive({
    read_csv(input$network_input$datapath)
    # read_csv('/Users/dlt/Downloads/soybean_counts.csv')
  })
  
  parse_subset <- function(raw_gene_list) {
    str_split(raw_gene_list, '\n')[[1]]
  }
  
  network_subset <- function(counts_data, gene_list) {
    counts_data %>%
      dplyr::filter(gene %in% gene_list) %>%
      dplyr::select(c(1,2,3,4,8,9,10)) %>%
      column_to_rownames('gene') %>%
      as.matrix()
  }
  
  network_heatmap <- function(gene_subset) {
    mypalette <- brewer.pal(6, 'PRGn')
    par(mar=c(5, 4, 4, 8), xpd=TRUE)
    heatmap(gene_subset, col = mypalette)
    legend(x='topleft', inset=c(-0.05, -0.125), legend = c("low", NA, NA, NA, NA, "high"), fill=mypalette)
  }
  
  create_network <- function(gene_subset, min_corr) {
    mat <- cor(t(gene_subset))
    mat[mat < min_corr] <- 0
    network <- graph_from_adjacency_matrix(mat, weighted = T, mode = "undirected", diag = F)
  }
  
  plot_network <- function(network, node_color, edge_color, label_color, background_color) {
    par(bg=glue('{background_color}'), mar=c(1.5,1.5,0,0))
    set.seed(1)
    plot(network,
         # === vertex
         vertex.color = alpha(glue('{node_color}'), alpha = 0.5),          # Node color
         vertex.frame.color = alpha(glue('{node_color}'), alpha = 0.75),    # Node border color
         vertex.shape = "circle",                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
         vertex.size = 14,                               # Size of the node (default is 15)
         
         # === vertex label
         vertex.label.color = glue('{label_color}'),
         vertex.label.family = "Helvetica",              # Font family of the label (e.g.“Times”, “Helvetica”)
         vertex.label.font = 1,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
         vertex.label.cex = 1,                           # Font size (multiplication factor, device-dependent)
         vertex.label.dist = 1.5,                        # Distance between the label and the vertex
         vertex.label.degree = pi/2,                     # The position of the label in relation to the vertex (use pi)
         
         # === edge
         edge.color = alpha(glue('{edge_color}'), alpha = 0.5), # Edge color
         edge.width = 3,                                 # Edge width, defaults to 1
         edge.arrow.size = 1,                            # Arrow size, defaults to 1
         edge.arrow.width = 1,                           # Arrow width, defaults to 1
         edge.lty = "solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
         edge.curved = 0.1                               # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
    )
  }
  
  network_metrics <- function(gene_subset, network) {
    tibble(Gene = rownames(gene_subset),
           Degree = degree(network),
           'Closeness centrality' = closeness(network),
           'Betweenness centrality' = betweenness(network))
  }
  
  # Network output
  output$network_heatmap <- renderPlot({
    gene_subset <- network_subset(load_network_data(), parse_subset(input$network_subset))
    network_heatmap(gene_subset)
  })
  
  output$network_plot <- renderPlot({
    gene_subset <- network_subset(load_network_data(), parse_subset(input$network_subset))
    network <- create_network(gene_subset, input$network_slider)
    plot_network(network, input$network_node_color, input$network_edge_color,
                 input$network_label_color, input$network_background_color)
  })
  
  output$network_table <- renderTable({
    gene_subset <- network_subset(load_network_data(), parse_subset(input$network_subset))
    network <- create_network(gene_subset, input$network_slider)
    network_metrics(gene_subset, network)
  })
  
  # Sample functions
  load_sample_data <- reactive({
    read_csv(input$sample_input$datapath)
    # read_csv('/Users/dlt/Downloads/soybean_metadata.csv')
  })
  
  summarize_metadata <- function(meta_data) {
    non_sample <- meta_data[-1]
    non_sample %>%
      reframe(column = colnames(non_sample),
              type = sapply(non_sample, class),
              value = unname(lapply(non_sample, base::unique)) #apply(non_sample, 2, unique)
      )
  }
  
  plot_metadata <- function(meta_data, column_x, column_y, column_group, plot_type) {
    x <- rlang::sym(column_x)
    y <- rlang::sym(column_y)
    group <- rlang::sym(column_group)
    
    meta_data %>%
      ggplot(aes(x = {{x}})) +
      theme_bw() +
      if (plot_type == 'Histogram') {
        geom_histogram(aes(fill = {{group}}), binwidth = 1)
      } else if (plot_type == 'Line') {
        geom_line(aes(y = {{y}}))
      } else if (plot_type == 'Point') {
        geom_point(aes(y = {{y}}, color = {{group}}))
      }
  }
  
  # Sample outputs
  output$sample_summary <- renderDataTable({
    summarize_metadata(load_sample_data())
  })
  output$sample_table <- renderDataTable({
    load_sample_data()
  })
  output$sample_plot <- renderPlot({
    plot_metadata(load_sample_data(), input$sample_x, input$sample_y,
                 input$sample_group, input$sample_plot)
  })
  
  # Counts functions
  load_counts_data <- reactive({
    read_csv(input$counts_input$datapath)
    # read_csv('/Users/dlt/Downloads/soybean_counts.csv')
  })
  
  filter_variance <- function(count_data, variance_cutoff) {
    variance_percent <- variance_cutoff/100
    filtered_gene_list <- count_data %>%
      tidyr::pivot_longer(-c(gene), names_to='sample') %>%
      dplyr::group_by(gene) %>%
      dplyr::summarize(var=var(value))
    
    quant <- quantile(filtered_gene_list$var, probs = variance_percent)
    
    filtered_gene_list <- filtered_gene_list %>%
      dplyr::filter(var > quant) %>%
      dplyr::pull(gene)
    
    count_data %>%
      dplyr::filter(!gene %in% filtered_gene_list) %>%
      dplyr::pull(gene)
  }
  
  filter_zero_expr <- function(count_data, zero_sample_cutoff) {
    count_data %>%
      dplyr::filter(!(rowSums(across(-gene) != 0) >= zero_sample_cutoff)) %>%
      dplyr::pull(gene)
  }
  
  double_filter <- function(count_data, variance_cutoff, zero_sample_cutoff) {
    filter_var <- count_data %>% filter_variance(variance_cutoff)
    filter_expr <- count_data %>% filter_zero_expr(zero_sample_cutoff)
    
    filter_levels <- c('Filters passed', 'Expression below threshold', 'Variance below threshold', 'Variance and expression below threshold')
    master_filter <- count_data %>%
      dplyr::mutate(filtered = case_when((gene %in% filter_var) & (gene %in% filter_expr) ~ 'Variance and expression below threshold',
                                         gene %in% filter_var ~ 'Variance below threshold',
                                         gene %in% filter_expr ~ 'Expression below threshold',
                                         (!gene %in% filter_var) & (!gene %in% filter_expr) ~ 'Filters passed'),
                    filtered = factor(filtered, levels = filter_levels)) 
  }
  
  summary_table <- function(count_data, filtered_counts) {
    count_dim <- dim(count_data %>% dplyr::select(where(is.numeric)))
    filter_dim <- dim(filtered_counts)
    
    num_sample <- count_dim[2]
    total_genes <- count_dim[1]
    num_pass <- filter_dim[1]
    percent_pass <- round(num_pass/total_genes*100, 2)
    num_fail <- total_genes - num_pass
    percent_fail <- round(num_fail/total_genes*100, 2)
    
    output <- tibble(num_sample, total_genes, num_pass, percent_pass, num_fail, percent_fail)
    colnames(output) <- c('Number of samples', 'Total genes', 'Number passed', 'Percent passed', 'Number filtered', 'Percent filtered')
    output
  }
  
  plot_variance_vs_median <- function(data, scale_x_axis=FALSE, scale_y_axis=FALSE, title="") {
    new_data <- data %>% dplyr::select(-filtered) %>% 
      tidyr::pivot_longer(-c(gene),names_to="sample") %>%
      tidyr::pivot_wider(names_from=c(gene))
    
    new_data %>%
      dplyr::reframe(gene = data$gene,
                     mean = colMeans(new_data[-1]),
                     median = apply(new_data[-1], 2, median),
                     var = colVars(as.matrix(new_data[-1])),
                     filtered = data$filtered) %>%
      dplyr::mutate(rank = rank(median, ties.method='average')) %>%
      
      ggplot(aes(x = median, y = var, color = filtered)) +
      theme_bw() +
      geom_point(alpha=0.3) +
      labs(title=title,
           x = 'Gene Median Expression',
           y = 'Gene Variance') +
      theme(legend.position = 'bottom') +
      (if(scale_x_axis){scale_x_log10()} else{scale_x_continuous()}) + 
      (if(scale_y_axis){scale_y_log10()} else{scale_y_continuous()})
  }
  
  plot_heatmap <- function(filtered_counts) {
    filtered_matrix <- filtered_counts %>%
      dplyr::select(!filtered) %>%
      tibble::column_to_rownames('gene') %>%
      as.matrix()
    #log_passed <- log10(passed_filter)
    mypalette <- brewer.pal(6, 'PRGn')
    
    par(mar=c(5, 4, 4, 8), xpd=TRUE)
    heatmap(filtered_matrix, col = mypalette) #log_passed
    legend(x='topleft', inset=c(-0.15, -0.125), legend = c("low", NA, NA, NA, NA, "high"), fill=mypalette)
  }
  
  conduct_pca <- function(filtered_counts) {
    counts_matrix <- filtered_counts %>%
      dplyr::select(where(is.numeric)) %>%
      as.matrix()
    
    pca_results <- prcomp(scale(t(counts_matrix)), center=FALSE, scale=FALSE)
  }
  
  calculate_variance_explained <- function(pca_results) {
    sdev <- pca_results$sdev**2
    sdev/sum(sdev)
  }
  
  variance_tibble <- function(pca_results) {
    var <- calculate_variance_explained(pca_results)
    
    pca_var <- tibble(
      PC=factor(str_c("PC",1:length(var)),
                levels=str_c("PC",1:length(var))
      ),
      variance = var,
      cumulative = cumsum(variance))
  }
  
  plot_pca <- function(pca_results, variance_tibble, x, y) {
    PC_x <- rlang::sym(glue('PC{x}'))
    PC_y <- rlang::sym(glue('PC{y}'))
    pca_df <- as.data.frame(pca_results$x[,c(x,y)])
    var_x <- round(variance_tibble$variance[x] * 100,1)
    var_y <- round(variance_tibble$variance[y] * 100,1)
    
    graph_tibble <- tibble(rownames_to_column(pca_df, 'sample'))
    graph_tibble %>%
      ggplot(aes(x={{PC_x}}, y={{PC_y}}, color=sample)) +
      theme_bw() +
      geom_point() +
      labs(x = glue('PC{x} ({var_x}% Variance)'),
           y = glue('PC{y} ({var_y}% Variance)'))
  }
  
  # Counts outputs
  output$counts_table <- renderTable({
    filtered_data <- double_filter(load_counts_data(), input$counts_variance, input$counts_expression) %>% 
      dplyr::filter(filtered == 'Filters passed')
    summary_table(load_counts_data(), filtered_data)
  })
  
  output$counts_scatter <- renderPlot({
    filtered_data <- double_filter(load_counts_data(), input$counts_variance, input$counts_expression)
    plot_variance_vs_median(filtered_data, input$counts_scatter_x, input$counts_scatter_y, '')
    
  })
  
  output$counts_heatmap <- renderPlot({
    filtered_data <- double_filter(load_counts_data(), input$counts_variance, input$counts_expression) %>% 
      dplyr::filter(filtered == 'Filters passed')
    plot_heatmap(filtered_data)
  })
  
  output$counts_pca <- renderPlot({
    filtered_data <- double_filter(load_counts_data(), input$counts_variance, input$counts_expression) %>% 
      dplyr::filter(filtered == 'Filters passed')
    pca_res <- conduct_pca(filtered_data)
    pca_var <- variance_tibble(pca_res)
    plot_pca(pca_res, pca_var, input$counts_pca_x, input$counts_pca_y)
  })
  
  # Differential Expression functions
  load_de_data <- reactive({
    read_csv(input$de_input$datapath)
    # read_csv('/Users/dlt/Downloads/soybean_1v3_deseq.csv')
  })
  
  filter_de_table <- function(dataf, slider) {
    output <- dataf %>%
      dplyr::filter(FDR < slider) %>%
      mutate(logFC = round(logFC, digits = 3),
             logCPM = round(logCPM, digits = 3),
             LR = round(LR, digits = 3),
             PValue = format(PValue, digits = 3, scientific = FALSE),
             FDR = format(FDR, digits = 3, scientific = FALSE))
  }
  
  plot_de <- function(dataf, x_name, y_name, slider, color1, color2) {
    x_column <- rlang::sym(x_name)
    plot_tib <- dataf %>%
      dplyr::mutate('-log10(FDR)' = -log10(FDR),
                    '-log10(PValue)' = -log10(PValue),
                    pval_cutoff = FDR < slider
      )
    if(y_name == 'FDR' & x_name == 'logFC') {y_name <- '-log10(FDR)'}
    else if(y_name == 'PValue' & x_name == 'logFC') {y_name <- '-log10(PValue)'}
    y_column <- rlang::sym(y_name)
    
    plot_tib %>%
      #drop_na(padj) %>%
      ggplot(aes(x = {{x_column}},y = {{y_column}})) +
      geom_point(aes(color = pval_cutoff)) +
      labs(x = x_name,
           y = y_name,
           title = '') +
      theme_bw() +
      scale_color_manual(name = 'FDR < cutoff', values=c(color1, color2)) +
      theme(legend.position = "bottom")
  }
  
  # Differential Expression outputs
  output$de_table <- renderDataTable({
    filter_de_table(load_de_data(), input$de_slider)
  })
  
  output$de_plot <- renderPlot({
    plot_de(load_de_data(), input$de_plot_x, input$de_plot_y,
            input$de_slider, input$de_color1, input$de_color2)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
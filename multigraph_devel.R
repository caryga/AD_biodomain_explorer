#
# This application plots information about scored targets from the TREAT-AD Emory-Sage-SGC center, including:
#   - Overall, Genetics, Omics, Literature, Neuropathology, and Druggability scores
#   - Meta-analysis predictions of RNA and Protein expression changes
#   - AD Biodomain annotations and Biodomain term network
#   - Cell-type specific gene expression from SEA-AD
#   - Phenotype data from the CRISPRbrain resource
#

# Load Packages -----------------------------------------------------------

require(network)
require(sna)
require(intergraph)
require(fst)
library(igraph)
library(GGally)
library(shiny)
library(tidyverse)

theme_set(theme_bw(base_size = 14))

# profvis::profvis({

# Load Pre-formatted Data -------------------------------------------------

##
# load from local data/ dir
##

biodom_genes = fst::read.fst('data/biodom_genes.fst') %>% as_tibble()

dom.cols = readRDS('data/dom_cols.rds')
term.graph = readRDS('data/term_graph.rds')
# trs.graph <- readRDS('data/trs_graph.rds')
# gen.graph <- readRDS('data/gen_graph.rds')
# txome.graph <- readRDS('data/txome_graph.rds')
# pxome.graph <- readRDS('data/pxome_graph.rds')

# Shiny App ---------------------------------------------------------------------

# [enter] key as button click ----
jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'


# User Interface ----
ui <- fluidPage(
  
  # Application title
  titlePanel("AD Biological Domain Explorer"),
  
  # Input row
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      div(
        tags$head(tags$script(HTML(jscode))),
        
        # TODO: filter term choices by biodomain
        # TODO: add node statistic to visualization
        
        radioButtons(
          inputId = 'term_nw_select',
          label = 'Select biodomain term network',
          choiceNames = c('term annotation',
                          'AD risk enriched',
                          'AD genetic risk enriched',
                          'AD transcriptome',
                          'AD proteome'),
          choiceValues = c('term.graph',
                           'trs.graph',
                           'gen.graph', 
                           'txome.graph', 
                           'pxome.graph'), 
          selected = 'term.graph'
        ),
        
        hr(),
        tagAppendAttributes(
          textInput(inputId = "term_accession", 
                    label = "Term Accession:",
                    placeholder = 'e.g. GO:0005743'),
          `data-proxy-click` = "update"),
        
        strong('- or -'),
        br(),
        
        tagAppendAttributes(
          selectizeInput(
            inputId = "term_name", 
            label = "Term Name:",
            choices = NULL, 
            selected = NULL, 
            multiple = FALSE,
            options = list(maxOptions = 5,
                           placeholder = 'e.g. mitochondrial inner membrane')
          ),
          `data-proxy-click` = "update"),
        
        hr(),
        selectInput(inputId = "edge_depth", 
                    label = "Edge Depth:",
                    choices = c(1,2,3),
                    selected = 1),
        br(),
        sliderInput(inputId = "edge_kappa_coefficient", 
                    label = "Edge Kappa Coefficient:",
                    min = 0.1, max = 1, value = 0.8, step = 0.01),
        
        hr(),
        actionButton("update", "Plot", class = "btn-primary"),
      
        )
    ),
    mainPanel(
      plotly::plotlyOutput('biodom_nw', height = '800px'), #
      p(strong("Legend:"),
        "Interactive network of terms topographically near the query term based on shared gene annotation. Colored nodes represent biological domains and grey nodes are gene ontology terms to which the gene is annotated, where the size of the node represents the node degree in this sub-plot. Edges are drawn between terms that have genes in common (Cohen's kappa coefficient above the specified value)."
        , style = "margin-left: 30px;"
      )
    )
  ),
  
  # For more information links
  fluidRow(
    column(width = 12, 
           hr(),
           uiOutput('linkOut')
    )
  )
)



# Server ----
# Define server logic required to draw plots for query gene
server <- function(input, output, session) {
  
  # updateSelectizeInput(session, 'term_id', choices = unique(biodom_genes$GO_ID), server = TRUE)
  updateSelectizeInput(session, 'term_name', 
                       choices = unique(biodom_genes$GOterm_Name), 
                       selected = NULL,
                       server = TRUE)
  
  # plot initialization ----
  rv = reactiveValues()
  
  rv$generate_biodom_plot = NULL
  
  # biodomain plot ----    
  observe({
    
    # pull GO ID from query inputs
    if( !is.null(input$term_name) ){
      
      tg.terms <- biodom_genes %>% 
        filter(GOterm_Name == input$term_name) %>% 
        select(GO_ID, GOterm_Name, Biodomain) %>% 
        distinct()
      
      tg_id = unique(tg.terms$GO_ID)
      tg_name = unique(tg.terms$GOterm_Name)
      
    } else if( !is.null(input$term_accession) ){
      
      tg.terms <- biodom_genes %>% 
        filter(GO_ID == input$term_accession) %>% 
        select(GO_ID, GOterm_Name, Biodomain) %>% 
        distinct()
      
      tg_id = unique(tg.terms$GO_ID)
      tg_name = unique(tg.terms$GOterm_Name)
      
    } 
    
    graph = eval(parse(text = input$term_nw_select))
    depth = input$edge_depth
    kappa_thresh = input$edge_kappa_coefficient
    
    # remove edges below kappa score threshold
    sub.g = delete.edges(graph ,  
                         edges = E(graph)[ E(graph)$kappa < kappa_thresh ] )
    
    if( 'n_le_genes' %in% vertex_attr_names(sub.g) ){
      sub.g = delete.vertices(sub.g,
                              v = V(sub.g)[ !is.na(V(sub.g)$none) & V(sub.g)$none == 1]
                              )
    }
    
    # identify terms N hops away from the query
    terms <- neighborhood(sub.g, order = depth, nodes = tg_id)[[1]]
    
    # remove terms not annotated to term
    sub.g <- delete.vertices(
      sub.g, 
      V(sub.g)[ !( V(sub.g) %in% unique(terms)
                   # | !is.na(igraph::V(sub.g)$label) 
                   )
                ]
      )
    
    # Remove the unconnected nodes.
    unconnected_nodes <- which( degree(sub.g) == 0 )
    sub.g <- delete.vertices(sub.g, unconnected_nodes)
    
    # remove redundant edges 
    simpl <- igraph::simplify(sub.g, edge.attr.comb = c(kappa = 'max'))
    
    # add node degree to graph
    V(simpl)$degree <- degree(simpl)
    
    V(simpl)$shape <- 20
    V(simpl)$shape[ which( names(V(simpl)) == tg_id) ] <- 18
    # V(simpl)$color[ which( names(V(simpl)) == tg_id) ] <- 'red'
    
    # TODO: color term nodes based on risk enrichment?
    
    bd.nw <- suppressWarnings(
      ggnet2(
        net = simpl,
        mode = 'kamadakawai',
        edge.color = 'grey50',
        edge.lty = 3,
        # node.size = 'degree',
        max_size = 9,
        node.color = 'NES',
        node.shape = 'shape',
        # node.alpha = 'degree',
        label = 'abbr' ,
        label.size = 3
      ) +
        theme(legend.position = 'none',
              plot.title = element_text(size = 10)) +
        geom_point( aes(text = paste0(V(simpl)$term, '\n', V(simpl)$n_symbol, ' genes')),
                    color = 'grey80', alpha = 0, size = 8) + #
        ggtitle(paste0(tg_name, ' (', tg_id,') network neighborhood'))
    )
    
    rv$generate_bd_nw = bd.nw 
    
  }) %>%  bindEvent(input$update, ignoreInit = F)
  
  
  # plot rendering ---- 
  output$biodom_nw <- plotly::renderPlotly( plotly::ggplotly(rv$generate_bd_nw, tooltip = 'text')  ) # %>% layout(width = 350)

  output$linkOut <- renderUI({ 
    list(
      HTML('<strong>Link Outs: </strong>'),
      HTML(paste0('<a href = "http://amigo.geneontology.org/amigo/search/ontology?q=', input$symbol, '" target="_blank">Gene Ontology</a>'))
    )
  }) 
  
}

# Run the application 
shinyApp(ui = ui, server = server)

# })

### EOF ###

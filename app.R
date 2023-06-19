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

map_viridis <- function(vec, num) {
  
  vector_expanded <-round(vec, 1) * 10 # expand to allow for decimal precision
  vector_exp_range <- max(vector_expanded) - min(vector_expanded) 
  
  colour_vector <- viridis(vector_exp_range + 1, option = 'viridis') # get vector of colour values for all possible decimals between min and max value
  value_to_colour <- colour_vector[num * 10 - min(vector_expanded) + 1] # retrieve colour value for number
  
  return(value_to_colour)
  
}

# Load Pre-formatted Data -------------------------------------------------

##
# load from local data/ dir
##

biodom_genes = fst::read.fst('data/biodom_genes.fst') %>% as_tibble()
dom.cols = readRDS('data/dom_cols.rds')
term.graph = readRDS('data/term_graph.rds')

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
        
        # radioButtons(
        #   inputId = 'term_nw_select',
        #   label = 'Select biodomain term network',
        #   choiceNames = c('term annotation',
        #                   'AD risk enriched',
        #                   'AD genetic risk enriched',
        #                   'AD transcriptome',
        #                   'AD proteome'),
        #   choiceValues = c('term.graph',
        #                    'trs.graph',
        #                    'gen.graph', 
        #                    'txome.graph', 
        #                    'pxome.graph'), 
        #   selected = 'term.graph'
        # ),
        
        # tagAppendAttributes(
        #   textInput(inputId = "term_accession", 
        #             label = "Term Accession:",
        #             placeholder = 'e.g. GO:0005743'),
        #   `data-proxy-click` = "update"),
        # 
        # strong('- or -'),
        # br(),br(),
        
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
        
        radioButtons(inputId = "node_color", 
                    label = "Color Nodes By:",
                    choiceNames = c('Genetics NES', 'Genetics p.adjust',
                                     'Omics NES', 'Omics p.adjust',
                                     'proteomics NES', 'proteomics p.adjust',
                                     'transcriptomics NES', 'transcriptomics p.adjust',
                                     'TRS NES', 'TRS p.adjust'),
                    choiceValues = c('gen_NES', 'gen_padj',
                                'omic_NES', 'omic_padj',
                                'pro_NES', 'pro_padj',
                                'rna_NES', 'rna_padj',
                                'trs_NES', 'trs_padj'),
                    selected = 'trs_padj'),
        
        br(),
        
        selectInput(inputId = "edge_depth", 
                    label = "Edge Depth:",
                    choices = c(1,2,3),
                    selected = 1),
        br(),
        sliderInput(inputId = "edge_kappa_coefficient", 
                    label = "Edge Kappa Coefficient:",
                    min = 0.1, max = 1, value = 0.5, step = 0.1),
        
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
    
    # # pull GO ID from query inputs
    # if( !is.null(input$term_name) ){
      
    tg.terms <- biodom_genes %>% 
      filter(GOterm_Name == input$term_name) %>% 
      select(GO_ID, GOterm_Name, Biodomain) %>% 
      distinct()
    
    tg_id = unique(tg.terms$GO_ID)
    tg_name = unique(tg.terms$GOterm_Name)
      
    # } else if( !is.null(input$term_accession) ){
    #   
    #   tg.terms <- biodom_genes %>% 
    #     filter(GO_ID == input$term_accession) %>% 
    #     select(GO_ID, GOterm_Name, Biodomain) %>% 
    #     distinct()
    #   
    #   tg_id = unique(tg.terms$GO_ID)
    #   tg_name = unique(tg.terms$GOterm_Name)
    #   
    # } 
    
    # graph = eval(parse(text = input$term_nw_select))
    graph = term.graph
    depth = input$edge_depth
    kappa_thresh = input$edge_kappa_coefficient
    node_attr = input$node_color
    
    # remove edges below kappa score threshold
    sub.g = igraph::delete.edges(graph ,  
                         edges = E(graph)[ E(graph)$kappa < kappa_thresh ] )
    
    # if( 'n_le_genes' %in% vertex_attr_names(sub.g) ){
    #   sub.g = delete.vertices(sub.g,
    #                           v = V(sub.g)[ !is.na(V(sub.g)$none) & V(sub.g)$none == 1]
    #                           )
    # }
    
    # identify terms N hops away from the query
    terms <- igraph::neighborhood(sub.g, order = depth, nodes = tg_id)[[1]]
    
    # remove terms not annotated to term
    term.g <- igraph::delete.vertices(
      sub.g, 
      V(sub.g)[ !( V(sub.g) %in% unique(terms) ) & V(sub.g)$bd == F ]
      )
    
    # Remove the unconnected nodes.
    unconnected_nodes <- which( igraph::degree(term.g) == 0 )
    term.g <- igraph::delete.vertices(term.g, unconnected_nodes)
    
    # remove redundant edges 
    simpl <- igraph::simplify(term.g, edge.attr.comb = c(kappa = 'max'))
    
    # add node degree to graph
    V(simpl)$degree <- igraph::degree(simpl)
    V(simpl)$degree[ which( V(simpl)$bd ) ] <- max(igraph::degree(simpl), na.rm = T) *1.5
    
    V(simpl)$shape <- 20
    V(simpl)$shape[ which( names(V(simpl)) == tg_id) ] <- 18
    
    # gen_NES gen_padj
    # omic_NES omic_padj
    # pro_NES pro_padj
    # rna_NES rna_padj
    # trs_NES trs_padj
    
    V(simpl)$plot_color <- map_chr( vertex_attr(simpl, node_attr) , ~ map_viridis( vertex_attr(sub.g, node_attr), .x))
    V(simpl)$plot_color[ which( V(simpl)$bd ) ] <- V(simpl)$color[ which( V(simpl)$bd ) ]
    
    bd.nw <- suppressWarnings(
      ggnet2(
        net = simpl,
        mode = 'spring',
        edge.color = 'grey85',
        # edge.lty = 3,
        node.size = 'degree',
        # max_size = 9,
        node.color = 'plot_color',
        node.shape = 'shape',
        # node.alpha = 'degree',
        label = 'abbr' ,
        label.size = 3
      ) +
        theme(legend.position = 'none',
              plot.title = element_text(size = 10)) +
        geom_point( aes(text = paste0(V(simpl)$term, 
                                      '\n', V(simpl)$n_symbol, ' genes',
                                      '\npadj = ', signif( 
                                        vertex_attr(simpl,
                                                    node_attr %>% str_split_fixed(., '_', 2) %>% .[1] %>% paste0(., '_padj') )), 
                                      '\nNES = ', signif( 
                                        vertex_attr(simpl,
                                                    node_attr %>% str_split_fixed(., '_', 2) %>% .[1] %>% paste0(., '_NES') ))
                                        ) ),
                    color = 'grey80', alpha = 0, size = 12) + #
        ggtitle(paste0(tg_name, ' (', tg_id,') network neighborhood'))
    )
    
    rv$generate_bd_nw = bd.nw
    
  }) %>%  bindEvent(input$update, ignoreInit = T)
  
  
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

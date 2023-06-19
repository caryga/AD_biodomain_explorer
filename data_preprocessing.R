# Load Packages -----------------------------------------------------------

# require(network)
# require(sna)
# require(intergraph)
require(fst)
library(igraph)
library(shiny)
library(tidyverse)

theme_set(theme_bw(base_size = 14))

# Load and format data -------------------------------------------

# 

# ## Gene-Biodomain mappings ----
# synapser::synLogin()
# biodom_genes <- readRDS(synapser::synGet('syn25428992')$path) %>%
#   select(Biodomain, GO_ID, hgnc_symbol) %>%
#   unnest_longer(hgnc_symbol) %>%
#   filter(!is.na(hgnc_symbol), hgnc_symbol != '') %>%
#   group_by(hgnc_symbol, Biodomain) %>%
#   summarise(n_term = length(unique(GO_ID))) %>%
#   ungroup() %>%
#   left_join(
#     .,
#     readRDS(synapser::synGet('syn25428992')$path) %>%
#       select(Biodomain, GO_ID, hgnc_symbol) %>%
#       group_by(Biodomain) %>%
#       summarise(bd_terms = length(unique(GO_ID))) %>%
#       ungroup(),
#     by = 'Biodomain') %>%
#   mutate(pct = 100*(n_term / bd_terms) ) %>%
#   left_join(
#     .,
#     readRDS(synapser::synGet('syn25428992')$path) %>%
#       select(Biodomain, hgnc_symbol, n_hgncSymbol, GO_ID, GOterm_Name ) %>%
#       unnest_longer(hgnc_symbol),
#     by = c('hgnc_symbol','Biodomain') ) %>%
#   full_join(
#     .,
#     # domain labels
#     read_csv(synapser::synGet('syn26856828')$path, col_types = cols()),
#     by = c('Biodomain'='domain')) %>%
#   mutate(Biodomain = case_when(Biodomain == 'none' ~ NA_character_, T ~ Biodomain)) %>%
#   rename(symbol = hgnc_symbol, n_symbol = n_hgncSymbol)
# 
# fst::write.fst(biodom_genes, 'data/biodom_genes.fst')
# 

## Biodom term NW ----

# synapser::synLogin()
# biodom = readRDS(synapser::synGet('syn25428992')$path)
# domains = biodom$Biodomain %>% unique() %>% sort()
# dom.cols = read_csv(synapser::synGet('syn26856828')$path, col_types = cols())
# 
# ks = read_csv(synapser::synGet('syn51163744')$path)
# ks1 <- ks %>% filter(goid_1 != goid_2, kappa > 0.5)
# 
#   # # look @ N overlap by kappa score level
#   # x = map_dfr(
#   #   seq(.1,1,.1),
#   #   ~ ks1 %>% arrange(desc(kappa)) %>% filter(kappa <= .x) %>% slice(1:5)
#   # ) %>%
#   #   rowwise() %>%
#   #   mutate(
#   #     n_genes_1 = biodom$n_symbol[biodom$GO_ID == goid_1] %>% unique(),
#   #     n_genes_2 = biodom$n_symbol[biodom$GO_ID == goid_2] %>% unique(),
#   #     n_overlap = length(intersect(
#   #       biodom %>% filter(GO_ID == goid_1) %>% slice(1) %>%
#   #         pull(symbol) %>% unlist() %>% unique(),
#   #       biodom %>% filter(GO_ID == goid_2) %>% slice(1) %>%
#   #         pull(symbol) %>% unlist() %>% unique()
#   #     ))
#   #   )
# 
# ks1 <- ks %>%
#   left_join(., biodom %>% select(goid_1 = GO_ID, n_gene1 = n_symbol), by = 'goid_1') %>%
#   left_join(., biodom %>% select(goid_2 = GO_ID, n_gene2 = n_symbol), by = 'goid_2') %>%
#   distinct() %>%
#   filter(goid_1 != goid_2,
#          kappa > 0.1,
#          n_gene1 > 3,
#          n_gene2 > 3)
# 
# gen <- read_csv(synapser::synGet('syn45824969')$path)
# omic <- read_csv(synapser::synGet('syn45824835')$path)
# pro <- read_csv(synapser::synGet('syn45824870')$path)
# rna <- read_csv(synapser::synGet('syn45824913')$path)
# trs <- read_csv(synapser::synGet('syn45824995')$path)
# 
# nw <- biodom %>%
#   select(n1 = GO_ID, n2 = Biodomain) %>%
#   mutate(kappa = 1) %>%
#   distinct() %>%
#   bind_rows(
#     ks1 %>% select(n1 = goid_1, n2 = goid_2, kappa),
#     . )
# 
# g = igraph::graph_from_data_frame(nw, directed = F)
# 
# v = tibble( n = igraph::V(g)$name ) %>%
#   left_join(., biodom %>% select(GO_ID, term = GOterm_Name, n_symbol), by = c('n'='GO_ID')) %>%
#   left_join(., dom.cols, by = c('n'='domain')) %>%
#   mutate(term = if_else(is.na(term), n, term),
#          term_size = if_else( is.na(n_symbol), max(n_symbol, na.rm =T), n_symbol ),
#          term_size = rank(term_size, ties.method = 'first'),
#          term_size = term_size / max(term_size),
#          term_size = if_else( n == term, term_size*3, term_size),
#          bd = if_else( n == term, T, F),
#          color = case_when( is.na(color) ~ 'grey80', T ~ color)
#   ) %>%
#   filter(!duplicated(n))
# 
# v <- v %>%
#   left_join(.,
#             gen %>% filter(padj <= 0.05) %>% select(n = ID, gen_NES = NES, gen_padj = padj)) %>%
#   left_join(.,
#             omic %>% filter(padj <= 0.05) %>% select(n = ID, omic_NES = NES, omic_padj = padj)) %>%
#   left_join(.,
#             pro %>% filter(padj <= 0.05) %>% select(n = ID, pro_NES = NES, pro_padj = padj)) %>%
#   left_join(.,
#             rna %>% filter(padj <= 0.05) %>% select(n = ID, rna_NES = NES, rna_padj = padj)) %>%
#   left_join(.,
#             trs %>% filter(padj <= 0.05) %>% select(n = ID, trs_NES = NES, trs_padj = padj))
# 
# v <- distinct(v)
# 
# v <- v %>%
#   mutate(
#     across(contains('padj'), ~if_else(is.na(.x), 1, .x)),
#     across(contains('NES'), ~if_else(is.na(.x), 0, .x))
#     )
# 
# nrow(v) == length(V(g))
# 
# for(i in 2:ncol(v)){
#   igraph::vertex_attr(g, names(v)[i], index = igraph::V(g)) <- v %>% pull(i)
# }
# 
# term_graph = g
# 
# igraph::write.graph(g, 'AD_biodomain_explorer/data/term_graph.graphml', format = 'graphml')
# saveRDS(g, 'AD_biodomain_explorer/data/term_graph.rds')



# # bd terms to add ---------------------------------------------------------
# 
# bd.e <- biodom_genes %>% filter(!is.na(GO_ID), !is.na(Biodomain)) %>% 
#   select(n1 = GO_ID, n2 = Biodomain) %>% distinct() 
# 
# bd.v <- biodom_genes %>% filter(!is.na(Biodomain)) %>% 
#   select(term = Biodomain) %>% distinct() %>% 
#   mutate(
#     ont = NA,
#     size = NA,
#     n_le_genes = NA,
#     NES = NA,
#     padj = NA,
#     `Mitochondrial Metabolism` = NA,
#     none = NA,
#     Synapse = NA,
#     `Structural Stabilization` = NA,
#     `Immune Response` = NA,
#     `Oxidative Stress` = NA,
#     `Lipid Metabolism` = NA,
#     `Cell Cycle` = NA,
#     Proteostasis = NA,
#     Endolysosome = NA,
#     Apoptosis = NA,
#     Vasculature = NA,
#     `Tau Homeostasis` = NA,
#     Myelination = NA,
#     `Metal Binding and Homeostasis` = NA,
#     `APP Metabolism` = NA,
#     Autophagy = NA,
#     `RNA Spliceosome` = NA,
#     Epigenetic = NA,
#     `DNA Repair` = NA,
#     name = term
#   ) %>% pivot_longer(cols = everything()) %>% 
#   group_by(name) %>% summarise(value = list(value)) %>% 
#   pull(value, name = name)

# y = tibble( na = vertex.attributes(trs.graph) ) %>%
#   t() %>% as_tibble(rownames = NA, .name_repair = 'unique') %>% unnest(everything()) %>%
#   rename_with(., ~names(vertex.attributes(trs.graph)), everything())

# # risk graph --------------------------------------------------------------
# 
# synapser::synLogin()
# trs.graph <- synGet('syn51387897')$path %>% igraph::read.graph(., format = 'graphml')
# trs.graph <- igraph::add.vertices(trs.graph, nv = 19, attr = bd.v)
# trs.graph <- igraph::add.edges( trs.graph, 
#                                 bd.e %>% filter(n1 %in% names(V(x))) %>% select(n1,n2) %>% 
#                                   pivot_longer(everything()) %>% pull(value) , 
#                                 attr = list('kappa'=1))
# saveRDS(trs.graph, 'trs_graph.rds')

# # genetics graph ----------------------------------------------------------
# 
# synapser::synLogin()
# gen.graph <- synGet('syn51317242')$path %>% igraph::read.graph(., format = 'graphml')
# gen.graph <- igraph::add.vertices(gen.graph, nv = 19, attr = bd.v)
# gen.graph <- igraph::add.edges( gen.graph, 
#                                 bd.e %>% filter(n1 %in% names(V(x))) %>% select(n1,n2) %>% 
#                                   pivot_longer(everything()) %>% pull(value) , 
#                                 attr = list('kappa'=1))
# saveRDS(gen.graph, 'gen_graph.rds')

# # txome graph -------------------------------------------------------------
# 
# synapser::synLogin()
# txome.graph <- synGet('syn51387896')$path %>% igraph::read.graph(., format = 'graphml')
# txome.graph <- igraph::add.vertices(txome.graph, nv = 19, attr = bd.v)
# txome.graph <- igraph::add.edges( txome.graph, 
#                                 bd.e %>% filter(n1 %in% names(V(x))) %>% select(n1,n2) %>% 
#                                   pivot_longer(everything()) %>% pull(value) , 
#                                 attr = list('kappa'=1))
# saveRDS(txome.graph, 'txome_graph.rds')

# # protx graph -------------------------------------------------------------
# 
# synapser::synLogin()
# pxome.graph <- synGet('syn51317244')$path %>% igraph::read.graph(., format = 'graphml')
# pxome.graph <- igraph::add.vertices(pxome.graph, nv = 19, attr = bd.v)
# pxome.graph <- igraph::add.edges( pxome.graph, 
#                                 bd.e %>% filter(n1 %in% names(V(x))) %>% select(n1,n2) %>% 
#                                   pivot_longer(everything()) %>% pull(value) , 
#                                 attr = list('kappa'=1))
# saveRDS(pxome.graph, 'pxome_graph.rds')


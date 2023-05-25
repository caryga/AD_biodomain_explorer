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

####
## TAD scores ----
####
# 
# synapser::synLogin()
# scores <- synapser::synTableQuery('select * from syn25575156')$filepath %>% read_csv() %>%
#   select(ENSG, GeneName, Overall, rank= Overall_rank, Omics = OmicsScore,
#          Genetics = GeneticsScore, Literature = LiteratureScore,
#          Neuropath = NeuropathScore,
#          Druggability = SM_Druggability_bucket,
#          Safety = safety_bucket,
#          Feasibility = feasibility_bucket) %>%
#   # rowwise() %>%
#   mutate(
#     Druggability = 2*(max(Druggability, na.rm=T) + 1 - Druggability) / max(Druggability, na.rm=T),
#     Safety = 2*(max(Safety, na.rm=T) + 1 - Safety) / max(Safety, na.rm=T),
#     Feasibility = 2*(max(Feasibility, na.rm=T) + 1 - Feasibility) / max(Feasibility, na.rm=T),
#   ) %>%
#   pivot_longer(cols = c(Overall, Omics, Genetics, Literature, Neuropath, Druggability, Safety, Feasibility)) %>%
#   mutate(
#     name = fct_relevel(name, c('Feasibility','Safety','Druggability','Neuropath','Literature','Omics', 'Genetics','Overall')),
#     value = case_when( value == 0 ~ NA_real_, T ~ value)
#     ) %>%
#   arrange(name)
# 
# # saveRDS(scores, 'data/scores.rds')
# fst::write.fst(scores, 'data/scores.fst')
# 
# ####
# ## DE meta-analyses ----
# ####
# 
# omics <- read_csv( synapser::synTableQuery('select * from syn22758536')$filepath ) %>%
#     select(GName, RNA_TE, RNA_fdr_CorPVal, Pro_TE, Pro_fdr_CorPVal) %>%
#     pivot_longer(
#         c(RNA_TE, Pro_TE), names_sep = '_', names_to = c('x', 'y'), values_to = 'TE'
#     ) %>%
#     pivot_longer(
#         c(RNA_fdr_CorPVal, Pro_fdr_CorPVal), names_sep = '_', names_to = c('a', 'b'), values_to = 'fdr'
#     ) %>%
#     distinct() %>% filter(x == a) %>% select(-y, -a, -b) %>%
#     mutate(x = str_replace_all(x, 'Pro','Protein'))
# 
# fst::write.fst(omics, 'data/omics.fst')
# 
# 
# ####
# ## Genetic Evidence ----
# ####
# 
# synapser::synLogin()
# gen <- synapser::synTableQuery('SELECT * FROM syn26844312')$filepath %>% read_csv() %>%
#     filter(!is.na(GeneName), !duplicated(GeneName)) %>%
#     mutate(
#       meanRank_qtlFDR = if_else( meanRank_qtlFDR == 0, NA_real_, meanRank_qtlFDR)
#        )%>%
#     select(GeneName, score_rank,
#            meanRank_gwasP, meanRank_qtlFDR,
#            coding_variant_summary, noncoding_variant_summary,
#            Hsap_pheno_score, Ortholog_pheno_score
#            ) %>%
#     pivot_longer(-c(GeneName,score_rank)) %>%
#     mutate( name = fct_relevel(name,
#                                c('Ortholog_pheno_score',
#                                  'Hsap_pheno_score',
#                                  'noncoding_variant_summary',
#                                  'coding_variant_summary',
#                                  'meanRank_qtlFDR',
#                                  'meanRank_gwasP')),
#             value = case_when( value == 0 ~ NA_real_, T ~ value)) %>%
#     arrange( name )
# 
# # saveRDS(gen, 'data/gen.rds')
# fst::write.fst(gen, 'data/genetics.fst')
# 
# ####
# ## Gene-Biodomain mappings ----
# ####
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
# ####
# ## Biodom term NW ----
# ####
# 
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
#          color = case_when( is.na(color) ~ 'grey80', T ~ color)
#   ) %>%
#   filter(!duplicated(n))
# 
# for(i in 2:ncol(v)){
#   igraph::vertex_attr(g, names(v)[i], index = igraph::V(g)) <- v %>% pull(i)
# }
# 
# term_graph = g
# 
# # igraph::write.graph(g, 'data/term_graph.graphml', format = 'graphml')
# saveRDS(g, 'data/term_graph.rds')
# 
# ####
# ## SEA-AD snRNAseq data ----
# ####
# 
# seaad <- read_csv( synapser::synGet('syn51197803')$path, col_types = cols() ) %>%
#   filter( mean_exp > 0 ) %>%
#   mutate(
#     group = factor(group, c( 'High',  'Intermediate','Low', 'Not AD', 'Dementia', 'No dementia')),
#     broad1 = as.factor(broad) %>% as.numeric(),
#     broad2 = case_when(group == 'Dementia' ~ broad1+.25,
#                         group == 'No dementia' ~ broad1-.25,
#                         group == 'High' ~ broad1+.25,
#                         group == 'Intermediate' ~ broad1+.125,
#                         group == 'Low' ~ broad1-.125,
#                         group == 'Not AD' ~ broad1-.25),
#     et = case_when( group %in% c('No dementia','Dementia') ~ 'dem',
#                     group %in% c('Not AD','Low','Intermediate','High')~'path' ),
#     label = NA_character_
#     ) %>%
#   arrange(desc(mean_exp), desc(fraction_expressed))
# 
# seaad <- seaad %>% 
#   arrange(gene) %>% 
#   group_by(gene) %>% 
#   arrange(desc(mean_exp), desc(fraction_expressed), .by_group = T)
# 
# gene_idx = tibble(gene = unique(seaad$gene)) %>% 
#   rowwise() %>% 
#   mutate(
#     start = which(seaad$gene == gene) %>% min(),
#     end = which(seaad$gene == gene) %>% max(),
#     f_idx = NA,
#     f_start = NA,
#     f_end = NA
#     )
# 
# # split up and save to local data/ dir
# idx = 1
# spl = 2000
# future::plan('multisession')
# for(i in seq(1, nrow(gene_idx), spl)){
#   
#   g_st = i
#   g_end = i+(spl-1)
#   if(g_end > nrow(gene_idx)){ g_end = nrow(gene_idx) }
#   
#   f_st = gene_idx$start[g_st]
#   f_end =  gene_idx$end[g_end]
#   if(f_end > nrow(seaad)){ f_end = nrow(seaad) }
#   
#   sub = seaad[f_st:f_end,]
#   gene_idx$f_idx[g_st:g_end] <- idx
#   gene_idx$f_start[g_st:g_end] <- furrr::future_map_dbl(
#     gene_idx$gene[g_st:g_end], ~ which(sub$gene == .x) %>% min())
#   gene_idx$f_end[g_st:g_end] <- furrr::future_map_dbl(
#     gene_idx$gene[g_st:g_end], ~ which(sub$gene == .x) %>% max())
#   
#   fst::write.fst(sub, paste0('data/seaad_',idx,'.fst'))
#   idx = idx+1
#   
# }
# 
# fst::write.fst(gene_idx, 'data/seaad_gene_idx.fst')
# 
# ####
# ## pre-calc SEA-AD violin plot ----
# ####
# 
# compute_density = function(x, from, to, bw = 'nrd0', adjust = 1, kernel = 'gaussian'){
#   nx <- length(x)
#   w <- rep(1 / nx, nx)
#   dens <- stats::density(x, weights = w, bw = 'nrd0', adjust = 5,
#                          kernel = kernel, n = 512, from = from, to = to)
#   tibble(
#     x = dens$x,
#     density = dens$y,
#     scaled =  dens$y / max(dens$y, na.rm = TRUE),
#     ndensity = dens$y / max(dens$y, na.rm = TRUE),
#     count =   dens$y * nx,
#     n = nx,
#     .size = length(dens$x)
#   )}
# 
# range <- range(seaad$mean_exp, na.rm = TRUE)
# modifier <- 0
# bw = stats::bw.nrd0(seaad$mean_exp)
# pdat <- seaad %>%
#   group_by(broad, broad1) %>%
#   do(compute_density(.$mean_exp,
#                      from = range[1] - modifier*bw,
#                      to = range[2] + modifier*bw,
#                      bw = bw,
#                      adjust = 5)) %>%
#   rename( mean_exp = x )
# 
# pdat <- rbind(pdat %>% mutate(broad1 = broad1 + scaled/2),
#               pdat %>% mutate(broad1 = broad1 - scaled/2))
# 
# base_vp = ggplot(pdat, aes(broad1, mean_exp)) +
#   geom_polygon(color = 'grey65', alpha = .3, aes(fill = broad)) +
#   theme(legend.position = 'none')+
#   # guides(fill = 'none', color = 'legend', size = 'legend' )+
#   viridis::scale_fill_viridis(discrete = T, guide ='none',
#                               begin = 0, end = .8,
#                               option = 'D') +
#   scale_x_continuous(breaks = 1:8,
#                      labels = seaad$broad %>% unique() %>% sort()) +
#   labs(x = '', y = 'mean expression, ln(UP10K+1)' )+
#   theme(axis.title.x = element_blank())
# 
# saveRDS(base_vp, 'data/base_violin.rds')
# 
# rm(pdat, range, modifier, bw, compute_density)
# 
# ####
# ## CRISPRbrain simple phenotype screens ----
# ####
# 
# f = list.files( paste0('~/treatAD_hypothesis','/data/crispr_brain_data/'), pattern = '.csv')
# simple <- map_dfr(
#     f %>% str_subset('RNA|CROP', negate = T),
#     ~ {
#         read_csv( paste0('~/treatAD_hypothesis','/data/crispr_brain_data/',.x)) %>%
#             mutate(cell = .x %>% str_split_fixed(.,'-(?!M)', 2) %>% as_tibble() %>% pull(1),
#                    mode = .x %>% str_extract('CRISPR[ain]'),
#                    expt = .x %>% str_remove_all('-CRISPR[ain].csv') %>%
#                        str_split_fixed('-(?!M)',2) %>% as_tibble() %>% pull(2) ) %>%
#             relocate(cell, mode, expt)
#     }
# ) %>%
#     rename(pval = `P Value`, pheno = `Phenotype`, gene_score = `Gene Score`) %>%
#     arrange(desc(abs(gene_score))) %>%
#     mutate(x = paste0(mode, ': ', TSS ))
# 
# simple <- tibble( expt = simple$expt %>% unique() %>% sort(),
#                   phenotype = c('Expansion (CD34 Staining)',
#                                 'Reactive Oxygen Species (CellRox Intensity)',
#                                 'Proliferation (CFSE Staining)',
#                                 'Survival (14 day)', 'Survival (21 day)', 'Survival (28 day)',
#                                 'Labile Iron (FeRhoNox Intensity)',
#                                 'Immune Activation (CD38 levels)',
#                                 'Lysosome Exocytosis (cell surface LAMP1), +IL1a +TNF +C1q',
#                                 'Lysosome Exocytosis (cell surface LAMP1), Vehicle',
#                                 'Peroxidized Lipids (Liperfluo intensity)',
#                                 'Lysosome (LysoTracker intensity)',
#                                 'Lysosome Mass or pH (LysoTracker staining), +IL1a +TNF +C1q',
#                                 'Lysosome Mass or pH (LysoTracker staining), Vehicle',
#                                 'Survival, no antioxidants',
#                                 'Phagocytosis (pHRodo-rat synaptosomes)',
#                                 'Phagocytosis (pHRodo-labeled synaptosomes), +IL1a +TNF +C1q',
#                                 'Phagocytosis (pHRodo-labeled synaptosomes), Vehicle',
#                                 'Survival, PSAP KO & no antioxidants',
#                                 'Survival, PSAP KO',
#                                 'Survival',
#                                 'Survival-Proliferation',
#                                 'Inflammatory activation (cell-surface VCAM1 levels), +IL1a +TNF +C1q')
# ) %>% left_join(simple, . , by = 'expt')
# 
# fst::write.fst(simple, 'data/crispr_brain_simple_phenos.fst')
# 
# 
# # deprecated --------------------------------------------------------------
# 
# ####
# ## synapse upload ----
# ####
# # 
# # synapser::synLogin()
# # files = list.files('data', full.names = T) %>% str_subset('allen|_v1', negate = T)
# # for(f in files){
# #   foo <- synapser::synStore( synapser::File(f, parent = 'syn51174314') )
# # }
# 
# ####
# ## populate data/ dir with synapse download ----
# ####
# # 
# # synapser::synLogin()
# # files = synapser::synGetChildren('syn51174314')$asList() %>%
# #   tibble(f = .) %>% unnest_wider(f)
# # syn.date = str_split_fixed(files$modifiedOn, 'T', 2)[,1] %>% sort(decreasing = T) %>% .[1]
# # if( !dir.exists('data') ){ dir.create('data')  }
# # if( length(list.files('data')) == 0 ){
# #   for(idx in 1:nrow(files)){
# #     foo = synapser::synGet( files$id[idx], downloadLocation = paste0(here::here(), '/data/'))
# #   }
# # } else {
# #   dir.date = file.info( list.files('data', full.names = T) ) %>% pull(mtime) %>% 
# #     str_split_fixed(string = ., pattern = ' ', n = 2) %>% .[,1] %>% sort(decreasing = T) %>% .[1]
# #   if( syn.date > dir.date ){
# #     for(idx in 1:nrow(files)){
# #       foo = synapser::synGet( files$id[idx], downloadLocation = paste0(here::here(), '/data/'))
# #     }
# #   }
# # }
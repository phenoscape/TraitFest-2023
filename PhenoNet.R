
#########################
### PhenoNet pipeline ###
#########################


# Brainstorming and tentative tutorial for reconstructing the evolution of complexity in fishes using data from the Phenoscape KB.

# The dataset used was that of Mirande et al (2019).


# STEP 1. Loading the packages.
library("rphenoscape")
library("rphenoscate")
library("igraph")


# STEP 2. Assembling the data set from a given study.

# Retrieve list of phylogenetic studies.
studies <- get_studies()

# Get a particular study # (Change this part to get a particular study).
study <- studies$id[studies$label == 'Mirande (2019)']

# Get NeXML data.
selected_study <- get_study_data(study)

# Build the original character matrix.
char_mat <- RNeXML::get_characters(selected_study[[1]])

# Get rownames and colnames from data set.
row_mat <- rownames(char_mat)
col_mat <- colnames(char_mat)

# Check rownames and colnames from data set.
row_mat[1:3]
col_mat[1:3]

# OTU labels are not available.
# Get metadata the original character matrix.
selected_study_meta <- get_char_matrix_meta(selected_study[[1]])

# Get new rownames.
row_mat <- selected_study_meta$id_taxa$label[
  match(rownames(char_mat), selected_study_meta$id_taxa$otu)]

# Remove duplicated species from matrix.
char_mat <- char_mat[-(which(row_mat == "Moenkhausia sanctaefilomenae")),]
row_mat <- row_mat[-(which(row_mat == "Moenkhausia sanctaefilomenae"))]
rownames(char_mat) <- row_mat


# STEP 3. Getting all semantic phenotypes from a given study.

# Retrieve all semantic phenotypes from the Phenoscape KB.
# Get all phenotype data from the study.
phenotypes <- get_phenotypes(study = study)

# Convert data to a phenotype object (Warning: ~10 min).
selected_study_obj <- as.phenotype(phenotypes, withTaxa = TRUE)

# Get all anatomical entities from the study.
ent <- unlist(sapply(selected_study_obj, function(x) x$eqs$entities ))


# STEP 4. Get dependency matrix and build a graph.

# Get dependency matrix.
ent2 <- unique(ent)
depmat <- pa_dep_matrix(ent2, .names = "label")

# Get a graph.
G <- depmat
diag(G) <- NA
G <- graph_from_adjacency_matrix(rphenoscate:::remove_indirect(t(as.matrix(G))))

# Plot graph.
plot(G, vertex.size = 3, edge.arrow.size = 0.5, vertex.label.cex = 0.5)

pdf(file = "graph.pdf")
plot(G, vertex.size = 3, edge.arrow.size = 0.5, vertex.label.cex = 0.5)
dev.off()

# Fancy graph.
x <- as.matrix(G)

# Remove indirect dependancies.
x <- rphenoscate:::remove_indirect(x)

# Remove NAs.
x[is.na(x)] <- 0

g1 <- graph_from_adjacency_matrix(x)

# Vertex size proportional to number of dependancies.
V(g1)$size <- (rowSums(x) + 1) / (max(rowSums(x))/8) # the divider here controls the relative size difference of the vertices. 

# Vertex colour based on shortest path distance from root.
pd <- shortest.paths(g1, 1)
V(g1)$color <- pd
g1$palette <- hcl.colors(max(pd))

# Plot.
plot(g1, edge.arrow.size=0.5, vertex.label = 1:nrow(x), vertex.label.cex = 0.4, vertex.label.color= "black")


# STEP 5a. Extract and plot subgraphs (better trying with OntoTrace!---see below).
char_mat2 <- matrix(0, nrow(char_mat), length(unique(ent)))
rownames(char_mat2) <- rownames(char_mat)
colnames(char_mat2) <- colnames(depmat)

for(i in 1:length(selected_study_obj)){
	v_no <- which(ent2 %in% selected_study_obj[[i]]$eqs$entities)
	tax <- which(rownames(char_mat2) %in% selected_study_obj[[i]]$taxa$label)
	char_mat2[tax, v_no] <- 1
}
	
taxon_x <- 1 
# Plot a subgraph.
g2 <- subgraph(g1, which(char_mat2[taxon_x,] == 1))
g2_m <- t(as_adjacency_matrix(g2, sparse = F))
V(g2)$size <- (rowSums(g2_m) + 4) / (max(rowSums(g2_m))/10) # the divider here controls the relative size difference of the vertices. 

plot(g2, edge.arrow.size=0.5, vertex.label = which(char_mat2[taxon_x,] == 1), vertex.label.cex = 0.4, vertex.label.color= "black")	
	
# Testing the connectivity of subgraphs.
g2 <- subgraph(g1, which(char_mat2[taxon_x,] == 1))
disconected <- setdiff(V(g2), subcomponent(g2, 1, mode = "all")) #ID of disconected nodes in g2
disconected_char <- which(colnames(char_mat2) %in% names(V(g2)[setdiff(V(g2), subcomponent(g2, 1, mode = "all"))])) #ID of disconected characters in g2

missing_nodes <- numeric()
for(j in disconected_char){
	missing_nodes <- c(missing_nodes, subcomponent(g1, j, mode = "out"))
}
missing_nodes <- unique(missing_nodes)

taxon_x_nodes <- which(char_mat2[i,] == 1)
missing_nodes <- setdiff(missing_nodes, taxon_x_nodes) # remove upstream nodes that are already coded as prsent (e.g. multicelularity)


# Add missing nodes.
char_mat2[taxon_x,missing_nodes] <- 1

# Replot.
g2 <- subgraph(g1, which(char_mat2[taxon_x,] == 1))
g2_m <- t(as_adjacency_matrix(g2, sparse = F))
V(g2)$size <- (rowSums(g2_m) + 4) / (max(rowSums(g2_m))/10) # the divider here controls the relative size difference of the vertices. 
plot(g2, edge.arrow.size=0.5, vertex.label = which(char_mat2[taxon_x,] == 1), vertex.label.cex = 0.4, vertex.label.color= "black")	
	

############################################################################

# STEP 5b. Rework with OntoTrace.

# Call ontotrace.
onto <- rphenoscape::get_ontotrace_data(taxon = rownames(char_mat), entity = colnames(depmat), subsume = F, variable_only = F)

# Get char matrix.
char_mat3 <-  RNeXML::get_characters(onto)

# Get list of indices for polymorphic cells.
poly <- which(char_mat3 == "1 and 0", arr.ind = T)
poly <- rbind(poly, which(char_mat3 == "0 and 1", arr.ind = T))

char_mat4 <- char_mat3 #backup

for(i in 1:nrow(poly)){

    char_mat4[poly[i,1], poly[i, 2]] <- 1
}

char_mat4 <- type.convert(char_mat4)
char_mat4 <- as.matrix(char_mat4)

# Lots of NAs. Could do ASR to guess presence/absence. We will just infer everything as present.
char_mat4[is.na(char_mat4)] <- 1

# Get subgraphs.
tax_subgraphs <- list()
for(i in 1:nrow(char_mat4)){
    g2 <- subgraph(g1, which(char_mat4[i,] == 1))
    g3_v <- subcomponent(g2, 1, mode = "all")
    tax_subgraphs[[i]] <- subgraph(graph = g2, v = g3_v)
}

v_per_taxa <- unlist(lapply(tax_subgraphs, function(x) vcount(x)))

# Get decent looking taxa.

tax_subgraphs2 <- tax_subgraphs[which(v_per_taxa > 140)]

##################################################################

# STEP 6. Plotting evolution of complexity.

# Get a timetree for the fish taxa. Only 130 or so species.
timetree <- fishtree_phylogeny(species = rownames(char_mat4)[which(v_per_taxa > 140)], type = "chronogram")

# Get graphs for species which are in the tree.

# List of species with > 140 vertices in graphs (names fixed to match timetree).
decent_species <- c("Acestrocephalus_sardina", "Acestrorhynchus_lacustris", "Acestrorhynchus_pantaneiro", 
"Acrobrycon_ipanquianus", "Agoniates_anchovia", "Agoniates_halecinus", 
"Apareiodon_affinis", "Aphyocharacidium_bolivianum", "Aphyocharax_anisitsi", 
"Aphyocharax_dentatus", "Aphyocharax_nattereri", "Aphyodite_grammica", 
"Argopleura_magdalenensis", "Astyanacinus_moorii", "Astyanax_abramis", 
"Astyanax_aurocaudatus", "Astyanax_bransfordii", "Astyanax_cf._rutilus", 
"Astyanax_chico", "Astyanax_correntinus", "Astyanax_cremnobates", 
"Astyanax_eigenmanniorum", "Astyanax_endy", "Astyanax_intermedius", 
"Astyanax_lacustris", "Astyanax_latens", "Astyanax_lineatus", 
"Astyanax_mexicanus", "Astyanax_paris", "Astyanax_pellegrini", 
"Astyanax_puka", "Astyanax_troya", "Atopomesus_pachyodus", "Attonitus_ephimeros", 
"Aulixidens_eugeniae", "Axelrodia_lindeae", "Bario_steindachneri", 
"Boulengerella_lateristriga", "Brachychalcinus_copei", "Brittanichthys_axelrodi", 
"Brycon_falcatus", "Brycon_gouldingi", "Brycon_orbignyanus", 
"Brycon_pesu", "Brycon_polylepis", "Bryconaethiops_macrops", 
"Bryconamericus_agna", "Bryconamericus_cf._iheringii", "Bryconamericus_exodon", 
"Bryconamericus_heteresthes", "Bryconamericus_iheringii", "Bryconamericus_ikaa", 
"Bryconamericus_mennii", "Bryconamericus_rubropictus", "Bryconamericus_stramineus", 
"Bryconamericus_uporas", "Bryconella_pallidifrons", "Bryconexodon_juruenae", 
"Bryconops_affinis", "Bryconops_melanurus", "Carlana_eigenmanni", 
"Carnegiella_strigata", "Chalceus_macrolepidotus", "Characidium_rachovii", 
"Characidium_zebra", "Charax_leticiae", "Charax_stenopterus", 
"Cheirodon_interruptus", "Compsura_heterura", "Creagrutus_anary", 
"Creagrutus_atrisignum", "Creagrutus_meridionalis", "Creagrutus_taphorni", 
"Ctenobrycon_hauxwellianus", "Cyanogaster_noctivaga", "Cynopotamus_argenteus", 
"Cyphocharax_spilotus", "Deuterodon_iguape", "Deuterodon_langei", 
"Deuterodon_longirostris", "Deuterodon_potaroensis", "Deuterodon_rosae", 
"Deuterodon_sp.", "Deuterodon_stigmaturus", "Diapoma_alburnum", 
"Diapoma_guarani", "Diapoma_obi", "Diapoma_speculiferum", "Diapoma_terofali", 
"Diapoma_uruguayense", "Distichodus_maculatus", "Ectrepopterus_uruguayensis", 
"Engraulisoma_taeniatum", "Eretmobrycon_emperador", "Eretmobrycon_peruanus", 
"Eretmobrycon_scleroparius", "Erythrocharax_altipinnis", "Exodon_paradoxus", 
"Galeocharax_gulo", "Galeocharax_humeralis", "Gnathocharax_steindachneri", 
"Grundulus_bogotensis", "Gymnocharacinus_bergii", "Gymnocorymbus_ternetzi", 
"Hasemania_crenuchoides", "Hasemania_nana", "Hemibrycon_caucanus", 
"Hemibrycon_dariensis", "Hemibrycon_surinamensis", "Hemigrammus_aguaruna", 
"Hemigrammus_bleheri", "Hemigrammus_cf._lunatus", "Hemigrammus_erythrozonus", 
"Hemigrammus_haraldi", "Hemigrammus_marginatus", "Hemigrammus_pulcher", 
"Hemigrammus_rodwayi", "Hemigrammus_ulreyi", "Hemigrammus_unilineatus", 
"Hemiodus_thayeria", "Heterocharax_macrolepis", "Heterocheirodon_yatai", 
"Hollandichthys_multifasciatus", "Hoplias_cf._malabaricus", "Hoplocharax_goethei", 
"Hyphessobrycon_anisitsi", "Hyphessobrycon_bifasciatus", "Hyphessobrycon_boulengeri", 
"Hyphessobrycon_compressus", "Hyphessobrycon_elachys", "Hyphessobrycon_eques", 
"Hyphessobrycon_erythrostigma", "Hyphessobrycon_herbertaxelrodi", 
"Hyphessobrycon_loweae", "Hyphessobrycon_megalopterus", "Hyphessobrycon_meridionalis", 
"Hyphessobrycon_moniliger", "Hyphessobrycon_montagi", "Hyphessobrycon_poecilioides", 
"Hyphessobrycon_pulchripinnis", "Hyphessobrycon_santae", "Hyphessobrycon_socolofi", 
"Hyphessobrycon_vanzolinii", "Hyphessobrycon_wajat", "Iguanodectes_geisleri", 
"Inpaichthys_kerri", "Jupiaba_acanthogaster", "Jupiaba_anteroides", 
"Jupiaba_mucronata", "Jupiaba_polylepis", "Jupiaba_scologaster", 
"Knodus_breviceps", "Knodus_gamma", "Knodus_meridae", "Knodus_moenkhausii", 
"Knodus_tanaothoros", "Leporinus_striatus", "Lonchogenys_ilisha", 
"Macropsobrycon_uruguayanae", "Markiana_nigripinnis", "Metynnis_mola", 
"Micralestes_stormsi", "Microgenys_minuta", "Microschemobrycon_casiquiare", 
"Mimagoniates_rheocharis", "Moenkhausia_bonita", "Moenkhausia_ceros", 
"Moenkhausia_collettii", "Moenkhausia_cotinho", "Moenkhausia_dichroura", 
"Moenkhausia_heikoi", "Moenkhausia_jamesi", "Moenkhausia_lata", 
"Moenkhausia_lepidura", "Moenkhausia_pirauba", "Moenkhausia_xinguensis", 
"Myxiops_aphos", "Nantis_indefessus", "Nematocharax_venustus", 
"Odontostilbe_fugitiva", "Odontostilbe_microcephala", "Odontostilbe_paraguayensis", 
"Odontostilbe_pequira", "Odontostilbe_pulchra", "Odontostoechus_lethostigmus", 
"Oligosarcus_bolivianus", "Oligosarcus_itau", "Oligosarcus_jenynsii", 
"Oligosarcus_longirostris", "Oligosarcus_menezesi", "Oligosarcus_paranensis", 
"Oligosarcus_pintoi", "Orthospinus_franciscensis", "Paracheirodon_axelrodi", 
"Parecbasis_cyclolepis", "Parodon_nasus", "Petitella_georgiae", 
"Phenacogaster_franciscoensis", "Phenacogaster_tegatus", "Phenagoniates_macrolepis", 
"Phycocharax_rasbora", "Piabarchus_analis", "Piabina_argentea", 
"Piabina_thomasi", "Piabucus_melanostoma", "Piaractus_mesopotamicus", 
"Poptella_paraguayensis", "Prionobrama_filigera", "Prionobrama_paraguayensis", 
"Pristella_maxillaris", "Probolodus_heterostomus", "Prochilodus_lineatus", 
"Prodontocharax_melanotus", "Psellogrammus_kennedyi", "Pseudochalceus_kyburzi", 
"Pseudocheirodon_terrabae", "Pseudocorynopoma_doriae", "Pseudocorynopoma_heterandria", 
"Pyrrhulina_australis", "Rhaphiodon_vulpinus", "Rhoadsia_altipinna", 
"Roeboexodon_guyanensis", "Roeboides_descalvadensis", "Roeboides_microlepis", 
"Salminus_brasiliensis", "Salminus_hilarii", "Serrapinnus_calliurus", 
"Serrapinnus_microdon", "Serrapinnus_notomelas", "Serrasalmus_maculatus", 
"Serrasalmus_rhombeus", "Spintherobolus_ankoseion", "Stethaprion_erythrops", 
"Stichonodon_insignis", "Tetragonopterus_argenteus", "Tetragonopterus_chalceus", 
"Thayeria_boehlkei", "Thayeria_obliqua", "Thoracocharax_stellatus", 
"Triportheus_angulatus", "Triportheus_nematurus", "Triportheus_pantanensis", 
"Xenagoniates_bondi")

species_ID <- which(decent_species %in% timetree$tip.label)
tax_subgraphs3 <- tax_subgraphs2[species_ID]
names(tax_subgraphs3) <- decent_species[species_ID]
tax_subgraphs3 <- tax_subgraphs3[match(timetree$tip.label, names(tax_subgraphs3))]

# Number of dependencies per taxon = crude measure of phenotype complexity.
n_edges <- unlist(lapply(tax_subgraphs3, function(x) ecount(x)))

# Number of verticies (anatomical entities) per taxon = crude measure of phenotype complexity.
n_vertices <- unlist(lapply(tax_subgraphs3, function(x) vcount(x)))

# Mean minimum path distance to the root character = measure of nested dependancies = ???
r_centrality <- unlist(lapply(tax_subgraphs3, function(x) (sum(distances(x, v = 1, to = V(x), weights = NULL))/(vcount(x) -1))))

# Mean minimum path distance to the root character = measure of nested dependancies = ???
max_centrality <- unlist(lapply(tax_subgraphs3, function(x) max(distances(x, v = 1, to = V(x), weights = NULL))))

# Turns out this is equivalent to n_edges!
# Mean degree distribution of subgraph. Measure of the number of nested dependencies = better measure of complexity.
# First get average in and out the degree distributions for the entire dataset (g1).
g1_in <- degree_distribution(g1, cumulative = F, mode = "in") * vcount(g1)
max_degree <- sum(0:(length(g1_in) - 1) * g1_in)

degree <- unlist(lapply(tax_subgraphs3, function(x){
	d_x <- degree_distribution(x, cumulative = F, mode = "in") * vcount(x)
	(sum(0:(length(d_x) - 1) * d_x))/max_degree
}))
	
# Plot on tree!!!
obj<-contMap(timetree,n_edges, plot = F)
plot(obj, fsize = 0.4)

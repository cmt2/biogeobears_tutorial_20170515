# Instalar el paquete
install.packages("BioGeoBEARS", dependencies=TRUE, repos="http://cran.rstudio.com")
install.packages("snow")
install.packages("roxygen2")

# Cargar el paquete y las dependencias
library(optimx)
library(FD)
library(snow)
library(parallel)
library(BioGeoBEARS)

# Cargar archivos necesarios para correr el paquete
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

# Cargar datos ejemplares 
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
tree_file_name = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep="")) # filogenia ejemplar
tr = read.tree(tree_file_name)

plot(tr) # plotear la filogenia ejemplar
title("Example Psychotria phylogeny from Ree & Smith (2008)")
axisPhylo()


geo_file_name = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep="")) # datos geográficos ejemplares
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geo_file_name)

tipranges # ver los datos geográficos 

####################

# Construir el modelo DEC (modelo 'default' del paquete)

BioGeoBEARS_run_object = define_BioGeoBEARS_run() # empezar un modelo

BioGeoBEARS_run_object$trfn = tree_file_name # dile al modelo el nombre de la filogenia que vas a utilizar
BioGeoBEARS_run_object$geogfn = geo_file_name # dile también el nombre de los datos geográficos

BioGeoBEARS_run_object$max_range_size = 4 # específica el número de zonas geográficas 
BioGeoBEARS_run_object$min_branchlength = 0.000001 # varios otros parámetros 
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Correr el modelo DEC 
results_DEC = bears_optim_run(BioGeoBEARS_run_object)

######################

# Construir el model DEC + J, empezando igual que antes 

BioGeoBEARS_run_object = define_BioGeoBEARS_run() 

BioGeoBEARS_run_object$trfn = tree_file_name # dile de nuevo los nombres de los archivos
BioGeoBEARS_run_object$geogfn = geo_file_name

BioGeoBEARS_run_object$max_range_size = 4 # de nuevo configurar los parámetros, igual que antes
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Ahora empezamos a cambiar la configuración para construir el DEC + J

dstart = results_DEC$outputs@params_table["d","est"] # usar las estimaciones de antes para dar inicio a la búesqueda huerística 
estart = results_DEC$outputs@params_table["e","est"]
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# añadir otro parámetro de J 

jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# correr el modelo

results_DECJ = bears_optim_run(BioGeoBEARS_run_object) 

#################

# Hacer gráficas de los resultados

pdffn = "Psychotria_DEC_vs_DEC+J.pdf"
pdf(pdffn, width=6, height=6) # empezar un pdf 

analysis_titletxt = "DEC on Psychotria" # los resultados del modelo DEC 
results_object = results_DEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text",
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                                cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                         label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

analysis_titletxt ="DEC+J on Psychotria" # los resultados del model DEC + J
results_object = results_DECJ
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text",
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                                cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                         label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off() # terminar el proceso de hacer el pdf y abrirlo
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#################

# Comparar el likelihood de los modelos

LnL_2 = get_LnL_from_BioGeoBEARS_results_object(results_DEC) # extraer el likelihood del DEC 
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(results_DECJ) # extraer el likelihood del DEC + J

numparams1 = 3 # definir el número de parámetros 
numparams2 = 2

stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2) # calcular el AIC 
stats$AIC1
stats$AIC2
stats$pval



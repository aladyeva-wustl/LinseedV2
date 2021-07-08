source('/app/MonteCarloLinseedSisal.R')
source('/app/LinseedMetadata.R')
library(optparse)
library(yaml)

option_list <- list(
  make_option(c("-c", "--config"), action="store", default=NA, type='character',
              help="YAML with configutation path"),
  make_option(c("-i", "--init"), action="store", default=NA, type='integer',
              help="Initialization number")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.na(opt$c)) {
  stop("YAML Configuration file path is not specified")
}

obj <- yaml.load_file(opt$c)
path_ <- file.path('/app', 'results', obj$path)
print(path_)
dir.create(path_, recursive = T)

analysis_name <- obj[['analysis_name']]
if (!is.na(opt$i)) {
  analysis_name <- paste0(analysis_name,"_",opt$i)
}

if (!is.null(obj[['lambda']])) {
  lambda <- obj[['lambda']]
} else {
  lambda <- 0.01
}

if (!is.null(obj[['beta']])) {
  beta <- obj[['beta']]
} else {
  beta <- 0.001
}

if (!is.null(obj[['data']])) {
  data_ <- readRDS(obj[['data']][['path']])
  tmp_mclin <- MonteCarloLinseedSisal$new(dataset = obj[['dataset']], path = path_, analysis_name = analysis_name,
                                                cell_types = obj[['cell_types']], 
                                                optimize_iterations = obj[['optimize_iterations']],
                                                search_iterations = obj[['search_iterations']], 
                                                data = data_, radius = obj[['radius']],
                                                lambda = lambda, beta = beta)
} else {
  tmp_mclin <- MonteCarloLinseedSisal$new(dataset = obj[['dataset']], path = path_, analysis_name = analysis_name,
                                                cell_types = obj[['cell_types']], 
                                                optimize_iterations = obj[['optimize_iterations']],
                                                search_iterations = obj[['search_iterations']], 
                                                radius = obj[['radius']],
                                                lambda = lambda, beta = beta)
}

  tmp_mclin$selectTopGenes(obj[['top_genes']])

  if (!is.null(obj[['svd_k']])) {
    k <- max(obj[['cell_types']],obj[['svd_k']])
  } else {
    k <- obj[['cell_types']]
  }

  print(tmp_mclin$path_)
  
  tmp_mclin$getSvdProjections(k=k)
  tmp_mclin$selectInit()
  tmp_mclin$optimizeInitProportions()
  tmp_mclin$runSearch()
  write.table(tmp_mclin$full_proportions,
            file=paste0(tmp_mclin$path_,"/","markers_",tmp_mclin$analysis_name,"_proportions.tsv"),
            sep="\t",col.names = NA, row.names = T, quote = F)
  colnames(tmp_mclin$W) <- paste0("Cell_type_",1:tmp_mclin$cell_types)
  toSave <- tmp_mclin$W
  toSave <- tmp_mclin$getFoldChange(toSave)
  toSave <- rbind(c(rep(NA,tmp_mclin$cell_types),round(apply(tmp_mclin$full_proportions,1,mean),4)),toSave)
  rownames(toSave) <- c("avg_proportions",rownames(tmp_mclin$filtered_dataset))
  write.table(toSave,file=paste0(tmp_mclin$path_,"/","markers_",tmp_mclin$analysis_name,"_basis_fc.tsv"),
            sep="\t",col.names = NA, row.names = T, quote = F)
  metadata_ <- LinseedMetadata$new(tmp_mclin)
  save(metadata_,file=paste0(metadata_$path_,"/",metadata_$analysis_name,".meta"))


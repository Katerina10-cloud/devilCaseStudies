rm(list = ls())

source("utils/utils.R")
source("utils/edgeR.R")
source("utils/limma.R")
source("utils/glmGamPoi.R")
source("utils/nebula.R")
source("utils/devil.R")
library(Seurat)
library(ggplot2)

list.func <- list(
  glmgp.mult,
  edger.mult,
  limma.mult,
  glmgp.cell.mult,
  nebula.mult,
  devil.base,
  devil.mixed,
  glmgp.cell.fixed
)

set.seed(12345)

MAX_GENE <- 1000
exp.cut <- 0.5
n.sim <- 5

args <- commandArgs(trailingOnly = TRUE)
author <- args[1]

if (!(author %in% c('bca', 'yazar', 'hsc', 'kumar'))){
  stop('author must be bca or yazar or hsc or kumar')
}

if (author == "bca") {
  seurat.obj <- readRDS("datasets/bca.seurat.rds")
} else if (author == "yazar") {
  seurat.obj <- readRDS('datasets/yazar.seurat.rds')
} else if (author == 'hsc') {
  seurat.obj <- readRDS('datasets/hsc.rds')
} else if (author == "kumar") {
  seurat.obj <- readRDS('datasets/kumar.rds')
} else {
  stop('author must be bca or yazar')
}

# prepare data
cols <- c('donor_id', 'cell_type')
col.data <- seurat.obj[[cols]]
col.data$cell_id <- rownames(col.data)
#cnt <- GetAssayData(object = seurat.obj, slot = "counts")
if (author %in% c("bca", "hsc", "kumar")) {
  cnt <- seurat.obj@assays$RNA$counts
} else if (author %in% c("yazar")) {
  cnt <- seurat.obj@assays$RNA$data
} else  {
  stop("author must be bca or yazar")
}

head(col.data)

# sort cell type by numbers, select top 6
ct.used <- col.data %>%
  group_by(cell_type) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  top_n(6) %>%
  pull(cell_type)

# select donors with more than 100 cells per selected cell types
if (author == "yazar") {
  donor.cells <- col.data %>%
    group_by(donor_id, cell_type, .drop=FALSE) %>%
    summarise(n=n()) %>%
    pivot_wider(names_from=cell_type, values_from=n) %>%
    select(ct.used) %>%
    filter(if_all(ct.used,~.>85))
} else {
  donor.cells <- col.data %>%
    group_by(donor_id, cell_type, .drop=FALSE) %>%
    summarise(n=n()) %>%
    pivot_wider(names_from=cell_type, values_from=n) %>%
    select(ct.used) %>%
    filter(if_all(ct.used,~.>50))  
}

donor.cells$cell.sum <- rowSums(donor.cells[,2:7])
donor.used <- donor.cells %>% arrange(-cell.sum) %>% pull(donor_id)

col.data %>%
  group_by(cell_type) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  top_n(6)

cell.used <- col.data$donor_id %in% donor.used
cnt.used <- cnt[rowMeans(cnt) > 0.1, cell.used]
col.data.used <- col.data[cell.used,]

# Simulations ####
timing_df <- dplyr::tibble()
IS.PB <- c(FALSE, TRUE)
INT.CT <- 1:length(ct.used)
I.ITER <- 1:n.sim
PROB.DE <- c(.005, .025, .05)
N.SAMPLE <- c(4, 20)
n.test <- prod(length(IS.PB), length(INT.CT), length(I.ITER), length(PROB.DE), length(N.SAMPLE))
curr.iter <- 1

for (is.pb in IS.PB) {
  for (int.ct in INT.CT) {
    for (i.iter in I.ITER) {
      
      s <- Sys.time()
      
      cell.type <- ct.used[[int.ct]]
      bool.ct <- col.data.used$cell_type == cell.type
      col.data.ct <- col.data.used[bool.ct,]

      # select genes with sufficient mean (nblmm doesn't work for small mean)
      cnt.ct <- cnt.used[,bool.ct]
      cnt.ct.bm <- cnt.ct[rowMeans(cnt.ct) > 0.1,]
      cnt.ct.bm <- cnt.ct.bm[sample(1:nrow(cnt.ct.bm), MAX_GENE),]

      # boolean index of cells
      bool.cell.donor <- lapply(
        unique(col.data.ct$donor_id),
        function(id.donor){
          return(col.data.ct$donor_id == id.donor)
        }
      )
      col.data.donor <- lapply(
        bool.cell.donor,
        function(bool.donor){
          return(col.data.ct[bool.donor,])
        }
      )
      cnt.donor <- lapply(
        bool.cell.donor,
        function(bool.donor){
          return(cnt.ct.bm[,bool.donor])
        }
      )

      for (prob_de in PROB.DE) {
        for (n.sample in N.SAMPLE) {
          message("Iteration = ", curr.iter, ", Percentage = ", round(100 * curr.iter / n.test, 2), "%")
          list.result.null <- list()
          list.result.pow <- list()

	        # sample individuals
          ind.select <- sample.int(
            n=length(unique(col.data.ct$donor_id)),
            size=n.sample,
            replace=TRUE
          )

          # construct col.data & cnt matrix
          col.data.select <- data.table::rbindlist(col.data.donor[ind.select], idcol="id")
          cnt.select <- do.call(cbind, cnt.donor[ind.select])
          colnames(cnt.select) <- rownames(col.data.select)

          # select cells randomly
          n.per.donor <- col.data.select %>% group_by(id, donor_id) %>% summarise(n=n()) %>% pull(n)
          #n.per.donor <- rep(min(n.per.donor), length(n.per.donor))
	        cell.select <- as.logical(unlist(sapply(n.per.donor, rbinom, size=1, p=min(n.per.donor) / n.per.donor)))
          col.data.select <- col.data.select[cell.select,]
          cnt.select <- cnt.select[,cell.select]

          # assign treatment
          col.data.select$tx_cell <- rbinom(n=nrow(col.data.select), size=1, p=0.5)
          if (is.pb){
            # assign treatment label
            n.tx <- as.integer(n.sample/2)
            urn <- c(rep(1,n.tx), rep(0,n.sample-n.tx))
            tx.ind <- sample(x=urn, size=n.sample, replace=FALSE)

            # assign tx to cells
            cell.per.ind <- col.data.select %>%
              group_by(id) %>%
              summarise(n=n())
            col.data.select$tx_cell <- rep(tx.ind, times=cell.per.ind$n)
          }

          # for power simulation, cut expression to half (+ force int)
          max_gene <- as.integer(MAX_GENE * prob_de)
          n.gene <- max_gene
          idx.de <- 1:max_gene
          idx.tx <- which(col.data.select$tx_cell == 1)
          downCells(cnt.select, idx.de, idx.tx, exp.cut)
          cnt.select.t <- t(cnt.select)

          # run tests
          message('start regression')
          # list.data <- list(
          #   cnt.select,
          #   cnt.select,
          #   cnt.select,
          #   cnt.select,
          #   cnt.select,
          #   cnt.select,
          #   cnt.select
          # )
          list.result.method <- list()
          timings <- c()
          cnt.input <- cnt.select %>% as.matrix()
          
          s <- Sys.time()
          for (int.test in 1:length(list.func)){
            print(int.test)
            list.result.method[[int.test]] <- list.func[[int.test]](
              cnt.input, 
              col.data.select
            )[,4:5]

            timings <- c(timings, unique(list.result.method[[int.test]][,2]))
            list.result.method[[int.test]] <- list.result.method[[int.test]][,1]
          }
          e <- Sys.time()

          timing_df <- dplyr::bind_rows(
            timing_df,
            dplyr::tibble(
              algo=c('glmGamPoi (Pb)', 'edgeR (Pb)', 'limma (Pb)', 'glmGamPoi (cell)', 'Nebula', 'Devil (base)', 'Devil (mixed)', "glmGamPoi (fixed)"),
	      timings = timings,
              author=author,
              is.pb=is.pb,
              n.sample=n.sample,
              n.gene=n.gene,
              int.ct=int.ct, n.cells=dim(cnt.select)[2], iter=i.iter))

          df.result <- do.call(rbind, list.result.method)
          rownames(df.result) <- c('glmGamPoi (Pb)', 'edgeR (Pb)', 'limma (Pb)', 'glmGamPoi (cell)', 'Nebula', 'Devil (base)', 'Devil (mixed)', "glmGamPoi (fixed)")
          list.result.null[[i.iter]] <- t(df.result)[(max_gene+1):MAX_GENE,]
          list.result.pow[[i.iter]] <- t(df.result)[1:max_gene,]

          # results
          df.null <- do.call(rbind, list.result.null)
          df.pow <- do.call(rbind, list.result.pow)

          # save name
          path.tail <- paste(author,'n',n.sample,'ngenes',n.gene,'ct',int.ct,'probde',prob_de,'iter',i.iter,'csv',sep='.')
          path.pow.head <- ifelse(is.pb, 'pow_subject/', 'pow_cell/')
          path.null.head <- ifelse(is.pb, 'null_subject/', 'null_cell/')

          write.csv(df.null, paste0(path.null.head, path.tail))
          write.csv(df.pow, paste0(path.pow.head, path.tail))

          saveRDS(timing_df, paste0("timing_results/", author, ".rds"))

          curr.iter <- curr.iter + 1
        }
      }
    }
  }
}

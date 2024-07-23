
method_colors = c(
  "glmGamPoi (Pb)" = "#A22E29",
  "edgeR" = "#7D629E",
  "limma" = "#B96461",
  "glmGamPoi (cell)" = "#EAB578",
  "Nebula" =  "#E4A6A7",
  "Devil (base)" = "#099668",
  "Devil (mixed)" = "#099668",
  "limma (Pb)" = "#B96461"
)


all_null_plots <- function(author, is.pb, algos = c("glmGamPoi (Pb)", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)"),
                                                    ct.indexes = NULL, genes.values=NULL, n_samples_vec=NULL, pde.values=NULL) {
  if (is.pb) {
    head_foler = "nullpower/null_subject"
  } else {
    head_foler = "nullpower/null_cell"
  }

  fl <- list.files(head_foler, full.names = TRUE)[grepl(author, list.files(head_foler, full.names = TRUE))]

  if (is.null(n_samples_vec)) {n_samples_vec <- str_extract_all(fl, "(?<=n\\.)\\d+(?=\\.ngenes)") %>% unlist() %>% unique()}
  n_samples <- 4
  plots <- lapply(n_samples_vec, function(n_samples) {
    print(n_samples)
    fl_samples <- fl[grepl(paste0("n.", n_samples), fl)]

    if (is.null(ct.indexes)) { ct.indexes <- str_extract_all(fl_samples, "(?<=ct.)\\d+(?=\\.prob)") %>% unlist() %>% unique() }
    if (is.null(genes.values)) { genes.values <- str_extract_all(fl_samples, "(?<=ngenes.)\\d+(?=\\.ct)") %>% unlist() %>% unique() }
    iter.values <- str_extract_all(fl_samples, "(?<=iter.)\\d+(?=\\.cs)") %>% unlist() %>% unique()
    if (is.null(pde.values)) {
      pde.values <- lapply(fl_samples, function(l) {
        unlist(strsplit(unlist(strsplit(l, "probde."))[2], ".iter"))[1]
      }) %>% unlist() %>% unique()
    }


    dd <- lapply(ct.indexes, function(ct.index) {
      dd_inner <- lapply(pde.values, function(pde) {
        dd_inner_inner <- lapply(iter.values, function(i.iter) {
          n_genes <- as.integer(as.numeric(pde) * 1000)

          file_name <- paste0(head_foler, "/", author ,".n.", n_samples,'.ngenes.',n_genes,".ct.",ct.index,".probde.",pde,".iter.",i.iter,".csv")
          if (file.exists(file_name)) {
            d <- read.delim(file_name, sep = ",")
            colnames(d) <- c("X", "glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)")

            mask <- colnames(d) %in% algos
            mask[1] <- TRUE
            d <- d[,mask]

            cols <- colnames(d)
            d <- lapply(2:ncol(d), function(c) {
              values = d[,c] %>% sort()
              x = seq(0,1,length = length(values))
              dplyr::tibble(x = x, observed_p_value = values, name = colnames(d)[c])
            }) %>% do.call("bind_rows", .) %>% dplyr::mutate(ct.index = ct.index, n.genes = n_genes)
            return(d)
          }
        }) %>% do.call("bind_rows", .)
        dd_inner_inner
      }) %>% do.call('bind_rows', .)
    }) %>% do.call("bind_rows", .)

    colnames(dd)

    dd <- dd %>%
      #group_by(name, ct.index, n.genes) %>%
      group_by(name, n.genes) %>%
      dplyr::arrange(observed_p_value) %>%
      dplyr::mutate(x = row_number() / n())
    xy_lines <- dd %>%
      group_by(ct.index, n.genes) %>%
      dplyr::summarise(min = 0, max = 1) %>%
      pivot_longer(!c(ct.index, n.genes))

    p <- dd %>%
      #dplyr::filter(grepl("devil", name) | grepl("Nebula", name)) %>%
      ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
      #geom_line() +
      geom_point(size = .1) +
      geom_line(data = xy_lines, mapping = aes(x=value, y=value), col="black", linetype="dashed") +
      #geom_abline(slope = 1, intercept = 0, col = "black", linetype="dashed") +
      #ggtitle(paste0("Cell type ", ct.index, " - ", n_samples, " patients"), subtitle = paste0("Patient hierarchy ", is.pb)) +
      theme_bw() +
      labs(x = "Uniform quantiles", y="Observed p-value", col="Algorithm") +
      #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
      scale_color_manual(values = method_colors) +
      #facet_wrap(~paste0(n.genes, " genes")) +
      ggh4x::facet_nested(~paste0(n_samples, " patient")+"N genes"+n.genes) +
      theme(legend.position = "bottom") +
      scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n=3))

    p #+ ggtitle(paste("N patients = ", n_samples))
  })

  plots
}


all_pow_plots <- function(author, is.pb, algos = c("glmGamPoi (Pb)", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)"),
                          ct.indexes = NULL, genes.values=NULL, n_samples_vec=NULL, pde.values=NULL) {
  if (is.pb) {
    head_foler = "nullpower/pow_subject/"
  } else {
    head_foler = "nullpower/pow_cell/"
  }

  fl <- list.files(head_foler, full.names = TRUE)[grepl(author, list.files(head_foler, full.names = TRUE))]

  if (is.null(n_samples_vec)) {n_samples_vec <- str_extract_all(fl, "(?<=n\\.)\\d+(?=\\.ngenes)") %>% unlist() %>% unique()}

  plots <- lapply(n_samples_vec, function(n_samples) {
    fl_samples <- fl[grepl(n_samples, fl)]

    # if (is.null(ct.indexes)) { ct.indexes <- str_extract_all(fl_samples, "(?<=ct.)\\d+(?=\\.fc)") %>% unlist() %>% unique() }
    # if (is.null(genes.values)) { genes.values <- str_extract_all(fl_samples, "(?<=ngenes.)\\d+(?=\\.ct)") %>% unlist() %>% unique() }
    #
    # dd <- lapply(ct.indexes, function(ct.index) {
    #   dd_inner <- lapply(genes.values, function(n_genes) {
    #     if (file.exists(paste0(head_foler, "/", author ,".n.", n_samples,'.ngenes.',n_genes,".ct.",ct.index,".fc.0.5.csv"))) {
    #       #d <- read.delim(paste0(head_foler, "/", author ,".n.", n_samples, ".ct.",ct.index,".fc.0.5.csv"), sep = ",")
    #       d <- read.delim(paste0(head_foler, "/", author ,".n.", n_samples,'.ngenes.',n_genes,".ct.",ct.index,".fc.0.5.csv"), sep = ",")
    #       colnames(d) <- c("X", "glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)")
    #
    #       mask <- colnames(d) %in% algos
    #       mask[1] <- TRUE
    #       d <- d[,mask]
    #
    #       cols <- colnames(d)
    #       d <- lapply(2:ncol(d), function(c) {
    #         values = d[,c] %>% sort(decreasing = TRUE)
    #         values[values <= 1e-300] = 1e-300
    #         x = seq(0,1,length = length(values))
    #         dplyr::tibble(x = x, observed_p_value = -log10(values), name = colnames(d)[c])
    #       }) %>% do.call("bind_rows", .) %>% dplyr::mutate(ct.index = ct.index, n.genes = n_genes)
    #       return(d)
    #     }
    #   }) %>% do.call("bind_rows", .)
    # }) %>% do.call("bind_rows", .)

    if (is.null(ct.indexes)) { ct.indexes <- str_extract_all(fl_samples, "(?<=ct.)\\d+(?=\\.prob)") %>% unlist() %>% unique() }
    if (is.null(genes.values)) { genes.values <- str_extract_all(fl_samples, "(?<=ngenes.)\\d+(?=\\.ct)") %>% unlist() %>% unique() }
    iter.values <- str_extract_all(fl_samples, "(?<=iter.)\\d+(?=\\.cs)") %>% unlist() %>% unique()
    if (is.null(pde.values)) {
      pde.values <- lapply(fl_samples, function(l) {
        unlist(strsplit(unlist(strsplit(l, "probde."))[2], ".iter"))[1]
      }) %>% unlist() %>% unique()
    }

    dd <- lapply(ct.indexes, function(ct.index) {
      dd_inner <- lapply(pde.values, function(pde) {
        dd_inner_inner <- lapply(iter.values, function(i.iter) {
          n_genes <- as.integer(as.numeric(pde) * 1000)
          file_name <- paste0(head_foler, "/", author ,".n.", n_samples,'.ngenes.',n_genes,".ct.",ct.index,".probde.",pde,".iter.",i.iter,".csv")
          if (file.exists(file_name)) {
            d <- read.delim(file_name, sep = ",")
            colnames(d) <- c("X", "glmGamPoi (Pb)", "edgeR", "limma", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)")

            mask <- colnames(d) %in% algos
            mask[1] <- TRUE
            d <- d[,mask]

            cols <- colnames(d)
            d <- lapply(2:ncol(d), function(c) {
              values = d[,c] %>% sort(decreasing = TRUE)
              values[values <= 1e-300] = 1e-300
              x = seq(0,1,length = length(values))
              dplyr::tibble(x = x, observed_p_value = -log10(values), name = colnames(d)[c])
            }) %>% do.call("bind_rows", .) %>% dplyr::mutate(ct.index = ct.index, n.genes = n_genes)
            return(d)
          }
        }) %>% do.call("bind_rows", .)
        dd_inner_inner
      }) %>% do.call('bind_rows', .)
    }) %>% do.call("bind_rows", .)

    colnames(dd)
    dd <- dd %>%
      #group_by(name, ct.index, n.genes) %>%
      group_by(name, n.genes) %>%
      dplyr::arrange(observed_p_value) %>%
      dplyr::mutate(x = row_number() / n())

    p <- dd %>%
      #dplyr::filter(grepl("devil", name) | grepl("Nebula", name)) %>%
      ggplot(mapping = aes(x=x, y=observed_p_value, col=name)) +
      #geom_line() +
      geom_point(size = .1) +
      #theme_minimal() +
      theme_bw() +
      scale_y_continuous(trans = "log10") +
      labs(x = "Uniform quantiles", y="Observed -log10(p-value)", col="Algorithm") +
      theme(legend.position = "bottom") +
      #facet_wrap(~ct.index) +
      #ggh4x::facet_nested("Cell type index"+ct.index~"N genes"+n.genes) +
      ggh4x::facet_nested(~paste0(n_samples, " patient")+"N genes"+n.genes) +
      #scale_color_manual(values = c("steelblue", "yellow", "indianred3", "orange", "purple", "forestgreen", "pink")) +
      scale_color_manual(values = method_colors) +
      scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n=3))

    p #+ ggtitle(paste("N patients = ", n_samples))
  })

  plots
}

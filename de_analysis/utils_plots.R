

custom_labels <- function(x) {
  sapply(x, function(xi) {
    if (xi == floor(xi)) {
      as.character(xi)       # No decimals for integers
    } else {
      format(round(xi, 2), nsmall = 2)  # 2 decimals otherwise
    }
  })
}

method_colors = c(
  #"glmGamPoi (Pb)" = "#A22E29",
  #"edgeR" = "#7D629E",
  "limma" = "#B96461",
  "limma (Pb)" = "#B96461",
  "glmGamPoi (cell)" = "#EAB578",
  "glmGamPoi" = "#EAB578",
  "Nebula" =  'steelblue', #"#B0C4DE",
  "NEBULA" =  'steelblue', #"#B0C4DE",
  "Devil (base)" = "#099668",
  "Devil (mixed)" = "#099668",
  "Devil" = "#099668",
  "devil" = "#099668"
)

method_levels = c("devil", "NEBULA", "glmGamPoi", "limma")


all_null_plots <- function(author, is.pb, algos = c("glmGamPoi (Pb)", "glmGamPoi (cell)", "Nebula", "Devil (base)", "Devil (mixed)"),
                                                    ct.indexes = NULL, genes.values=NULL, n_samples_vec=NULL, pde.values=NULL, only_tibble=FALSE) {
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

    if (only_tibble) {
      return(dd %>% dplyr::mutate(n.samples = n_samples))
    }

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
                          ct.indexes = NULL, genes.values=NULL, n_samples_vec=NULL, pde.values=NULL, only_tibble=FALSE) {
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

    if (only_tibble) {
      return(dd %>% dplyr::mutate(n.samples = n_samples))
    }

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


plot_pvalues = function(author, method_cellwise, method_patientwise) {
  d1 <- all_null_plots(author, FALSE, algos =method_cellwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20), only_tibble=TRUE)[[1]] %>%
    dplyr::mutate(ytype = "p-value", xtype="Cell-wise")
  d2 <- all_null_plots(author, TRUE, algos = method_patientwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20), only_tibble = T)[[1]] %>%
    dplyr::mutate(ytype = "p-value", xtype="Patient-wise")
  d3 <- all_pow_plots(author, FALSE, algos =method_cellwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20), only_tibble = T)[[1]] %>%
    dplyr::mutate(ytype = "-log10 p-value", xtype="Cell-wise")
  d4 <- all_pow_plots(author, TRUE, algos = method_patientwise, ct.indexes = c(1), pde.values = c(.05), n_samples_vec = c(20), only_tibble = T)[[1]] %>%
    dplyr::mutate(ytype = "-log10 p-value", xtype="Patient-wise")


  method_levels <- c("limma", "glmGamPoi", "glmGamPoi (cell)", "Nebula", "NEBULA", "Devil (mixed)", "Devil (base)", "Devil", "devil")

  dplyr::bind_rows(d1, d2) %>%
    dplyr::mutate(name = dplyr::if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = dplyr::if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = dplyr::if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::group_by(xtype, ytype) %>%
    dplyr::mutate(name = factor(name, levels = method_levels)) %>%
    ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
    ggridges::geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = F, alpha = .7) +
    scale_color_manual(values = sort(method_colors)) +
    scale_fill_manual(values = sort(method_colors)) +
    facet_grid(~xtype, scales = "free") +
    theme_bw() +
    labs(x = "p-value", y="", col="Algorithm") +
    scale_color_manual(values = method_colors) +
    #facet_wrap(~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
    theme(legend.position = "none") +
    theme()

  dplyr::bind_rows(d3, d4) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::group_by(xtype, ytype) %>%
    dplyr::mutate(name = factor(name, levels = method_levels)) %>%
    ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
    ggridges::geom_density_ridges(stat = "binline", bins = 100, scale = 0.95, draw_baseline = F, alpha = .7) +
    scale_color_manual(values = sort(method_colors)) +
    scale_fill_manual(values = sort(method_colors)) +
    facet_grid(~xtype, scales = "free") +
    theme_bw() +
    labs(x = bquote(-log[10]~ "(p-value)"), y="", col="Algorithm") +
    scale_color_manual(values = method_colors) +
    #facet_wrap(~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(transform = "log10") +
    #scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1)) +
    scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
    theme(legend.position = "none") +
    theme()


  pB <- dplyr::bind_rows(d1, d2) %>%
    dplyr::mutate(name = dplyr::if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = dplyr::if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = dplyr::if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::group_by(xtype, ytype) %>%
    dplyr::mutate(name = factor(name, levels = method_levels)) %>%
    ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
    ggridges::geom_density_ridges(alpha = .7, scale = 1) +
    #ggridges::geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = F) +
    scale_color_manual(values = sort(method_colors)) +
    scale_fill_manual(values = sort(method_colors)) +
    facet_grid(~xtype, scales = "free") +
    theme_bw() +
    labs(x = "p-value", y="", col="Algorithm") +
    scale_color_manual(values = method_colors) +
    #facet_wrap(~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1)) +
    scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
    theme(legend.position = "none") +
    theme()

  pC <- dplyr::bind_rows(d3, d4) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::group_by(xtype, ytype) %>%
    dplyr::mutate(name = factor(name, levels = method_levels)) %>%
    ggplot(mapping = aes(x=observed_p_value, col=name, fill=name, y=name)) +
    ggridges::geom_density_ridges(alpha = .7, scale = 1) +
    scale_color_manual(values = sort(method_colors)) +
    scale_fill_manual(values = sort(method_colors)) +
    facet_grid(~xtype, scales = "free") +
    theme_bw() +
    labs(x = bquote(-log[10]~ "(p-value)"), y="", col="Algorithm") +
    scale_color_manual(values = method_colors) +
    #facet_wrap(~paste0(n.genes, " genes")) +
    theme(legend.position = "bottom") +
    scale_x_continuous(transform = "log10") +
    #scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0,1)) +
    scale_y_discrete(expand = expand_scale(mult = c(0.01, .25))) +
    theme(legend.position = "none") +
    theme()

  list(null_pvalue = pB, de_pvalue = pC)
}

plot_MCCs = function(author, method_cellwise, method_patientwise) {
  res <- readRDS("nullpower/final_res/results.rds")
  a <- author

  pfalse <- res %>%
    na.omit() %>%
    dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::filter(author == a) %>%
    dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = patients, Ngenes = paste0(ngenes, " genes")) %>%
    #dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
    dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
    dplyr::group_by(Npatients, ngenes) %>%
    ggplot(mapping = aes(x=as.factor(ngenes), y=MCC, col=name)) +
    geom_boxplot() +
    scale_color_manual(values = method_colors) +
    labs(x = "Number of DE genes", y = "MCC", col="Model", linetype = "N patients", shape = "N patients") +
    ggh4x::facet_nested(~"Number of samples"+Npatients, scales = "free_y") +
    theme_bw() +
    theme(text = element_text(size = 12), legend.position = "bottom")

  ptrue <- res %>%
    na.omit() %>%
    dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::filter(author == a) %>%
    dplyr::mutate(Dataset = paste0("Dataset : ", author), Npatients = patients, Ngenes = paste0(ngenes, " genes")) %>%
    #dplyr::mutate(Npatients = factor(Npatients, levels = c("4 patients", "20 patients"))) %>%
    dplyr::mutate(Ngenes = factor(Ngenes, levels = c("5 genes", "25 genes", "50 genes"))) %>%
    dplyr::group_by(Npatients, ngenes) %>%
    ggplot(mapping = aes(x=as.factor(ngenes), y=MCC, col=name)) +
    geom_boxplot() +
    scale_color_manual(values = method_colors) +
    labs(x = "Number of DE genes", y = "MCC", col="Model", linetype = "N patients", shape = "N patients") +
    ggh4x::facet_nested(~"Number of samples"+Npatients, scales = "free_y") +
    theme_bw() +
    theme(text = element_text(size = 12), legend.position = "bottom")

  list(cellwise = pfalse, patientwise = ptrue)
}

plot_ks = function(author, method_cellwise, method_patientwise) {
  res <- readRDS("nullpower/final_res/results.rds")
  a = author

  {
    dx <- .15

    r_ks <- res %>%
      na.omit() %>%
      dplyr::filter(author == a) %>%
      dplyr::filter(is.pb == FALSE, name %in% method_cellwise)

    r_ks_tot <- lapply(unique(r_ks$name), function(n) {
      MCCs <- r_ks %>%
        dplyr::filter(name==n) %>%
        dplyr::pull(MCC)
      ECDFs <- lapply(MCCs, function(mcc) {
        sum(MCCs <= mcc) / length(MCCs)
      }) %>% unlist()

      dplyr::tibble(MCCs =MCCs, ECDFs=ECDFs, name=n)
    }) %>% do.call("rbind", .)

    ks_pvals <- lapply(method_cellwise[method_cellwise!="Devil (base)"], function(m) {
      print(m)
      pval = ks.test(filter(r_ks_tot, name=="Devil (base)")$MCCs, filter(r_ks_tot, name==m)$MCCs)$p.value
      if (pval == 0) {
        pval <- "< 2e-16"
      } else {
        pval <- paste0("= ", round(pval, 2))
      }
      dplyr::tibble(m=m, pval=pval)
    }) %>% do.call("bind_rows", .)

    ks_false <- r_ks_tot %>%
      dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
      dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
      dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
      dplyr::mutate(name = factor(name, levels = method_levels)) %>%
      ggplot(mapping = aes(x=MCCs, y=ECDFs, col=name)) +
      geom_line(linewidth = .8, position = position_dodge(width = .02)) +
      scale_color_manual(values = sort(method_colors)) +
      #ggtitle(paste0("Cell-wise")) +
      labs(x = "MCC", y = "Empirical CDF", col="Model") +
      theme_bw()

    L = ggplot_build(ks_false)$layout$panel_params[[1]]
    Lx = (abs(L$x.range[2] - L$x.range[1]) * .2) + L$x.range[1]
    Ly = (abs(L$y.range[2] - L$y.range[1]) * .9) + L$y.range[1]

    ks_false <- ks_false +
      annotate(geom='label', x=Lx, y=Ly, label=paste0('p ', ks_pvals$pval[1]), color=method_colors[ks_pvals$m[1]]) +
      annotate(geom='label', x=Lx, y=Ly - dx, label=paste0('p ', ks_pvals$pval[2]), color=method_colors[ks_pvals$m[2]]) +
      annotate(geom='label', x=Lx, y=Ly - 2*dx, label=paste0('p ', ks_pvals$pval[3]), color=method_colors[ks_pvals$m[3]]) +
      theme(legend.position = "left")
  }

  # kolmogorov smirnof plots - TRUE ####
  {
    r_ks <- res %>%
      na.omit() %>%
      dplyr::filter(author == a) %>%
      dplyr::filter(is.pb == TRUE, name %in% method_patientwise)

    r_ks_tot <- lapply(unique(r_ks$name), function(n) {
      MCCs <- r_ks %>%
        dplyr::filter(name==n) %>%
        dplyr::pull(MCC)
      ECDFs <- lapply(MCCs, function(mcc) {
        sum(MCCs <= mcc) / length(MCCs)
      }) %>% unlist()

      dplyr::tibble(MCCs =MCCs, ECDFs=ECDFs, name=n)
    }) %>% do.call("rbind", .)

    ks_pvals <- lapply(method_patientwise[method_patientwise!="Devil (mixed)"], function(m) {
      print(m)
      pval = ks.test(filter(r_ks_tot, name=="Devil (mixed)")$MCCs, filter(r_ks_tot, name==m)$MCCs)$p.value
      true_pval <- pval
      if (pval <= 1e-6) {
        pval <- "< 1e-6"
      } else {
        pval <- paste0("= ", round(pval, 3))
      }
      dplyr::tibble(m=m, pval=pval, true_pval=true_pval)
    }) %>% do.call("bind_rows", .)

    ks_true <- r_ks_tot %>%
      dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
      dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
      dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
      dplyr::mutate(name = factor(name, levels = method_levels)) %>%
      ggplot(mapping = aes(x=MCCs, y=ECDFs, col=name)) +
      geom_line(linewidth = .8, position = position_dodge(width = .02)) +
      scale_color_manual(values = method_colors) +
      #ggtitle(paste0("Patient-wise")) +
      #ggtitle(paste0("Cell-wise")) +
      labs(x = "MCC", y = "Empirical CDF", col="Model") +
      theme_bw()

    L = ggplot_build(ks_true)$layout$panel_params[[1]]
    Lx = (abs(L$x.range[2] - L$x.range[1]) * .2) + L$x.range[1]
    Ly = (abs(L$y.range[2] - L$y.range[1]) * .9) + L$y.range[1]

    ks_true <- ks_true +
      annotate(geom='label', x=Lx, y=Ly, label=paste0('p ', ks_pvals$pval[1]), color=method_colors[ks_pvals$m[1]]) +
      annotate(geom='label', x=Lx, y=Ly - dx, label=paste0('p ', ks_pvals$pval[2]), color=method_colors[ks_pvals$m[2]]) +
      annotate(geom='label', x=Lx, y=Ly - 2*dx, label=paste0('p ', ks_pvals$pval[3]), color=method_colors[ks_pvals$m[3]]) +
      theme(legend.position = "left")
  }

  list(cellwise = ks_false, patientwise = ks_true)
}

plot_timing = function(author, competitor = "NEBULA") {
  a = author
  timing_plot <- readRDS(paste0("nullpower/timing_results/", a,".rds")) %>%
    dplyr::filter(algo %in% c("Devil (base)", "glmGamPoi (cell)", "Nebula")) %>%
    dplyr::rename(name = algo) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
    dplyr::rename(algo = name) %>%
    dplyr::group_by(author, is.pb, n.sample, n.gene, int.ct, iter) %>%
    dplyr::mutate(time_fold = timings[algo == competitor] / timings) %>%
    dplyr::mutate(cell_order = ifelse(n.cells < 1000, "< 1k", if_else(n.cells > 20000, "> 20k", "1k-20k"))) %>%
    dplyr::mutate(cell_order = factor(cell_order, levels = c("< 1k", "1k-20k", "> 20k"))) %>%
    #dplyr::filter(algo %in% c("glmGamPoi", "NEBULA")) %>%
    #ggplot(mapping = aes(x=cell_order, y=timings, col=algo)) +
    ggplot(mapping = aes(x=cell_order, y=time_fold, col=algo)) +
    geom_boxplot() +
    scale_color_manual(values = method_colors) +
    labs(x = "Number of cells", y=paste0("Speedup (vs. ",competitor,")"), col = "Model") +
    theme_bw() +
    geom_hline(yintercept = 1, color = "darkslategray", linetype = 'dashed') +
    coord_flip() +
    scale_y_continuous(transform = "log10")
  timing_plot
}


plot_all_models = function() {
  method_cellwise <- c("glmGamPoi (cell)", "Devil (base)", "limma", "Nebula", "glmGamPoi (fixed)", "edgeR", "edgeR (Pb)", "limma (Pb)")
  method_patientwise <- c("Nebula", "Devil (mixed)", "limma", "glmGamPoi (cell)", "glmGamPoi (fixed)", "edgeR", "edgeR (Pb)", "limma (Pb)")
  method_colors = c(
    "glmGamPoi (fixed)" = "#A22E29",
    "edgeR" = "#7D629E",
    "edgeR (Pb)" = "#7D629E",
    "limma" = "#B96461",
    "limma (Pb)" = "#B96461",
    "glmGamPoi (cell)" = "#EAB578",
    "glmGamPoi" = "#EAB578",
    "Nebula" =  'steelblue', #"#B0C4DE",
    "NEBULA" =  'steelblue', #"#B0C4DE",
    "Devil (base)" = "#099668",
    "Devil (mixed)" = "#099668",
    "Devil" = "#099668",
    "devil" = "#099668"
  )

  res = readRDS("nullpower/final_res/results.rds") %>%
    dplyr::filter((is.pb == TRUE & name %in% method_patientwise) | (is.pb == FALSE & name %in% method_cellwise))

  all_timing = lapply(list.files("nullpower/timing_results/"), function(p) {
    readRDS(file.path("nullpower/timing_results/", p)) %>% dplyr::mutate(dataset = p)
  }) %>% do.call("bind_rows", .) %>%
    dplyr::filter((is.pb == TRUE & algo %in% method_patientwise) | (is.pb == FALSE & algo %in% method_cellwise))

  all_MCC_boxplots = res %>%
    dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
    dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
    dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
    ggplot(mapping = aes(x=name, y=MCC, col=name)) +
    geom_boxplot() +
    ggh4x::facet_nested("Dataset"+author~is.pb) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    labs(color = "Model") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )

  all_time_boxplots = all_timing %>%
    dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
    dplyr::mutate(name = algo) %>%
    dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
    dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
    ggplot(mapping = aes(x=name, y=timings, col=name)) +
    geom_boxplot() +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    labs(color = "Model") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_y_continuous(transform = "log10") +
    coord_flip() +
    facet_wrap(~author) +
    labs(y = "Runtime (seconds)")

  failure_rate_plot = all_timing %>%
    dplyr::filter(is.pb) %>%
    dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
    dplyr::mutate(name = algo) %>%
    dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
    dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
    dplyr::mutate(is_bad = is.na(timings)) %>%
    dplyr::group_by(name, is.pb, author) %>%
    dplyr::summarise(`Failure rate` = sum(is_bad) / n()) %>%
    ggplot(mapping = aes(x=author, y=`Failure rate`, fill=name, col=name)) +
    geom_col(position = "dodge") +
    ylim(c(0,1)) +
    ggh4x::facet_nested(~is.pb) +
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_colors) +
    labs(color = "", fill="", x = "Dataset name") +
    theme_bw()

  # failure_rate_plot = all_timing %>%
  #   dplyr::filter(is.pb) %>%
  #   dplyr::mutate(is.pb = if_else(is.pb, "Patient-wise", "Cell-wise")) %>%
  #   dplyr::mutate(name = algo) %>%
  #   dplyr::mutate(name = ifelse(grepl("Devil", name, fixed = TRUE), "devil", name)) %>%
  #   dplyr::mutate(name = ifelse(grepl("(cell)", name, fixed = TRUE), "glmGamPoi", name)) %>%
  #   dplyr::mutate(is_bad = is.na(timings)) %>%
  #   dplyr::group_by(name, is.pb, author) %>%
  #   dplyr::summarise(`Failure rate` = sum(is_bad) / n()) %>%
  #   ggplot(mapping = aes(x=name, y=`Failure rate`, fill=name, col=name)) +
  #   geom_col(col="black") +
  #   ylim(c(0,1)) +
  #   ggh4x::facet_nested("Dataset"+author~is.pb) +
  #   scale_fill_manual(values = method_colors) +
  #   labs(color = "") +
  #   theme_bw() +
  #   coord_flip() +
  #   theme(
  #     legend.position = "none",
  #     axis.title.y = element_blank()
  #   )

  list(MCC = all_MCC_boxplots, timing = all_time_boxplots, failure_rate = failure_rate_plot)
}

plot_MCCs_boxplots = function(a = NULL) {
  res = readRDS("nullpower/final_res/results.rds")
  if (!is.null(a)) {
    res = res %>% dplyr::filter(author == a)
  }

  df_cellwise = res %>%
    na.omit() %>%
    dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  df_patientwise = res %>%
    na.omit() %>%
    dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
    dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
    dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
    dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name))

  df_all = dplyr::bind_rows(df_patientwise, df_cellwise)

  df_all %>%
    dplyr::mutate(name = factor(name, levels = method_levels)) %>%
    dplyr::group_by(name, ct.index, is.pb, author, patients) %>%
    dplyr::summarise(MCC = median(MCC)) %>%
    dplyr::mutate(is.pb = ifelse(is.pb, "Patient-wise", "Cell-wise")) %>%
    ggplot(mapping = aes(x = name, y=MCC, col=name)) +
    geom_boxplot() +
    geom_point() +
    coord_flip() +
    ggh4x::facet_nested(~is.pb+paste0(patients, " patients")) +
    scale_y_continuous(labels = custom_labels) +
    scale_color_manual(values = method_colors) +
    theme_bw() +
    labs(x = "", col = "")
}

#
#
# df_cellwise %>%
#   dplyr::group_by(name, is.pb, author) %>%
#   dplyr::arrange(MCC) %>%
#   dplyr::mutate(idx = row_number()) %>%
#   ggplot(mapping = aes(x = idx, y=MCC, col=name)) +
#   #geom_point() +
#   geom_line() +
#   scale_color_manual(values = method_colors) +
#   facet_wrap(~author) +
#   theme_bw()
#
#
# model_name = df_cellwise$name[1]
# d = lapply(unique(df_cellwise$name), function(model_name) {
#   MCCs = df_cellwise %>% dplyr::filter(name == model_name)  %>% dplyr::pull(MCC)
#   lapply(MCC_cuts, function(cut) {
#     n = as.integer(f * length(MCCs))
#     pobserved = lapply(1:1000, function(i) {
#       sum(sample(MCCs, n) <= cut) / n
#     }) %>% unlist()
#     dplyr::tibble(name = model_name,
#                   cut = cut,
#                   y = median(pobserved),
#                   ylow = stats::quantile(pobserved, .05),
#                   yhigh = stats::quantile(pobserved, .95))
#   }) %>% do.call("bind_rows", .)
# }) %>% do.call("bind_rows", .)
#
# d %>%
#   ggplot(mapping = aes(x = cut, y=y, ymin=ylow, ymax=yhigh, fill=name, col=name)) +
#   geom_line() +
#   geom_ribbon(alpha = .5, linewidth = 0) +
#   scale_color_manual(values = method_colors) +
#   scale_fill_manual(values = method_colors) +
#   theme_bw()
#
#
#
# df_patientwise %>%
#   dplyr::group_by(name, is.pb, author, patients) %>%
#   dplyr::arrange(MCC) %>%
#   dplyr::summarise(n = n()) %>%
#   dplyr::mutate(idx = row_number()) %>%
#   ggplot(mapping = aes(x = idx, y=MCC, col=name)) +
#   #geom_point() +
#   geom_line() +
#   scale_color_manual(values = method_colors) +
#   facet_grid(patients~author) +
#   theme_bw()
#
# model_name = df_patientwise$name[1]
# d = lapply(unique(df_patientwise$name), function(model_name) {
#   MCCs = df_patientwise %>% dplyr::filter(name == model_name)  %>% dplyr::pull(MCC)
#   lapply(MCC_cuts, function(cut) {
#     n = as.integer(f * length(MCCs))
#     pobserved = lapply(1:1000, function(i) {
#       sum(sample(MCCs, n) <= cut) / n
#     }) %>% unlist()
#     dplyr::tibble(name = model_name,
#                   cut = cut,
#                   y = median(pobserved),
#                   ylow = stats::quantile(pobserved, .05),
#                   yhigh = stats::quantile(pobserved, .95))
#   }) %>% do.call("bind_rows", .)
# }) %>% do.call("bind_rows", .)
#
# d %>%
#   ggplot(mapping = aes(x = cut, y=y, ymin=ylow, ymax=yhigh, fill=name, col=name)) +
#   geom_line(linewidth = 1) +
#   geom_ribbon(alpha = .25, linewidth = 0)
#
#
#
# MCC_cuts = seq(0, 1, by = .05)
# f = .5
#
#
#
#
# author = "hsc"
# res <- readRDS("nullpower/final_res/results.rds")
# a <- author
#
# df_cellwise = res %>%
#   na.omit() %>%
#   dplyr::filter(is.pb == FALSE, name %in% method_cellwise) %>%
#   dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
#   dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
#   dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
#   dplyr::filter(author == a) %>%
#   dplyr::group_by(name, is.pb, author, patients) %>%
#   dplyr::summarise(y = mean(MCC))
#
# n = df_cellwise$name[1]
# lapply(unique(df_cellwise$name), function(n) {
#   MCCs = df_cellwise %>% dplyr::filter(name == n)  %>% dplyr::pull(MCC)
# })
#
# MCC_cuts = seq(0, 1, by = .05)
# f = .5
#
#
#
#
#
#
# df_patientwise <- res %>%
#   na.omit() %>%
#   dplyr::filter(is.pb == TRUE, name %in% method_patientwise) %>%
#   dplyr::mutate(name = if_else(grepl("Devil", name), "devil", name)) %>%
#   dplyr::mutate(name = if_else(grepl("glmGamPoi", name), "glmGamPoi", name)) %>%
#   dplyr::mutate(name = if_else(grepl("Nebula", name), "NEBULA", name)) %>%
#   dplyr::filter(author == a) %>%
#   dplyr::group_by(name, is.pb, author, patients) %>%
#   dplyr::summarise(y = mean(MCC))
#
# df = dplyr::bind_rows(df_cellwise, df_patientwise)
# df %>%
#   dplyr::mutate(y_name = paste0(is.pb, "_", patients)) %>%
#   dplyr::group_by(y_name) %>%
#   dplyr::mutate(yf = y / max(y)) %>%
#   ggplot(mapping = aes(x=name, y=y_name, fill=yf)) +
#   geom_raster()
#
#
#
#
#
#
# df_cellwise %>%
#   dplyr::group_by(name, is.pb, author, patients, ct.index) %>%
#   dplyr::summarise(y = median(MCC)) %>%
#   ggplot(mapping = aes(x = ct.index, y=y, col=name, linetype = as.factor(patients))) +
#   geom_point() +
#   geom_line() +
#   scale_color_manual(values = method_colors) +
#   facet_wrap(~author) +
#   theme_bw()
#
# df_patientwise %>%
#   dplyr::group_by(name, is.pb, author, patients, ct.index) %>%
#   dplyr::summarise(y = mean(MCC)) %>%
#   ggplot(mapping = aes(x = ct.index, y=y, col=name, linetype = as.factor(patients))) +
#   geom_point() +
#   geom_line() +
#   scale_color_manual(values = method_colors) +
#   facet_grid(patients~author) +
#   theme_bw()
#
#

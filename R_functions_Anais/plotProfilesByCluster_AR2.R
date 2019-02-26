plotProfilesByCluster_AR <-function (riskProfObj, whichCovariates = NULL, rhoMinimum = NULL,
                                     useProfileStar = TRUE, rho_in_xlabels = FALSE, covariate_info = list(title = NULL,
                                                                                                          levels = NULL, labels = NULL, split = NULL))
{
  profileDF <- tabulateCovariateProfiles_AR(riskProfObj = riskProfObj,
                                            whichCovariates = whichCovariates, rhoMinimum = rhoMinimum,
                                            useProfileStar = useProfileStar)
  empirical_proportions <- profileDF %>% dplyr::group_by(category,
                                                         covname) %>% dplyr::summarise(x = mean(mean)) %>% dplyr::group_by(covname) %>%
    dplyr::arrange(category) %>% dplyr::mutate(emp_propn = cumsum(x)) %>%
    dplyr::filter(emp_propn < max(emp_propn)) %>% dplyr::ungroup() %>%
    dplyr::select(-x)
  covtab <- profileDF %>% dplyr::left_join(empirical_proportions,
                                           by = c("covname", "category")) %>% dplyr::group_by(cluster,
                                                                                              category, covname, fillColor, rhoMean, rhoRank, emp_propn) %>%
    dplyr::summarise(prop = mean(est)) %>% dplyr::mutate(covlab = ifelse(rho_in_xlabels,
                                                                         sprintf("%s (%.2f)", covname, rhoMean), covname)) %>%
    dplyr::ungroup()
  if (!is.null(covariate_info$levels) & !is.null(covariate_info$labels)) {
    covtab <- covtab %>% dplyr::mutate(category = factor(category,
                                                         levels = covariate_info$levels, labels = covariate_info$labels))
  }
  if (!is.null(covariate_info$split)) {
    if (length(covariate_info$split) < 3) {
      covariate_info$split <- rep(covariate_info$split,
                                  3 - length(covariate_info$split))
      covariate_info$split <- c(covariate_info$split, sprintf("not %s",
                                                              covariate_info$split[2]))
    }
    covtab <- covtab %>% dplyr::mutate(type = stringr::str_detect(covname,
                                                                  covariate_info$split[1]), type = dplyr::recode(as.numeric(type),
                                                                                                                 `1` = covariate_info$split[2], `0` = covariate_info$split[3]))
    facetting_layer <- list(ggplot2::facet_grid(cluster ~
                                                  type, scales = "free_x", space = "free_x"))
  }
  else {
    facetting_layer <- list(ggplot2::facet_grid(cluster ~
                                                  .))
  }
  rhoMinimum <- ifelse(is.null(rhoMinimum), min(covtab$rhoMean),
                       rhoMinimum)
  covariate_info$title <- ifelse(is.null(covariate_info$title),
                                 "Covariate\ncategory", covariate_info$title)
  xaxistitle = ifelse(rho_in_xlabels, sprintf("Covariate (rho >= %.2f, top %i)",
                                              rhoMinimum, length(unique(covtab$covlab))), sprintf("Covariate (top %i)",
                                                                                                  length(unique(covtab$covlab))))
  ggplot2::ggplot(covtab, ggplot2::aes(x = reorder(covlab,
                                                   rhoRank), y = prop, fill = factor(category), alpha = fillColor ==
                                         "high")) + ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::geom_point(aes(y = emp_propn, group = category),
                        col = "black", fill = "white", alpha = 1, shape = 18) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1, vjust = 0.5)) + ggplot2::labs(x = xaxistitle,
                                                                                                y = "Proportion (by cluster)", title = "Covariate profiles",
                                                                                                subtitle = "Covariate categories are plotted with dark fill if they are more prevalent within the cluster than overall.\nOverall (empirical) proportions in each covariate category are indicated by black diamonds.") +
    ggplot2::scale_fill_discrete(name = covariate_info$title) +
    ggplot2::scale_alpha_discrete(guide = FALSE, range = c(0.25,
                                                           1)) + facetting_layer
}

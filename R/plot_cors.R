make_cortb <- function(mat_cor) {
  mat_cor %>% 
    as.data.frame() %>% tibble::rownames_to_column("feature1") %>% 
    tidyr::pivot_longer(cols = -feature1,
                        names_to = "feature2",
                        values_to = "cor") %>% 
    dplyr::mutate(feature1 = forcats::as_factor(feature1),
                  feature2 = factor(feature2, levels = levels(feature1))) %>% 
    dplyr::filter(as.numeric(feature1) <= as.numeric(feature2))
}

make_cortb2 <- function(cor1, cor2, labels = c("cor1", "cor2")) {
  list(cor1, cor2) %>% 
    purrr::map2_dfr(
      labels,
      ~ .x %>% 
        as.data.frame() %>% tibble::rownames_to_column("feature1") %>% 
        tidyr::pivot_longer(cols = -feature1,
                            names_to = "feature2",
                            values_to = "cor") %>% 
        dplyr::mutate(feature1 = forcats::as_factor(feature1),
                      feature2 = factor(feature2, levels = levels(feature1))) %>% 
        dplyr::filter(as.numeric(feature1) < as.numeric(feature2)) %>% 
        dplyr::mutate(label = .y)
    ) %>% 
    dplyr::mutate(label = factor(label, levels = labels))
}
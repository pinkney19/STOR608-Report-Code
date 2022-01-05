my_functional_boxplot <- function (dts, depth_method = c("mbd", "tvd", "extremal", 
                                "dirout", "linfinity", "bd", "erld", 
                                "dq"), depth_values = NULL, emp_factor = 1.5, central_region = 0.5, 
          erld_type = NULL, dq_quantiles = NULL) 
{
  dm <- dim(dts)
  n <- dm[1]
  p <- dm[2]
  if (is.data.frame(dts)) {
    dt <- as.matrix(dts)
  }
  if (any(!is.finite(dts))) {
    stop("Missing or infinite values are not allowed in argument \"dts\"")
  }
  if (!is.array(dts) || !is.numeric(dts)) 
    stop("Argument \"dts\" must be a numeric matrix or dataframe.")
  if (length(dm) != 2) 
    stop("Dimension of 'dts' must be of length 2. Only univariate functional data is supported.")
  if (is.null(depth_values)) {
    depth_method <- match.arg(depth_method)
    if (depth_method == "mbd") {
      depth_values <- modified_band_depth(dts)
    }
    else if (depth_method == "tvd") {
      depth_values <- total_variation_depth(dts)$tvd
    }
    else if (depth_method == "extremal") {
      depth_values <- extremal_depth(dts)
    }
    else if (depth_method == "dirout") {
      depth_values <- -dir_out(dts, return_distance = T)$distance
    }
    else if (depth_method == "linfinity") {
      depth_values <- linfinity_depth(dts)
    }
    else if (depth_method == "bd") {
      depth_values <- band_depth(dts)
    }
    else if (depth_method == "erld") {
      if (is.null(erld_type)) {
        warning("The 'type' argument for extreme rank length depth not specified. Using the default type of 'two_sided'. ")
        depth_values <- extreme_rank_length(dts)
      }
      else {
        depth_values <- extreme_rank_length(dts, type = erld_type)
      }
    }
    else if (depth_method == "dq") {
      if (is.null(dq_quantiles)) {
        warning("Using the default quantile probabilites of 0.025 and 0.975 for directional quantile.")
        depth_values <- -directional_quantile(dts)
      }
      else {
        depth_values <- -directional_quantile(dts, quantiles = dq_quantiles)
      }
    }
  }
  else {
    if (length(depth_values) != n) {
      stop("Length of argument 'depth_values' must be equal to the number of rows in 'dts'.")
    }
  }
  if (central_region >= 1 || central_region <= 0) {
    stop("Argument 'central_region' must be greater than 0 and less than 1.")
  }
  sorted_depths <- sort(depth_values, decreasing = T, index.r = T)
  index_sorted_depth <- sorted_depths$ix
  sorted_depths <- sorted_depths$x
  median_curve <- index_sorted_depth[1]
  n_obs_central <- ceiling(n * central_region)
  center <- dts[index_sorted_depth[1:n_obs_central], ]
  inf <- apply(center, 2, min)
  sup <- apply(center, 2, max)
  distt <- emp_factor * (sup - inf)
  upper <- sup + distt
  lower <- inf - distt
  tdts <- t(dts)
  outlier_test <- (tdts < lower) + (tdts > upper)
  outliers <- which(colSums(outlier_test) > 0)
  if( length(outliers) ){
    upper <- apply(dts[-outliers,],2,max)
    lower <- apply(dts[-outliers,],2,min)
  } else {
    upper <- apply(dts,2,max)
    lower <- apply(dts,2,min)
  }
  return(list(outliers = unname(outliers), depth_values = depth_values, 
              median_curve = median_curve, lower=lower,inf=inf,sup=sup, upper=upper ))
}
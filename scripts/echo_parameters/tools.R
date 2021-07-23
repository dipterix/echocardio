findpeaks <- function(s, min_val = NA, min_distance = 4, min_width = 2){
  df1 <- c(1, diff(s))
  df2 <- c(1, 1, diff(s, differences = 2))

  # necessary condition to be peaks (1st & 2nd derivatives)
  idx <- (df1 * c(df1[-1], 0)) <= 0 & c(df2[-1], 0) < 0

  # Get peaks that are beyond min_val
  if( is.na(min_val) ){
    # automatically decide min_val from standard error
    min_val <- sd(s)
    if(all(s >= 0)){
      min_val <- min(s) + min_val
    } else {
      min_val <- median(s) + min_val
    }
  }
  idx <- which(idx & (s > min_val))

  ord <- order(s[idx], decreasing = TRUE)
  idx_desc <- idx[ord]

  # merge peaks that are less than min_distance away
  for(ii in seq_along(idx_desc)){
    elem <- idx_desc[[ii]]
    idx_desc[idx_desc > 0 & abs(idx_desc - elem) < min_distance] <- elem
  }
  idx_desc <- unique(idx_desc)

  # check +- min_width to see if the it is still a peak
  p_idx <- idx_desc + min_width
  p_idx[p_idx > length(s)] <- length(s)
  m_idx <- idx_desc - min_width
  m_idx[m_idx < 1] <- 1
  sign <- (s[p_idx] - s[idx_desc]) * (s[idx_desc] - s[m_idx])
  idx_desc <- idx_desc[sign <= 0]
  list(
    index = idx_desc,
    values = s[idx_desc]
  )
}


preview_array <- function(z, x, y, ...){
  if(missing(y)){
    y <- seq_len(nrow(z))
  }
  if(missing(x)){
    x <- seq_len(ncol(z))
  }
  if(y[2] < y[1]){
    y <- rev(y)
  }
  image(t(z[rev(seq_len(nrow(arr))), ]),
        x = x,
        y = y, ...)

}

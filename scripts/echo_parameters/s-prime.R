# install.packages('oro.dicom')
library(oro.dicom)
source('./scripts/echo_parameters/tools.R')
root_dir <- 'data/'

scode <- "Pt57"

subject_dcm <- function(name){
  # A4C_RV_S
  normalizePath(file.path(
    root_dir, scode, 'deid-dicom', sprintf('%s_%s.dcm', scode, name)),
    mustWork = TRUE)
}

# ravepy::configure()
# ravepy::add_packages(c('deid', 'pydicom', 'numpy', 'matplotlib'))
ravepy::ensure_ravepy()
# reticulate::repl_python()

pydicom <- reticulate::import("pydicom")
dcm <- pydicom$dcmread(subject_dcm('A4C_RV_S'))
img <- dcm$pixel_array
shape <- dim(img)

# get sub-region
len <- reticulate::py_len(dcm$SequenceOfUltrasoundRegions)
tmp <- dcm$SequenceOfUltrasoundRegions[[1]]
# (0018, 6012) Region Spatial Format               US: 3
# (0018, 6014) Region Data Type                    US: 3
# (0018, 6016) Region Flags                        UL: 2
# (0018, 6018) Region Location Min X0              UL: 50
# (0018, 601a) Region Location Min Y0              UL: 262
# (0018, 601c) Region Location Max X1              UL: 913
# (0018, 601e) Region Location Max Y1              UL: 721
# (0018, 6020) Reference Pixel X0                  SL: 864
# (0018, 6022) Reference Pixel Y0                  SL: 181
# (0018, 6024) Physical Units X Direction          US: 4
# (0018, 6026) Physical Units Y Direction          US: 7
# (0018, 6028) Reference Pixel Physical Value X    FD: 2.9448220001726026
# (0018, 602a) Reference Pixel Physical Value Y    FD: 0.0
# (0018, 602c) Physical Delta X                    FD: 0.0023148148148148147
# (0018, 602e) Physical Delta Y                    FD: 0.10736213133370061
# (0018, 6032) Pulse Repetition Frequency          UL: 1393

# Check data type and spatial format (should be 3)
# 0003H PW Spectral Doppler
tmp$RegionDataType
# 0003H Spectral (CW or PW Doppler)
tmp$RegionSpatialFormat

x1 <- tmp$RegionLocationMaxX1 + 1
y1 <- tmp$RegionLocationMaxY1 + 1
x0 <- tmp$RegionLocationMinX0 + 1
y0 <- tmp$RegionLocationMinY0 + 1
x_origin = tmp$ReferencePixelX0 + 1
y_origin = tmp$ReferencePixelY0 + 1
slice <- img[y0:y1, x0:x1, ]
# preview_array(slice[,,1])

# get physical unit value
# https://dicom.innolitics.com/ciods/us-image/us-region-calibration/00186011/00186024

# code: 4 - seconds
x_unit_code <- tmp$PhysicalUnitsXDirection
x_unit <- tmp$PhysicalDeltaX  # one pixel is 0.002 seconds

# code: 7 - cm/seconds
y_unit_code <- tmp$PhysicalUnitsYDirection
y_unit <- tmp$PhysicalDeltaY  # one pixel is 0.107 cm/s

# which color channel to get?
slice_idx <- which.max(apply(slice, 3, quantile, 0.99))[[1]]
arr <- slice[,,slice_idx]
# par(mfrow = c(1,2))
# preview_array(arr)

min_val <- 20
# preview_array(arr > min_val)

# get the top values of arr
mask <- arr > min_val
wrapper_idx <- apply(mask, 2, function(x){
  re <- dipsaus::deparse_svec(which(x), concatenate = FALSE)
  if( isTRUE(re == '') ){ return(c(y_origin, y_origin)) }
  l <- sapply(re, function(y){
    z <- dipsaus::parse_svec(y)
    length(z)
  })
  re <- re[which.max(l)]
  if(!length(re)){ return(c(y_origin, y_origin)) }
  re <- paste(re, collapse = ',')
  range(dipsaus::parse_svec(re))
})
# reticulate::repl_python()
top <- (y_origin - wrapper_idx[1,]) * y_unit
bottom <- (y_origin - wrapper_idx[2,]) * y_unit
time <- seq_len(ncol(mask)) * x_unit
yaxis <- -(seq_len(nrow(mask)) - y_origin) * y_unit
s <- top + bottom

# plot
par(mfrow = c(1,1))
plot(time, top + bottom, type = 'l')


# find reverse peaks
valley <- findpeaks(-s, min_distance = 0.2 / x_unit)
split_value <- -min(valley$values)

# split according to split_value
chunks <- lapply(
  dipsaus::deparse_svec(which(s > split_value),
                        concatenate = FALSE),
  function(x){
    x <- dipsaus::parse_svec(x)
    if(length(x) < 0.2 / x_unit){ return(NULL) }
    return(x)
  }
)
chunks <- dipsaus::drop_nulls(chunks)
length(chunks)

# for each chunks, get peak value, and second peak
peak_idx <- sapply(seq_along(chunks), function(ii){
  chunk <- chunks[[ii]]
  x <- s[chunk] - split_value
  peaks <- findpeaks(x, min_distance = 0.05 / x_unit)
  if(!length(peaks$index)){ return(NA) }
  peak_idx <- peaks$index[[1]]
  if(length(peaks$index) >= 2){
    dist <- peaks$index - peak_idx
    sel0 <- dist < 0 & abs(dist) < 0.15 / x_unit
    if(!any(sel0)){
      sel <- dist > 0 & dist < 0.15 / x_unit
      if(any(sel)){
        peak_idx <- peaks$index[sel][[1]]
      }
    }
  }
  chunk[peak_idx]
})
peak_idx <- peak_idx[!is.na(peak_idx)]

# peak is unlikely to be <2 ?
peak_idx <- peak_idx[s[peak_idx] > 1]

par(mfrow = c(1,1))
preview_array(arr, y = yaxis, zlim = c(0, 255))
lines(top, col = 'green', lty = 3)
lines(s, col = 'cyan')
abline(v = peak_idx, col = 'blue')


median(top[peak_idx])




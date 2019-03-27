# DDX6-foci-analysis
# This program extracts DDX6 or DCP1A foci in DAZL positive oocytes and outputs their total size as the number of pixel.
# Setting: Objective lens (100x), zoom setting (1.5x), pixel size (1024x1024)
# channel1: an oocyte marker (DAZL)
# channel2: DDX6 or DCP1A

kern <- makeBrush(3, shape="diamond")
rmask <- NULL

ddx6area <- function(x){
  # Creating the mask of DAZL
  imgg <- x[, , 1] 
  display(normalize(imgg))
  gauss <- gblur(imgg, 99, radius=2*ceiling(3*99)+1)
  blur <- imgg - gauss
  mask1 <- blur > 0
  meanint <- computeFeatures.basic(mask1, imgg)[, 1]
  mask2 <- blur > meanint/10
  mask3 <- closing(opening(mask2, kern),kern)
  label1 <- bwlabel(mask3)
  area1 <- computeFeatures.shape(label1)[,1]
  d100 <- which(area1 < 100)
  remove1 <- rmObjects(label1, d100)
  dazlmask <- remove1 > 0
  area2 <- computeFeatures.shape(dazlmask)[,1]
  display(dazlmask)
  # Creating the mask of DDX6
  imgr <- x[, , 2] 
  display(normalize(imgr))
  gaussr <- gblur(imgr, 50, radius=2*ceiling(3*50)+1)
  blurr <- imgr - gaussr
  blurr2 <- blurr > 0
  meanint <- computeFeatures.basic(blurr2, imgr)[, 1]
  localth <- thresh(blurr, 50, 50, meanint/10)
  gaussr2 <- gblur(imgr, 3, radius=2*ceiling(3*3)+1)
  blurr3 <- imgr - gaussr2
  meanint2 <- computeFeatures.basic(localth, blurr3)[, 1]
  pbodythresh <- 10*meanint2
  pbodymask <- thresh(blurr3, 50, 50, pbodythresh)
  pbodymask2 <- fillHull(pbodymask)
  pbodyindazl <- pbodymask2*dazlmask
  pbodyindazl2 <- pbodyindazl > 0
  pbodyindazl3 <- bwlabel(pbodyindazl2)
  area3 <- computeFeatures.shape(pbodyindazl3)[, 1]
  d1 <- which(area3 < 5) # This value should be changed to 1 in the case of DCP1A foci
  pbodyindazl4 <- rmObjects(pbodyindazl3, d1)
  pbodyindazl5 <- pbodyindazl4 > 0
  area4 <- computeFeatures.shape(pbodyindazl5)[, 1]
  display(pbodyindazl5)
  output <- c(area4, area2)
  mynames <- c("pixels of DDX6 foci", "pixels of oocytes")
  names(output) <- mynames
  return(output)
}

# This program extracts DDX6 foci in DAZL positive oocytes and outputs their size as the number of pixel.
# Setting: Objective lens (100x), zoom (1.5), pixel size (1024x1024)
# channel1: an oocyte marker (DAZL)
# channel2: DDX6

kern <- makeBrush(3, shape="diamond")
rmask <- NULL

ddx6foci <- function(x){
    #Creating the mask of DAZL
    imgg <- x[, , 1] 
    display(normalize(imgg))
    gauss <- gblur(imgg, 99, radius=2*ceiling(3*99)+1)
    blur <- imgg - gauss
    mask1 <- blur > 0
    meanint <- computeFeatures.basic(mask1, imgg)[, 1]
    mask2 <- blur > meanint/10
    mask3 <- closing(opening(mask2, kern),kern)
    label1 <- bwlabel(mask3)
    area1 <- computeFeatures.shape(label1)[,1]
    d100 <- which(area1 < 100)
    remove1 <- rmObjects(label1, d100)
    dazlmask <- remove1 > 0
    display(dazlmask)
    #Creating the mask of DDX6 foci in oocytes
    imgr <- x[, , 2] 
    display(normalize(imgr))
    gaussr <- gblur(imgr, 50, radius=2*ceiling(3*50)+1)
    blurr <- imgr - gaussr
    blurr2 <- blurr > 0
    meanint <- computeFeatures.basic(blurr2, imgr)[, 1]
    localth <- thresh(blurr, 50, 50, meanint/10)
    gaussr2 <- gblur(imgr, 3, radius=2*ceiling(3*3)+1)
    blurr3 <- imgr - gaussr2
    meanint2 <- computeFeatures.basic(localth, blurr3)[, 1]
    pbodythresh <- 10*meanint2
    pbodymask <- thresh(blurr3, 50, 50, pbodythresh)
    pbodymask2 <- fillHull(pbodymask)
    pbodyindazl <- pbodymask2*dazlmask
    pbodyindazl2 <- pbodyindazl > 0
    pbodyindazl3 <- bwlabel(pbodyindazl2)
    for(var in 1: max(pbodyindazl3)){ 
      coord <- which(pbodyindazl3 == var, arr.ind=TRUE)
      for(var2 in 1: nrow(coord)){
        xy <- coord[var2,  ]
        if((xy[1] == 1) || (xy[1] == 1024) || (xy[2] == 1) || (xy[2] == 1024)){
          rmask <- append(rmask, var)
          break} 
      }
    }
    pbodyindazl4 <- rmObjects(pbodyindazl3, rmask)
    area2 <- computeFeatures.shape(pbodyindazl4)[, 1]
    d1 <- which(area2 < 5)
    pbodyindazl5 <- rmObjects(pbodyindazl4, d1)
    pbodyindazl6 <- pbodyindazl5 > 0
    pbodyindazl7 <- bwlabel(pbodyindazl6)
    display(pbodyindazl7)
    pbodysize <- computeFeatures.shape(pbodyindazl7)[, 1]
    return(pbodysize)
}
    
    

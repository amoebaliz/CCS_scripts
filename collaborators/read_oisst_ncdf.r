library(ncdf4)


# fxns
    replace.fillvalues <- function(arr, ncfile, dimname)
    {
        # change fill value from the ncdf file to NA in R
        fillvalue <- ncatt_get(ncfile, dimname, "_FillValue")$value
        arr[arr == fillvalue] <- NA
        arr
    }


    zmatrix.to.cols <- function (x, y=NULL, zmat=NULL, remove.nas=TRUE) 
    {
        # convert matrix to dataframe. 
        
        nmes <- c("x", "y", "z")
        if (is.matrix(x)) 
            x <- data.frame(x)
        if (is.list(x)) {
            if (length(intersect(c("x", "y", "z"), names(x))) != 
                3) {
                nmes <- names(x)[1:3]
                names(x)[1:3] <- c("x", "y", "z")
            }
            zmat <- x$z
            y <- x$y
            x <- x$x
        }
        dat <- expand.grid(x = x, y = y)
        dat$z <- as.vector(zmat)
        if (remove.nas) 
            dat <- dat[!is.na(dat$z), ]
        dat <- dat[order(dat$x, dat$y), ]
        row.names(dat) <- 1:nrow(dat)
        names(dat) <- nmes
        dat
    }


    cols.to.zmatrix <- function (x, y = NULL, z = NULL) 
    {
        # df -> matrix, opposite of zmatrix.to.cols
        nmes <- c("x", "y", "z")
        if (is.matrix(x)) 
            x <- data.frame(x)
        if (is.data.frame(x)) {
            if (length(intersect(c("x", "y", "z"), names(x))) != 
                3) {
                nmes <- names(x)[1:3]
                names(x)[1:3] <- c("x", "y", "z")
            }
            z <- x$z
            y <- x$y
            x <- x$x
        }
        if (any(duplicated(cbind(x, y)))) {
            stop("x-y positions are not unique")
        }
        if ((length(x) != length(y)) | (length(x) != length(z)) | 
            (length(y) != length(z))) {
            stop("x, y and z must be the same length")
        }
        i <- as.numeric(as.factor(x))
        j <- as.numeric(as.factor(y))
        x <- sort(unique(x))
        y <- sort(unique(y))
        mat <- matrix(NA, length(x), length(y))
        for (k in 1:length(z)) {
            mat[i[k], j[k]] <- z[k]
        }
        lst <- list(x, y, mat)
        names(lst) <- nmes
        lst
    }


# read data
    sstfile <- nc_open('oisst_01.nc')

    dteindex <- 1
    dte <- ncvar_get(sstfile, 'time')[dteindex]

    sst <- ncvar_get(sstfile, 'tos', start=c(1, 1, dteindex),
                     count=c(-1, -1, 1))

    # oisst geographics
    lons <- ncvar_get(sstfile, 'lon')
    lats <- ncvar_get(sstfile, 'lat')
 
    #image(lons, lats, sst, col=rev(heat.colors(10)))
    dat <- zmatrix.to.cols(lons, lats, sst, remove.nas=F)
    names(dat) <- c('lon', 'lat', 'sst')
    
    nc_close(sstfile)
 

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
        
        # I think there must be a standard more efficient function
        # than this somewhere but I haven't seen it. This is what I
        # have been using for years.

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
    grdfil <- nc_open('CCS_grd_high_res_bathy_jerlov.nc')
    tfile <- nc_open('SST_10yr_his_clim.nc')
    # sfile <- nc_open('SSS_10yr_his_clim.nc')
    # ufile <- nc_open('SUV_10yr_his_clim.nc')


    # Months are 1 through 12 in the time index: 1 = Jan etc...
    dteindex <- 6
    dte <- ncvar_get(tfile, 'ocean_time')[dteindex]

    # TEMPERATURE FILES
    temp <- ncvar_get(tfile, 'temp', start=c(1, 1, 1, dteindex),
                     count=c(-1, -1, 1, 1))
    # SALINITY FILES
    # salt <- ncvar_get(sstfile, 'salt', start=c(1, 1, 1, dteindex),
    #                 count=c(-1, -1, 1, 1))    

    # U VELOCITY FILES
    # u <- ncvar_get(ufile, 'u', start=c(1, 1, 1, dteindex),
    #                 count=c(-1, -1, 1, 1))

    # V VELOCITY FILES
    # v <- ncvar_get(vfile, 'v', start=c(1, 1, 1, dteindex),
    #                 count=c(-1, -1, 1, 1))

    # read in lat/lons from grid file
    # FOR TEMPERATURE OR SALINITY
    lons <- ncvar_get(grdfil, 'lon_rho')
    lats <- ncvar_get(grdfil, 'lat_rho')

    # FOR U VELOCITY:
    # lons <- ncvar_get(grdfil, 'lon_u')
    # lats <- ncvar_get(grdfil, 'lat_u') 

    # FOR V VELOCITY:
    # lons <- ncvar_get(grdfil, 'lon_v')
    # lats <- ncvar_get(grdfil, 'lat_v') 
    
    # BATHYMETRY (uses lon_rho, lat_rho coordinates)
    bathy  <- ncvar_get(grdfil, 'h')
    # image(bathy)
    dat <- zmatrix.to.cols(lons, lats, sst, remove.nas=F)
    names(dat) <- c('lon', 'lat', 'sst')
    
    nc_close(sstfile)
 

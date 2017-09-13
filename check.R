requiredR <- '3.4.0'

cranPackages <- c('ggplot2'='2.2.1',
                  'dplyr'='0.7.3')

biocPackages <- c('gage'='2.26.1',
                  'edgeR'='3.18.1',
                  'pcaMethods'='1.68.0',
                  'limma'='3.32.6',
                  'org.EcK12.eg.db'='3.4.1')
                  
check <- function() {
    rVersion <- paste(version$major, version$minor, sep=".")
    if(compareVersion(rVersion, requiredR) < 0)
        stop('Your R version is too old. Please download and install R>',
             requiredR, ' from https://cran.r-project.org/mirrors.html')
    fixCran <- packagesToFix(cranPackages)
    fixBioc <- packagesToFix(biocPackages)
    allOk <- TRUE
    for(pkg in fixCran) {
        allOk <- FALSE
        message('version of ', pkg,
                ' is too low/missing. Run install.packages("', pkg, '")')
    }
    for(pkg in fixBioc) {
        allOk <- FALSE
        source('https://bioconductor.org/biocLite.R')
        message('version of ', pkg,
                ' is too low. Run biocLite("', pkg, '")')
    }
    if(!allOk)
        stop('one or more required packages are too old, please update ',
             paste(c(fixCran, fixBioc), collapse=', '))
}

packagesToFix <- function(pkgs) {
    ok <- sapply(seq_along(pkgs), function(i) {
        current <- tryCatch({
            packageVersion(names(pkgs[i]))
        }, error=function(e) '0.0.0')
        compareVersion(as.character(current), pkgs[i]) >= 0
    })
    return(names(pkgs)[!ok])
}

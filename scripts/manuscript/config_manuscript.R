metadatadir = '/central/groups/carnegie_poc/meixilin/grenenet/metadata/data/'
require(ggplot2)
theme_set(cowplot::theme_cowplot())

loadSomeRData <- function(x, file) {
    E <- new.env()
    load(file=file, envir=E)
    return(mget(x, envir=E, inherits=F))
}


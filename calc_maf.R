calc_maf <- function(x)
{
    x <- as.matrix(x)
    
    # Optimize: Vectorized replacement of invalid values
    # Keep only 0, 1, 2; set others to NA
    x[!x %in% c(0, 1, 2)] <- NA
    
    ## calc_n using vectorized rowSums (much faster than apply)
    n0 <- rowSums(x == 0, na.rm = TRUE)
    n1 <- rowSums(x == 1, na.rm = TRUE)
    n2 <- rowSums(x == 2, na.rm = TRUE)
    
    n <- n0 + n1 + n2
    
    ## calculate allele frequencies
    p <- ((2*n0)+n1)/(2*n)
    q <- 1 - p
    maf <- pmin(p, q)

    # MODIFIED 21 Oct 2012:  prior to this version, we had "mono=(mgf<0)" instead of "mono<(maf<0)"
    res <- data.frame( n=n, n0=n0, n1=n1, n2=n2, p=p, maf=maf,
                       mono=(maf<=0), loh=(n1<=0), 
                       stringsAsFactors=F )
    row.names(res) <- row.names(x)
    res
}

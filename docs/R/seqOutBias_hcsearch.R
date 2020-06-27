#
# Hill-Climbing Search
#
source("https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_evalmask.R")
require(parallel)

#' Convert a mask in string form to a kmer-mask in integer form
#' Map:
#' -1 = X :: positions that do not contribute to k-mer
#'  0 = C :: marker for cut-site position in k-mer mask
#'  1 = N :: positions that contribute to k-mer
#' @return integer vector
mask.to.vec <- function(mask) {
    mask = toupper(mask)
    sapply(charToRaw(mask), function(c) {
        if (c == 88) { -1 } # X
        else if (c == 78) { 1 } # N
        else if (c == 67) { 0 } # C
        else stop(paste("Unknown caracter in mask:", rawToChar(c)))
    })
}

#' Convert a mask in integer form to a kmer-mask in string form
#' Map:
#' -1 = X :: positions that do not contribute to k-mer
#'  0 = C :: marker for cut-site position in k-mer mask
#'  1 = N :: positions that contribute to k-mer
#' @return character string
vec.to.mask <- function(maskvec) {
    rawToChar(as.raw(sapply(maskvec, function(code) {
        if (code == -1) { 88 } # X
        else if (code == 1) { 78 } # N
        else if (code == 0) { 67 } # C
        else stop(paste("Unknown code in maskvec:", code))
    })))
}

# Start with a fully opaque kmer-mask and add one useful base per step
# picking the one that gives the greatest decrease in score.
hc.search <- function(sites, start_mask, seqoutbias.args, prefix = "run_", bw.split = FALSE, cleanup = TRUE, sqcmd = "seqOutBias", mc.cores = 2) {
    neighbors <- function(vecmask) {
         # generate list of neighboring masks that differ
         # by the addition of a single unmasked position
         result <- vector(mode="list", length=length(vecmask))
         n <- 0
         for (k in 1:length(vecmask)) {
             if (vecmask[k] == -1) { # X
                n <- n + 1
                vi = vecmask
                vi[k] = 1 # N
                result[[n]] <- vec.to.mask(vi)
             }
         }

         if (n == 0) return(NULL)
         result[1:n]
    }

    #
    mask = character()
    score = numeric()
    rtable = data.frame(c(score, mask), stringsAsFactors=FALSE)
    nextlst = neighbors(mask.to.vec(start_mask))
    while (!is.null(nextlst)) {
        # score ecah neighbor in parallel
        scores = mclapply(nextlst, function(cutmask) {
            bw.paths = run.cutmask(cutmask, seqoutbias.args, sqcmd=sqcmd, prefix=prefix)
            bw.plus = load.bigWig(bw.paths[1])
            bw.minus = load.bigWig(bw.paths[2])
            score = eval.cutmask(sites, bw.plus, bw.minus)
            #score = runif(1)*1000
            data.frame(score = score, mask = cutmask, stringsAsFactors=FALSE)
        }, mc.cores = mc.cores)
        scores.values = sapply(scores, function(pair) pair$score)

        # pick the neighbor with the best score
        k = which.min(scores.values)[1]
        mask = scores[[k]]$mask
        cat(scores.values[k], mask, "\n")
        rtable = rbind(rtable, scores[[k]])

        # next
        nextlst <- neighbors(mask.to.vec(mask))
    }

    rtable
}

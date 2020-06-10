#' Get posteriors
#'
#' @param dta data frame that contains group identifier and binary outcome
#' @param method \code{powerp} for power prior; \code{{mixp} for mixture prior
#' @param var_grp column name in dta for the group identifier
#' @param var_outcome column name in dta for the binary outcome
#' @param vague_ab alpha and beta parameter for the vague beta prior. Defaults
#'     are 0.5, 0.5
#' @param grp_hist group identifiers for historical studies
#' @param grp_cur group identifiers for current study
#' @param weights Depending on the method, power parameters or mixture weights.
#'     Length should be the same as grp_hist. All weights are non-negative. For
#'     mixture weights, the sum has to be less than or equal to 1.
#' @param ... Parameters for STAN sampling
#' @return Posterior samples of theta
#'
#' @examples
#' \dontrun{
#' dta_1 <- data.frame(group = 1, y = rbinom(100, 1, 0.8))
#' dta_2 <- data.frame(group = 2, y = rbinom(100, 1, 0.6))
#' dta_3 <- data.frame(group = 3, y = rbinom(100, 1, 0.7))
#' dta   <- rbind(dta_1, dta_2, dta_3)
#'
#' synSample(dta, method = "powerp",
#'           var_grp = "group", var_outcome = "y",
#'           vague_ab = c(0.5, 0.5),
#'           grp_hist = c(1, 2), grp_cur = 3, weights = c(0, 0))
#'
#' synSample(dta, method = "mixp",
#'           var_grp = "group", var_outcome = "y",
#'           vague_ab = c(0.5, 0.5),
#'           grp_hist = c(1, 2), grp_cur = 3, weights = c(0, 0.9))}
#'
#'
#' @export
#'
#'
synSample <- function(dta, method = c("powerp", "mixp"), var_grp = "group",
                      var_outcome = "Y", vague_ab = c(0.5, 0.5),
                      grp_hist = NULL, grp_cur = NULL,
                      weights = NULL, ...) {

    method <- match.arg(method)

    stopifnot(!is.null(dta[[var_grp]]) &
              !is.null(dta[[var_outcome]]))

    if (!is.null(grp_hist) & !is.null(weights))
        stopifnot(length(grp_hist) == length(weights))

    grps <- sort(unique(dta[[var_grp]]))

    if (is.null(grp_hist)) {
        grp_hist <- grps[-length(grps)]
    } else {
        stopifnot(all(grp_hist %in% grps))
    }

    if (is.null(grp_cur)) {
        grp_cur <- grps[length(grps)]
    } else {
        stopifnot(all(grp_cur %in% grps))
    }


    ## prepare data
    S        <- 1 + length(grp_hist)

    ## weights
    if ("powerp" == method) {
        if (is.null(weights)) {
            weights <- rep(0, length(grp_hist))
        }
        WEIGHTS <- c(1, weights)
    } else {
        if (is.null(weights)) {
            WEIGHTS <- rep(1/S, S)
        } else {
            stopifnot(sum(weights) <= 1)
            WEIGHTS <- c(1 - sum(weights), weights)
        }
    }


    ## vague
    N     <- sum(vague_ab)
    YSUM  <- vague_ab[1]

    if (length(grp_hist) > 0) {
        sn <- NULL
        sy <- NULL
        for (i in seq_len(length(grp_hist))) {
            inx  <- which(grp_hist[i] == dta[[var_grp]])
            sn   <- c(sn, length(inx))
            sy   <- c(sy, sum(dta[[var_outcome]][inx]))
        }

        if ("mixp" == method) {
            ## get posterior
            sn <- sn + sum(vague_ab)
            sy <- sy + vague_ab[1]
        }

        N    <- c(N, sn)
        YSUM <- c(YSUM, sy)
    }

    inx     <- which(dta[[var_grp]] %in% grp_cur)
    NCUR    <- length(inx)
    YSUMCUR <- sum(dta[[var_outcome]][inx])

    stan_data <- list(S = S, N = N, YSUM = YSUM, WEIGHTS = WEIGHTS,
                      NCUR = NCUR, YSUMCUR = YSUMCUR)

    print(stan_data)

    ## calll
    synSTAN(stan_data, stan.mdl = method, ...)
}

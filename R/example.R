#' Get posteriors
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

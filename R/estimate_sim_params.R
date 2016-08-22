estGammaParams <- function(counts) {

    # Check counts
    counts <- as.matrix(counts)

    # Normalise for differing library sizes
    lib.sizes <- colSums(counts)
    lib.med <- median(lib.sizes)
    norm.counts <- t(t(counts) / lib.sizes * lib.med)

    means <- rowMeans(norm.counts)
    means <- means[means != 0]

    z <- log(means)

    mean.z <- mean(z)
    var.z <- var(z)

    alpha <- limma::trigammaInverse(var.z)
    beta <- exp(digamma(alpha) - mean.z)

    estimates <- list("Shape" = alpha, "Rate" = beta)

    return(estimates)
}

estLogisticParams <- function(counts) {
    # Check counts
    counts <- as.matrix(counts)

    logistic <- function(x, x0, k) {1 / (1 + exp(-k * (x - x0)))}

    lib.sizes <- colSums(counts)
    lib.med <- median(lib.sizes)
    norm.counts <- t(t(counts) / lib.sizes * lib.med)

    means <- rowMeans(norm.counts)

    x <- log(means)

    y <- rowSums(counts == 0) / ncol(counts)

    df <- data.frame(x, y)

    fit <- nls(y ~ logistic(x, x0 = x0, k = k), data = df,
               start = list(x0 = 0, k = -1))

    mid <- summary(fit)$coefficients["x0", "Estimate"]
    shape <- summary(fit)$coefficients["k", "Estimate"]

    estimates <- list("Midpoint" = mid, "Shape" = shape)

    return(estimates)
}


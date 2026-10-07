library(bnlearn)
library(arrow)
library(parallel)

num_boots <- 100
num_cores <- 24
subsample <- 0.75
threshold <- 0.5

entropy <- function(x) {
    p <- table(x) / length(x)
    p <- p[p > 0]
    return(-sum(p * log(p)))
}

test <- function(center, set, data) {
    sapply(set, function(g) {
        tt <- ci.test(x = center, y = g, z = setdiff(set, g), data = data, test = "mi")
        mi = unname(tt$statistic / (2 * nrow(data)))
        nmi = mi / (entropy(data[[center]]) * entropy(data[[g]]))^0.5
        return(c(statistic = unname(tt$statistic), df = unname(tt$parameter[1]), mi = mi, nmi = nmi, pval = tt$p.value))
    })
}

process <- function(name, source_path) {
    set.seed(0)

    data <- read_feather(source_path)
    data <- as.data.frame(lapply(data, factor))

    indices <- lapply(seq_len(num_boots), function(i) sample.int(nrow(data), round(subsample * nrow(data))))
    boots <- mclapply(
        indices,
        function(index) {message(1); learn.mb(x = data[index, , drop = FALSE], node = "gravity", method = "iamb.fdr", test = "mi")},
        mc.cores = num_cores
    )

    candidates <- unique(unlist(boots))
    stabilities <- sapply(candidates, function(g) mean(vapply(boots, function(b) g %in% b, logical(1))))
    consensus <- names(stabilities)[stabilities >= threshold]

    results <- data.frame(gene = consensus, stability = stabilities[consensus], t(test("gravity", consensus, data)))
    write.csv(results, paste("data/blanket_", name, ".csv", sep=""), row.names=FALSE)
    saveRDS(list(indices = indices, boots = boots), paste("data/boots_", name, ".rds", sep=""))
}

process("402", "data/discretized_402.feather")
process("403", "data/discretized_403.feather")
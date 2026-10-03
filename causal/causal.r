library(bnlearn)
library(arrow)
library(parallel)
library(tidyverse)

num_boots <- 100
num_cores <- 2
subsample_proportion <- 0.8

process <- function(name, source_path) {
    set.seed(42)

    data <- read_feather(source_path)
    data <- as.data.frame(lapply(data, factor))
    group_size <- count(data, gravity, age, name="count") %>% select(count) |> min()

    bootstrap <- mclapply(1:num_boots, function(i) {
        message("Bootstrapping ", i, "/", num_boots)
        subset <- group_by(data, gravity, age) %>% slice_sample(n=round(group_size * subsample_proportion)) %>% ungroup() |> as.data.frame()
        boot <- learn.mb(x = subset, node = "gravity", method = "iamb.fdr", test = "mi")
        return(boot)
    }, mc.cores = num_cores)

    mb <- unique(unlist(bootstrap))

    message("Calculating statistics")
    results <- sapply(mb, function(g) {
        test <- ci.test(x = "gravity", y = g, z = setdiff(mb, g), data = data, test = "mi")
        stability <- sum(sapply(bootstrap, function(x) g %in% x)) / length(bootstrap)
        return(c(mi = test$statistic[["mi"]], pval = test$p.value, stability = stability))
    })

    write.csv(t(results), paste("data/", name, ".csv", sep=""))
}

process("402", "data/discretized_402.feather")
process("403", "data/discretized_403.feather")
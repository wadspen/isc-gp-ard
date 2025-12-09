library(dplyr)
library(stringr)

files <- list.files(pattern = ".rds")
ress <- list()
for (i in 1:length(files)) {
	ress[[i]] <- readRDS(files[i])
	ress[[i]]$proc <- str_replace(files[i], ".rds", "")
}

res <- do.call(rbind, ress)
saveRDS(res, "comb_res.rds")

# Build gh-pages docs and copy image files

pkgdown::build_site()

vigDir <- "vignettes"
artDir <- "docs/articles"
pFiles <- list.files(vigDir, pattern = "page")
file.copy(file.path(vigDir, pFiles), file.path(artDir, pFiles))

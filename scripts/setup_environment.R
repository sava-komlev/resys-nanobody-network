# Run from the project directory. Copies dependencies; never changes global libraries.
dir.create("environment/library", recursive = TRUE, showWarnings = FALSE)
private <- normalizePath("environment/library", mustWork = TRUE)
db <- installed.packages()
roots <- c("Rcpp", "scales", "igraph", "Biostrings", "rpart", "MASS", "jsonlite")
deps <- unique(c(roots, unlist(tools::package_dependencies(roots, db = db,
                which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE))))
deps <- intersect(deps, rownames(db))
deps <- deps[is.na(db[deps, "Priority"]) | db[deps, "Priority"] != "base"]
for (p in deps) {
  dest <- file.path(private, p)
  if (!dir.exists(dest)) stopifnot(file.copy(file.path(db[p, "LibPath"], p), private, recursive = TRUE))
}
write.csv(db[deps, c("Package", "Version", "LibPath", "Built"), drop = FALSE],
          "environment/copied_packages.csv", row.names = FALSE)
cat("Copied", length(deps), "packages into", private, "\n")

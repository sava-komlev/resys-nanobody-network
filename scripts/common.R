options(stringsAsFactors = FALSE)
private <- normalizePath("environment/library", mustWork = TRUE)
.libPaths(c(private, .Library))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
set.seed(20261030)
for (d in c("data", "results", "figures", "logs")) dir.create(d, showWarnings = FALSE)
log_line <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush.console() }
save_csv <- function(x, path) write.csv(x, path, row.names = FALSE)
assert_private <- function(pkgs) {
  for (p in pkgs) stopifnot(startsWith(normalizePath(find.package(p)), private))
}
edge_key <- function(a,b) ifelse(a < b, paste(a,b,sep="--"), paste(b,a,sep="--"))
jaccard <- function(a,b) if (!length(union(a,b))) 1 else length(intersect(a,b))/length(union(a,b))

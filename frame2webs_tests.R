source("./bipartite/R/frame2webs.R")
library(microbenchmark)

# Use the DBIF data base for a list of species names and web IDs
# See https://data-package.ceh.ac.uk/data/33a825f3-27cb-4b39-b59c-0f8182e8e2e4
# to download the data.
dbif <- read.csv("./DBIF/data/DBIFInteractionFrequencies.csv")

# Extract the names from the data base and convert them to factor-levels
lower_species_names <- levels(factor(dbif$Plant.Species.Name))
higher_species_names <- levels(factor(dbif$Insect.Species.Name))
web_id_names <- levels(factor(dbif$DBIF.Source))

# Free space by removing the DBIF data.frame
rm(dbif)

frame2webs_test <- function(lower_species_names,
                            higher_species_names,
                            web_id_names,
                            n_lower,
                            n_higher,
                            n_ids,
                            n_total) {
  t_lower  <- sample(lower_species_names, n_lower)
  t_higher <- sample(higher_species_names, n_higher)
  t_id     <- sample(web_id_names, n_ids)
  t_frame <- data.frame(higher = sample(t_higher, n_total, replace = TRUE),
                        lower  = sample(t_lower, n_total, replace = TRUE),
                        webID  = sample(t_id, n_total, replace = TRUE),
                        freq   = sample.int(10, n_total, replace = TRUE,
                                            prob = 1 / (1:10)))
  web_old  <- frame2webs(t_frame, type.out = "array", emptylist = FALSE)
  web_fast <- frame2webs_fast(t_frame, type.out = "array", emptylist = FALSE)

  check_v <- identical(web_old, web_fast)
  check_r <- identical(rownames(web_old), rownames(web_fast))
  check_c <- identical(colnames(web_old), colnames(web_fast))

  rm(web_old, web_fast, t_frame)

  return(check_v & check_r & check_c)
}

frame2webs_test_old <- function(lower_species_names,
                             higher_species_names,
                             web_id_names,
                             n_lower,
                             n_higher,
                             n_ids,
                             n_total) {
  t_lower  <- sample(lower_species_names, n_lower)
  t_higher <- sample(higher_species_names, n_higher)
  t_id     <- sample(web_id_names, n_ids)
  t_frame <- data.frame(higher = sample(t_higher, n_total, replace = TRUE),
                        lower  = sample(t_lower, n_total, replace = TRUE),
                        webID  = sample(t_id, n_total, replace = TRUE),
                        freq   = sample.int(10, n_total, replace = TRUE,
                                            prob = 1 / (1:10)))
  web  <- frame2webs(t_frame, type.out = "array", emptylist = FALSE)
  rm(web, t_frame)
}
frame2webs_test_fast <- function(lower_species_names,
                             higher_species_names,
                             web_id_names,
                             n_lower,
                             n_higher,
                             n_ids,
                             n_total) {
  t_lower  <- sample(lower_species_names, n_lower)
  t_higher <- sample(higher_species_names, n_higher)
  t_id     <- sample(web_id_names, n_ids)
  t_frame <- data.frame(higher = sample(t_higher, n_total, replace = TRUE),
                        lower  = sample(t_lower, n_total, replace = TRUE),
                        webID  = sample(t_id, n_total, replace = TRUE),
                        freq   = sample.int(10, n_total, replace = TRUE,
                                            prob = 1 / (1:10)))
  web  <- frame2webs_fast(t_frame, type.out = "array", emptylist = FALSE)
  rm(web, t_frame)
}

n_tests <- 100

####################
#### EDGE CASES ####
####################
print("Testing edge cases...")

### FIRST EDGE CASE
cat("Only one higher and lower species as well as web ID and one entry:")
pb <- txtProgressBar(min = 0, max = n_tests, initial = 0,
                     style = 3, title = "Just progress")
all_res <- array(FALSE, n_tests)
for (i in 1:n_tests) {
  setTxtProgressBar(pb, i)
  res <- frame2webs_test(lower_species_names,
                         higher_species_names,
                         web_id_names,
                         n_lower = 1,
                         n_higher = 1,
                         n_ids = 1,
                         n_total = 1)
  all_res[i] <- res
}
print(paste0("Pass: ", all(all_res)))

### SECOND EDGE CASE
cat("Only one higher and lower species,",
    "as well as web ID and up to 100.000 entries:")
pb <- txtProgressBar(min = 0, max = n_tests, initial = 0,
                     style = 3, title = "Just progress")
all_res <- array(FALSE, n_tests)
for (i in 1:n_tests) {
  setTxtProgressBar(pb, i)
  res <- frame2webs_test(lower_species_names,
                         higher_species_names,
                         web_id_names,
                         n_lower = 1,
                         n_higher = 1,
                         n_ids = 1,
                         n_total = sample.int(1e5, 1))
  all_res[i] <- res
}
print(paste0("Pass: ", all(all_res)))

### THIRD EDGE CASE
cat("Only one higher species up to 1000 lower species,",
    "one web ID and up to 100.000 entries:")
pb <- txtProgressBar(min = 0, max = n_tests, initial = 0,
                     style = 3, title = "Just progress")
all_res <- array(FALSE, n_tests)
for (i in 1:n_tests) {
  setTxtProgressBar(pb, i)
  res <- frame2webs_test(lower_species_names,
                         higher_species_names,
                         web_id_names,
                         n_lower = sample.int(1e3, 1),
                         n_higher = 1,
                         n_ids = 1,
                         n_total = sample.int(1e5, 1))
  all_res[i] <- res
}
print(paste0("Pass: ", all(all_res)))


### FOURTH EDGE CASE
cat("Only lower species up to 1000 higher species,",
    "one web ID and up to 100.000 entries:")
pb <- txtProgressBar(min = 0, max = n_tests, initial = 0,
                     style = 3, title = "Just progress")
all_res <- array(FALSE, n_tests)
for (i in 1:n_tests) {
  setTxtProgressBar(pb, i)
  res <- frame2webs_test(lower_species_names,
                         higher_species_names,
                         web_id_names,
                         n_lower = 1,
                         n_higher = sample.int(1e3, 1),
                         n_ids = 1,
                         n_total = sample.int(1e5, 1))
  all_res[i] <- res
}
print(paste0("Pass: ", all(all_res)))

### FIFTH EDGE CASE
cat("Only one higher and lower species,",
    "up to 600 web IDs and up to 100.000 entries:")
pb <- txtProgressBar(min = 0, max = n_tests, initial = 0,
                     style = 3, title = "Just progress")
all_res <- array(FALSE, n_tests)
for (i in 1:n_tests) {
  setTxtProgressBar(pb, i)
  res <- frame2webs_test(lower_species_names,
                         higher_species_names,
                         web_id_names,
                         n_lower = 1,
                         n_higher = 1,
                         n_ids = sample.int(600, 1),
                         n_total = sample.int(1e5, 1))
  all_res[i] <- res
}
print(paste0("Pass: ", all(all_res)))

### FULL TEST
print("Full Test:")
pb <- txtProgressBar(min = 0, max = n_tests, initial = 0,
                     style = 3, title = "Just progress")
all_res <- array(FALSE, n_tests)
for (i in 1:n_tests) {
  setTxtProgressBar(pb, i)
  res <- frame2webs_test(lower_species_names,
                         higher_species_names,
                         web_id_names,
                         n_lower = sample.int(200, 1),
                         n_higher = sample.int(200, 1),
                         n_ids = sample.int(200, 1),
                         n_total = sample.int(1e5, 1))
  all_res[i] <- res
}
print(paste0("Pass: ", all(all_res)))

##################
### BENCHMARKS ###
##################
benchmark <- microbenchmark(
    frame2webs_test_old(lower_species_names, higher_species_names,
                        web_id_names, 1, 1, 1, 1e5),
    frame2webs_test_fast(lower_species_names, higher_species_names,
                        web_id_names, 1, 1, 1, 1e5),
    times = 100, unit = "ms")
print(benchmark)
benchmark <- microbenchmark(
    frame2webs_test_old(lower_species_names, higher_species_names,
                        web_id_names, 200, 200, 10, 5e4),
    frame2webs_test_fast(lower_species_names, higher_species_names,
                        web_id_names, 200, 200, 10, 5e4),
    times = 100, unit = "ms")
print(benchmark)
benchmark <- microbenchmark(
    frame2webs_test_old(lower_species_names, higher_species_names,
                        web_id_names, 200, 200, 100, 5e4),
    frame2webs_test_fast(lower_species_names, higher_species_names,
                        web_id_names, 200, 200, 100, 5e4),
    times = 100, unit = "ms")
print(benchmark)

rm(lower_species_names, higher_species_names, web_id_names)
## code to prepare `DATASET` dataset goes here

dat.Boucher2016.pairwise <- cbind.data.frame("study" = c("Edwards (2000)", "Storey (2001)",
                                                         "Brandes (2004)", "Diener (2004)",
                                                         "Silberstein (2004)",
                                                         "Silberstein (2006)"),
                                             "duration" = c("short","long")[c(1,1,2,2,2,1)],
                                             "r1" = c(4, 8, 9, 5, 4, 4),
                                             "n1" = c(73, 116, 143, 113, 21, 15),
                                             "r2" = c(63, 53, 81, 57, 13, 9),
                                             "n2" = c(140, 113, 144, 117, 19, 15))


usethis::use_data(dat.Boucher2016.pairwise,
                  overwrite = TRUE, version = 3)

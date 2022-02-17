#remotes::install_github("rqthomas/cronR")
#remotes::install_deps()
library(cronR)

home_dir <- "/home/rstudio/Documents/scripts"
log_dir <- "/efi_neon_challenge/log/cron"

phenology_repo <- "neon4cast-phenology"

#Go to healthchecks.io. Create a project.  Add a check. Copy the url and add here.  
health_checks_url <- "https://hc-ping.com/a5a12c66-6e38-4415-aa40-4eb24881a949"

cmd <- cronR::cron_rscript(rscript = file.path(home_dir, phenology_repo,"phenology-nullmodels.R"),
                           rscript_log = file.path(log_dir, "phenology-null.log"),
                           log_append = FALSE,
                           workdir = file.path(home_dir, phenology_repo),
                           trailing_arg = paste0("curl -fsS -m 10 --retry 5 -o /dev/null ", health_checks_url))
cronR::cron_add(command = cmd, frequency = 'daily', at = "10PM", id = 'phenology-nullmodel')

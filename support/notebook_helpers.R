# Source from the lab folder. Python runs only after the editor has closed.
.lab_support <- normalizePath("support", mustWork=TRUE)
source(file.path(.lab_support, "landscape_editor.R"))
source(file.path(.lab_support, "plots.R"))

gnx_parameters <- function(phi=0.05, n_loci=1, K_factor=5,
                           movement_distance=0.5, mating_radius=4, seed=154) {
  values <- c(phi,n_loci,K_factor,movement_distance,mating_radius,seed)
  if (length(values)!=6 || any(!is.finite(values)) ||
      phi<0 || phi>1 || !(n_loci %in% 1:100) ||
      any(values[3:5]<=0) || seed<1 || seed!=floor(seed))
    stop("Check the simulation parameters.")
  list(phi=phi, n_loci=n_loci, K_factor=K_factor,
       movement_distance=movement_distance, mating_radius=mating_radius, seed=seed)
}

run_geonomics <- function(carrying_capacity, connectivity, environment,
                           params=gnx_parameters(), steps=200, name="simulation",
                           rerun=interactive(), background=FALSE) {
  if (length(name)!=1 || !grepl("^[A-Za-z0-9_-]+$", name))
    stop("Use only letters, digits, underscores, and hyphens in the run name.")
  folder <- file.path("runs", name)
  result_file <- file.path(folder,"result.rds")
  if (!rerun && file.exists(result_file)) return(readRDS(result_file))
  if (!interactive() && !rerun)
    stop("Run this simulation interactively before knitting: ", name)
  if (!length(steps) || any(!is.finite(steps)) || any(steps<1 | steps!=floor(steps)))
    stop("steps must contain positive whole numbers.")
  layers <- list(capacity=as.matrix(carrying_capacity),
                 connectivity=as.matrix(connectivity), environment=as.matrix(environment))
  dims <- lapply(layers,dim)
  if (!identical(dims[[1]],dims[[2]]) || !identical(dims[[1]],dims[[3]]) ||
      any(vapply(layers,function(x) !is.numeric(x) || any(!is.finite(x)) ||
                   any(x<0 | x>1),logical(1))))
    stop("All three landscapes must have matching dimensions and values from 0 to 1.")
  if (!any(layers$capacity>0)) stop("Carrying capacity cannot be zero everywhere.")
  # A new directory prevents a failed rerun from leaving stale successful results.
  dir.create("runs",showWarnings=FALSE)
  if (dir.exists(folder)) folder <- tempfile(paste0(name,"_"), tmpdir="runs")
  dir.create(folder)
  for (nm in names(layers))
    write.table(layers[[nm]], file.path(folder,paste0(nm,".csv")),
                sep=",",row.names=FALSE,col.names=FALSE)
  py <- Sys.getenv("EVODOODLE_PYTHON", "/srv/conda/envs/notebook/bin/python")
  if (!file.exists(py)) stop("Python not found: ",py)
  vals <- unlist(params[c("phi","n_loci","K_factor","movement_distance","mating_radius","seed")],
                 use.names=FALSE)
  if (length(vals)!=6) stop("Use gnx_parameters() to create params.")
  lib <- Sys.getenv("EVODOODLE_PYTHON_LIB",file.path(dirname(dirname(py)),"lib"))
  env <- "MPLBACKEND=Agg"
  if (.Platform$OS.type!="windows") env <- c(env,paste0("LD_LIBRARY_PATH=",shQuote(lib)))
  log <- file.path(folder,"python.log")
  folder <- normalizePath(folder)
  log <- file.path(folder,"python.log")
  answer <- structure(list(folder=file.path("runs",basename(folder)),layers=layers,stages=NULL,params=params),
                      class="lab_geonomics_run")
  args <- c("-u",file.path(.lab_support,"run_geonomics.py"),folder,
            paste(steps,collapse=","),as.character(vals))
  if (background) stop("Run simulations sequentially with background=FALSE.")
  message("Running Geonomics. Log: ",log)
  exit <- system2(py,args=shQuote(args),env=env,stdout=log,stderr=log)
  writeLines(as.character(exit),file.path(folder,"exit_status"))
  answer <- completed_run(answer)
  message("Simulation complete.")
  saveRDS(answer,file.path(folder,"result.rds"))
  dir.create(file.path("runs",name),showWarnings=FALSE)
  saveRDS(answer,file.path("runs",name,"result.rds"))
  answer
}

simulation_status <- function(run) {
  sentinel <- file.path(run$folder,"exit_status")
  if (file.exists(sentinel)) {
    code <- suppressWarnings(as.integer(readLines(sentinel,warn=FALSE)[1]))
    if (is.na(code) || code!=0) {
      cat("Simulation failed. Last log lines:\n")
      cat(tail(readLines(file.path(run$folder,"python.log"),warn=FALSE),25),sep="\n")
      return(invisible("failed"))
    }
    cat("Simulation complete. You can run the plotting chunk.\n")
    return(invisible("complete"))
  }
  # Older saved runs did not have an exit-status file.
  if (!is.null(run$stages) && file.exists(file.path(run$folder,"stages.csv"))) {
    cat("Simulation complete.\n")
    return(invisible("complete"))
  }
  progress <- file.path(run$folder,"progress.txt")
  cat(if(file.exists(progress)) readLines(progress,warn=FALSE) else "Waiting for Python to start…",sep="\n")
  cat("\nRun this status chunk again to check progress.\n")
  invisible("running")
}

completed_run <- function(run) {
  status <- simulation_status(run)
  if (status!="complete") stop("Simulation is ",status,"; check simulation_status(run).",call.=FALSE)
  run$stages <- read.csv(file.path(run$folder,"stages.csv"))
  run
}

plot_popgen <- function(run, stage=NULL) {
  run <- completed_run(run)
  if (is.null(stage)) stage <- nrow(run$stages)-1L
  row <- match(stage,run$stages$stage)
  if (is.na(row)) stop("stage must be one of: ",paste(run$stages$stage,collapse=", "))
  plot_snapshot(read_snapshot(run$folder,stage),run$layers,run$stages$steps[row])
  invisible(run)
}

plot_landscapes <- function(carrying_capacity,connectivity,environment) {
  old <- par(mfrow=c(1,3),mar=c(2,2,3,1))
  on.exit(par(old))
  layers <- list("Carrying capacity"=carrying_capacity,
                 "Connectivity"=connectivity,"Environment"=environment)
  palettes <- list(
    c("#f7fbff","#c6dbef","#6baed6","#2171b5","#08306b"),
    c("#f7fcf5","#c7e9c0","#74c476","#238b45","#00441b"),
    c("#542788","#998ec3","#f7f7f7","#f1a340","#b35806"))
  for (i in seq_along(layers)) {
    nm <- names(layers)[i]
    x <- layers[[nm]]
    image(seq_len(ncol(x))-.5,seq_len(nrow(x))-.5,t(x),
      col=colorRampPalette(palettes[[i]])(257),zlim=c(0,1),ylim=c(nrow(x),0),
      xlim=c(0,ncol(x)),asp=1,axes=FALSE,xlab="",ylab="",main=nm)
    box()
  }
}

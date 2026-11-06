# Landscape simulation how-to

A standalone R Markdown guide for students considering a simulation-based final
project. This folder contains the reusable tools and a walkthrough, independent
of the Lab 9 scenario exercise.

Upload the contents to GitHub, preserving support/. Students open the .Rproj
in DataHub RStudio, then Landscape_Simulation_HowTo.Rmd and run chunks in order.
The guide covers three-layer Shiny drawing, editable parameters, sequential
Geonomics runs, snapshots, VCFs, coordinates, replicate design, and reopening work.
Optional editing, vcfR import, and replicate chunks have eval=FALSE.

Requires R packages shiny, jsonlite, knitr, rmarkdown and the existing DataHub
Geonomics Python installation. vcfR is optional for importing exported VCFs.
The template uses 1,000 genomic loci and starts with 500 individuals before burn-in;
recombination settings remain unchanged. No completed simulation data are included.

Code syntax and snapshot export wiring were checked locally. The core 200-step
configuration was tested by Ian on DataHub; the added native VCF export still
needs a fresh end-to-end DataHub run before student distribution.

See support/LICENSE-parameter-template.txt for the original Evodoodle template license.

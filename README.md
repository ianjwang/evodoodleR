# Landscape simulation using evodoodleR and geonomics

A standalone R Markdown guide for anyone interested in learning how to use evodoodleR. 
This repo contains the reusable tools and a walkthrough.

The guide covers three-layer Shiny drawing, editable parameters, sequential
Geonomics runs, snapshots, VCFs, coordinates, replicate design, and reopening work.
Optional editing, vcfR import, and replicate chunks have eval=FALSE.

Requires R packages shiny, jsonlite, knitr, rmarkdown and the a Geonomics
Python installation. vcfR is optional for importing exported VCFs.

The template uses 1,000 genomic loci and starts with 500 individuals before burn-in.

See support/LICENSE-parameter-template.txt for the original Evodoodle template license.

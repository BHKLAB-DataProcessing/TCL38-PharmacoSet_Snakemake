#!env bash
set -o pipefail

R -e "pak::pkg_install(c('bhklab/CoreGx','bhklab/AnnotationGx'), dependencies = TRUE)"
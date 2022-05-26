# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:4.2.0

#Install packages
RUN RUN R -e 'install.packages("Seurat",repos=("http://cran.rstudio.com")); install.packages("spatstat.sparse",repos=("http://cran.rstudio.com"));install.packages('argparser', repos=("http://cran.rstudio.com"))'
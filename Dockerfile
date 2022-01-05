FROM rocker/rstudio:4.1.2
ENV RENV_VERSION 0.14.0


RUN sudo apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    libpng-dev \
    libgsl-dev \ 
    apt-transport-https \
    build-essential \
    curl \
    gfortran \
    libatlas-base-dev \
    libbz2-dev \
    libcairo2 \
    libcurl4-openssl-dev \
    libicu-dev \
    liblzma-dev \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libpcre3-dev \
    libtcl8.6 \
    libtiff5 \
    libtk8.6 \
    libx11-6 \
    libxt6 \
    locales \
    tzdata \
    zlib1g-dev \ 
    libcurl4-openssl-dev \ 
    libssl-dev \
    pandoc \
    libxml2-dev

RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"

RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY renv.lock renv.lock
RUN R -e 'renv::restore()'

# Base Image: Contains R version 4.2.2 and Rstudio Server
FROM rocker/r-ver:4.2.2

# Use latest version of Renv (version 0.16.0) from Github
ENV RENV_VERSION 0.16.0

# Install the packages needed for many of the libraries (including RENV)
# Afterwords we need to do some cleaning
RUN --mount=type=cache,target=/var/cache/apt \
    apt-get update && \
    apt-get install -y --no-install-recommends "$@" \
    procps \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev \
    default-libmysqlclient-dev \
    libpq-dev \
    libsasl2-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    libxtst6 \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libglpk-dev \
    unixodbc-dev \
    libpcre3-dev \
    zlib1g-dev \
    libbz2-dev && \
    rm -rf /var/lib/apt/lists/* && \
    Rscript -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    Rscript -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')" && \
    mkdir -p renv

# Repare to install packages by environment variables
COPY renv.lock renv.lock
ENV RENV_PATHS_LIBRARY renv/library

# Install packages and do some cleanup
RUN R -e 'renv::restore()' && \
    rm -rf /tmp/downloaded_packages && \
    strip /usr/local/lib/R/site-library/*/libs/*.so

# Set the working directory for the project
WORKDIR /home/rstudio

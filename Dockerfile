# 1. Base Image: Use the latest RStudio image built on a Debian-based Linux.
FROM rocker/rstudio:4.3.2

## --------------------------------------------------------------------------
## 2. Install System Dependencies and Compilers
## --------------------------------------------------------------------------
# Added development headers (zlib1g-dev, libxml2-dev, etc.) necessary for R packages
# like 'multitaper', 'waveslim', and 'devtools' dependencies to compile on Linux.

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    gfortran \
    gcc \
    make \
    cmake \
    libblas-dev liblapack-dev \
    # Add critical development headers:
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

## --------------------------------------------------------------------------
## 3. Apply Custom R Compiler Configurations
## --------------------------------------------------------------------------
# Moved Makevars to the system-wide /etc/R/ to ensure it's respected by all users (like 'test').
# This prevents compiler optimization conflicts with RcppEigen (used by gradfps).

RUN mkdir -p /etc/R && \
    echo "CXXFLAGS = -O0 -fno-tree-vectorize -fno-optimize-sibling-calls -fPIC -std=c++14" >> /etc/R/Makevars && \
    echo "FC = /usr/bin/gfortran" >> /etc/R/Renviron

## --------------------------------------------------------------------------
## 4. Install R Packages (UPDATED and Finalized)
## --------------------------------------------------------------------------
# All CRAN dependencies are installed first.
RUN install2.r --error \
    # Dependencies required by LSPCA or your scripts:
    RSpectra \
    mvtnorm \
    Matrix \
    lpSolve \
    multitaper \
    waveslim \
    lattice \
    dplR \
    # General useful packages:
    ggplot2 \
    astsa \
    complexplus \
    devtools \
    && rm -rf /var/lib/apt/lists/*

# Install the critical GitHub package: gradfps (was commented out)
RUN R -e 'devtools::install_github("yixuan/gradfps", upgrade = "never", \
          lib = "/usr/local/lib/R/site-library", force = TRUE)'
RUN R -e 'devtools::install_github("jamnamdari/LSPCA", upgrade = "never", \
          lib = "/usr/local/lib/R/site-library", force = TRUE)'

# 5. Define the user/password as ARG for the build stage
ARG R_USER=test
ARG R_PASSWORD=1234

# 6. RUN command: Create the user and set the password using ARG values.
RUN useradd --create-home --shell /bin/bash $R_USER && \
    echo "$R_USER:$R_PASSWORD" | chpasswd

# 7. ENV command: Define RSTUDIO_USER environment variable.
ENV RSTUDIO_USER=$R_USER


# 8. Final configuration and execution
WORKDIR /home/test/project
EXPOSE 8787

# Final CMD
CMD ["/init"]
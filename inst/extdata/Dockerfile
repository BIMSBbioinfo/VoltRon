# Start from the RStudio base image
FROM rocker/rstudio:latest

# dependencies
RUN apt-get update
RUN apt-get install -y libgdal-dev libfftw3-dev libmagick++-dev cmake libhdf5-dev git libopencv-dev libopencv-features2d-dev  
RUN apt-get install -y libssl-dev libcurl4-openssl-dev libgit2-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libz-dev

# OpenCV
RUN git clone https://github.com/opencv/opencv.git
RUN git clone https://github.com/opencv/opencv_contrib.git
RUN mkdir build
WORKDIR "build"
RUN cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=/usr/local -D INSTALL_C_EXAMPLES=ON -D INSTALL_PYTHON_EXAMPLES=ON  -D OPENCV_GENERATE_PKGCONFIG=ON -D OPENCV_EXTRA_MODULES_PATH=../opencv_contrib/modules/  -D BUILD_opencv_xfeatures2d=ON ../opencv/
RUN make -j5
RUN make install
RUN sh -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/opencv.conf'
RUN ldconfig

# Install required R packages
RUN R -e "install.packages(c('shiny', 'devtools', 'BiocManager'), repos='http://cran.rstudio.com/')"

# Install VoltRon dependencies
RUN R -e "install.packages(c('grDevices', 'data.table', 'RcppAnnoy', 'RANN', 'Matrix', 'dplyr', 'ggplot2', 'ggrepel', 'igraph', 'irlba', 'rjson', 'magick', 'ids', 'sp', 'reshape2', 'rlang', 'ggpubr', 'shinyjs'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('stringr', 'uwot', 'RCDT', 'reticulate'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('EBImage', 'basilisk', 'S4Arrays'))"

# Install Suggested dependencies
RUN R -e "BiocManager::install(c('DelayedArray'))"
RUN R -e "BiocManager::install(c('HDF5Array'))"
RUN R -e "remotes::install_github('bnprks/BPCells/r@v0.3.0')"
RUN R -e "remotes::install_github('BIMSBbioinfo/ImageArray')"
RUN R -e "remotes::install_github('BIMSBbioinfo/HDF5DataFrame')"
RUN R -e "remotes::install_github('BIMSBbioinfo/ZarrDataFrame')"

# Install VoltRon dependencies
RUN R -e "devtools::install_github('BIMSBbioinfo/VoltRon')"

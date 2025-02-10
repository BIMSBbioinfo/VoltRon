# Start from the RStudio base image
FROM rocker/rstudio:latest

# dependencies
RUN apt-get update
RUN apt-get install -y libgdal-dev libfftw3-dev libmagick++-dev cmake libhdf5-dev git libopencv-dev libopencv-features2d-dev  
RUN apt-get install -y libssl-dev libcurl4-openssl-dev libgit2-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libz-dev

# OpenCV
RUN wget https://github.com/opencv/opencv/archive/refs/tags/4.8.1.zip
RUN unzip 4.8.1.zip
RUN rm 4.8.1.zip
RUN wget https://github.com/opencv/opencv_contrib/archive/refs/tags/4.8.1.zip
RUN unzip 4.8.1.zip
RUN rm 4.8.1.zip
RUN mkdir build
WORKDIR "build"
RUN cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=/usr/local -D INSTALL_C_EXAMPLES=ON -D INSTALL_PYTHON_EXAMPLES=ON  -D OPENCV_GENERATE_PKGCONFIG=ON -D OPENCV_EXTRA_MODULES_PATH=../opencv_contrib-4.8.1/modules/  -D BUILD_opencv_xfeatures2d=ON ../opencv-4.8.1/
RUN make -j5
RUN make install
RUN sh -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/opencv.conf'
RUN ldconfig

# Install required R packages
RUN R -e "install.packages(c('shiny', 'devtools', 'BiocManager'), repos='http://cran.rstudio.com/')"

# Install VoltRon dependencies
RUN R -e "install.packages(c('grDevices', 'data.table', 'RcppAnnoy', 'RANN', 'Matrix', 'dplyr', 'ggplot2', 'ggrepel', 'igraph', 'irlba', 'rjson', 'magick', 'ids', 'sp', 'reshape2', 'rlang', 'ggpubr', 'shinyjs'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('stringr', 'uwot', 'RCDT'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('EBImage', 'S4Arrays'))"

# Install Suggested dependencies
RUN R -e "BiocManager::install(c('DelayedArray'))"
RUN R -e "BiocManager::install(c('HDF5Array'))"
RUN R -e "remotes::install_github('bnprks/BPCells/r@v0.3.0')"
RUN R -e "remotes::install_github('BIMSBbioinfo/ImageArray')"
RUN R -e "remotes::install_github('BIMSBbioinfo/HDF5DataFrame')"
RUN R -e "remotes::install_github('BIMSBbioinfo/ZarrDataFrame')"
RUN R -e "install.packages('Seurat')"
RUN R -e "BiocManager::install('glmGamPoi')"
RUN R -e "install.packages('hdf5r')"
RUN R -e "install.packages('arrow')"
RUN R -e "BiocManager::install('RBioFormats')"
RUN R -e "BiocManager::install('ComplexHeatmap')"
RUN R -e "devtools::install_github('dmcable/spacexr')"
RUN R -e "devtools::install_github('xuranw/MuSiC')"
RUN R -e "BiocManager::install('SingleCellExperiment')"
RUN R -e "BiocManager::install('SpatialExperiment')"
RUN R -e "install.packages('dplyr')"
RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "install.packages('ggnewscale')"
RUN R -e "install.packages('patchwork')"
RUN R -e "install.packages('anndata')"

# Install VoltRon dependencies
RUN R -e "devtools::install_github('BIMSBbioinfo/VoltRon')"

# Install basilisk and setup environment
USER rstudio
RUN R -e "BiocManager::install('basilisk')"
RUN R -e "basilisk::obtainEnvironmentPath(VoltRon::getBasilisk())"

# return to root
USER root
RUN touch temp.txt

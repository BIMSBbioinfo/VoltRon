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
RUN R -e "install.packages(c('grDevices', 'data.table', 'RcppAnnoy', 'RANN', 'Matrix', 'dplyr', 'ggplot2', 'ggrepel', 'igraph', 'rjson', 'magick', 'ids', 'sp', 'reshape2', 'rlang', 'ggpubr', 'shinyjs'), repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('stla/RCDT')"
RUN R -e "install.packages(c('stringr', 'uwot'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('EBImage', 'S4Arrays', 'BiocSingular'))"

# set up java
USER root
RUN apt-get update -y
RUN apt upgrade -y
RUN apt-get install -y openjdk-21-jdk
RUN export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-arm64/
RUN R CMD javareconf -e
USER rstudio

# Install Suggested dependencies
RUN R -e "options(timeout = 600000000); remotes::install_github('BIMSBbioinfo/VoltRonStore')"
RUN R -e "options(timeout = 600000000); install.packages('Seurat')"
RUN R -e "BiocManager::install('glmGamPoi')"
RUN R -e "install.packages('arrow')"
RUN R -e "BiocManager::install('ComplexHeatmap')"
RUN R -e "options(timeout = 600000000); devtools::install_github('xuranw/MuSiC')"
RUN R -e "BiocManager::install('SingleCellExperiment')"
RUN R -e "BiocManager::install('SpatialExperiment')"
RUN R -e "install.packages('dplyr')"
RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "install.packages('ggnewscale')"
RUN R -e "install.packages('patchwork')"
RUN R -e "install.packages('anndata')"
RUN R -e "install.packages('R.utils')"
RUN R -e "devtools::install_github('immunogenomics/presto')"

# VoltRon
RUN R -e "devtools::install_github('BIMSBbioinfo/VoltRon')"

# Install basilisk and setup environment
USER rstudio
RUN R -e "BiocManager::install('basilisk')"
RUN R -e "basilisk::obtainEnvironmentPath(VoltRon::getBasilisk())"
RUN sh -c 'echo "options(voltron.python.path = \"/home/rstudio/.cache/R/basilisk/1.18.0/VoltRon/0.2.0/VoltRon_basilisk_env/bin/python\")" > /home/rstudio/.Rprofile'

# Install java based packages
USER root
RUN R -e "install.packages('rJava')"
RUN R -e "BiocManager::install('RBioFormats')"
RUN sh -c 'echo "options(java.parameters = \"-Xmx10g\")" >> /home/rstudio/.Rprofile'
USER rstudio
RUN R -e "options(timeout = 600000000); library(RBioFormats)"
USER root

# Install spacexr
RUN apt-get install -y libgsl-dev
# RUN R -e "options(timeout = 600000000); devtools::install_github(\"dmcable/spacexr\")"
RUN R -e "options(timeout = 600000000); BiocManager::install(\"spacexr\")"

# increase cache disk size for ImageMagick
RUN sed -i 's/2GiB/10GiB/g' /etc/ImageMagick-6/policy.xml

# vitessceR
RUN apt-get update -y
RUN apt upgrade -y
RUN apt-get install -y libsodium-dev 
RUN R -e "options(timeout = 600000000); devtools::install_github(\"vitessce/vitessceR\")"
RUN sh -c 'echo "options(timeout = 600000000)">> /home/rstudio/.Rprofile'

# SimpleITK
# TODO: for now dont install SimpleITK
# RUN R -e "devtools::install_github('SimpleITK/SimpleITKRInstaller', configure.vars=c('MAKEJ=1', 'ADDITIONAL_SITK_MODULES=-DSimpleITK_USE_ELASTIX=ON'))"

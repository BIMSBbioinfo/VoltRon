![](https://bimsbstatic.mdc-berlin.de/landthaler/VoltRon/Package/images/voltron_framework.png)

<br>

**Website and Tutorials**: <a href="https://bioinformatics.mdc-berlin.de/VoltRon">https://bioinformatics.mdc-berlin.de/VoltRon</a>

**VoltRon**  is a spatial omic analysis toolbox for multi-omics integration using spatial image registration. VoltRon is also capable of analyzing multiple types of spatially-aware data modalities.

   <ul class="maintext2">
     <li style="padding-bottom: 10px">
      <strong> Unique Data Structure </strong> of VoltRon allows users to seamlessly define tissue blocks, layers and multiple assay types in one R object.
     </li>
     <li style="padding-bottom: 10px">
      <strong> End-to-end downstream data analysis </strong> for distinct spatial biology technologies are supported. VoltRon visualizes and analyzes regions of interests (ROIs), spots, cells and even molecules <strong>(under development)</strong>.
     </li>
    <li style="padding-bottom: 10px">
      <strong> Automated Image Registration </strong> incorporates <a href="https://opencv.org/">OpenCV</a> (fully embedded into the package using <a href="https://www.rcpp.org/">Rcpp</a>) to detect common features across images and achieves registration. Users may interact with built-in mini shiny apps to change alignment parameters and validate alignment accuracy. 
    </li>
    <li style="padding-bottom: 10px">
      <strong> Manual Image Registration </strong> helps users to select common features across spatial datasets using reference images stored in VoltRon objects. In case automated image registration doesn't work, you can still align images by manually picking landmark points.
    </li>
    <li style="padding-bottom: 10px">
      <strong> Niche Clustering </strong> allows integration to single cell RNA analysis datasets using <a href="https://satijalab.org/seurat/">Seurat</a>, <a href="https://www.bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html">SingleCellExperiment</a> and <a href="https://github.com/dmcable/spacexr">spacexr</a> for spot deconvolution. Estimated cell type abundances are then used to cluster spots into groups of cell type niches which are defined as spots with distinct composition of cell types.  
    </li>
   </ul>

## Citation

VoltRon preprint is now available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.12.15.571667v1).


## Staying up-to-date

To ask questions please use VoltRon discussion forum on google groups.

- https://groups.google.com/forum/#!forum/voltron_discussion

<!--
## Configuration

Before installing VoltRon, you have to install [OpenCV](https://opencv.org/) library first.

On **Windows** and **MacOS**, OpenCV will be downloaded automatically upon installation. However, [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html) may be required to be downloaded too, hence this may take some time!

# To install on **MacOS**, you need to first install [homebrew](https://brew.sh/): 
# 
# ```sh
# /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# ```
# 
# Then, OpenCV is installed using **brew install** command.
# 
# ```sh
# brew install pkg-config
# brew install opencv
# ```
-->

<!--
On **Ubuntu** or **Fedora** you need [`libopencv-dev`](https://packages.debian.org/testing/libopencv-dev) or [`opencv-devel`](https://src.fedoraproject.org/rpms/opencv):

```sh
sudo apt-get install libopencv-dev
```
-->

## Installation

Install from the GitHub repository using devtools (with R version 4.3.0 or higher):

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("BIMSBbioinfo/VoltRon")
```

On **Windows** and **MacOS**, OpenCV will be downloaded automatically upon installation. However, [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html) may be required to be downloaded too, hence this may take some time!

On **Ubuntu** or **Fedora** you may need [`libopencv-dev`](https://packages.debian.org/testing/libopencv-dev) or [`opencv-devel`](https://src.fedoraproject.org/rpms/opencv):

```sh
sudo apt-get install libopencv-dev
```

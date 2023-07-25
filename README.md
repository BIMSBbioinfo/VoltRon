![](https://bimsbstatic.mdc-berlin.de/landthaler/VoltRon/Package/images/main_small.png)

<br>

**VoltRon**  is a spatial data analysis toolbox for spatial data analysis, multi-omics integration using spatial image registration. VoltRon is capable of analyzing multiple types and modalities of spatially-aware datasets.

   <ul class="maintext2">
     <li style="padding-bottom: 10px">
      <strong> Unique Data Structure </strong> of VoltRon allows users to seamlessly define tissue blocks, layers and multiple assay types.
     </li>
     <li style="padding-bottom: 10px">
      <strong> End-to-end Analysis </strong> for distinct spatial biology technologies are supported. VoltRon visualizes and analyzes regions of interests (ROIs), spots and cells **(even molecules, under development)**.
     </li>
    <li style="padding-bottom: 10px">
      <strong> Automated Spatial Registration </strong> incorporates <a href="https://opencv.org/">OpenCV</a> (fully embedded into the package using <a href="https://www.rcpp.org/">Rcpp</a>) to detect common features across images and achieves registration. Users may interact with built-in mini shiny apps to change alignment parameters and validate accuracy. 
    </li>
    <li style="padding-bottom: 10px">
      <strong> Manual Spatial Registration </strong> helps users to select common features across spatial datasets using reference images stored in VoltRon objects. In case, automated registration doesnt work, you can still align images.
    </li>
    <li style="padding-bottom: 10px">
      <strong> Niche Clustering </strong> allows integration to single cell RNA analysis platforms such as <a href="https://satijalab.org/seurat/">Seurat</a> and <a href="https://github.com/dmcable/spacexr">spacexr</a> to deconvolute spots using reference single cell datasets. Estimated cell type abundances are then used to partition spots into groups of cell type niches, defined as spots with distinct composition of cell types.  
    </li>
   </ul>

## Configuration

Before installing VoltRon, you have to install [OpenCV](https://opencv.org/) library first. 

To install on **MacOS**, you need to install the OpenCV library from homebrew:

```sh
brew install opencv
```

On **Ubuntu** or **Fedora** you need [`libopencv-dev`](https://packages.debian.org/testing/libopencv-dev) or [`opencv-devel`](https://src.fedoraproject.org/rpms/opencv):

```sh
sudo apt-get install libopencv-dev
```

On **Windows**, OpenCV will be downloaded automatically upon installation. 
However, [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html) may be required to be downloaded too, hence this may take some time!

## Installation

Install from the GitHub repository using devtools:

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("BIMSBbioinfo/VoltRon")
```

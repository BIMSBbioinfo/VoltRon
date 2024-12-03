![](https://bimsbstatic.mdc-berlin.de/landthaler/VoltRon/Package/images/voltron_framework_box_io.png)

<br>

**Website and Tutorials**: <a href="https://bioinformatics.mdc-berlin.de/VoltRon">https://bioinformatics.mdc-berlin.de/VoltRon</a>

**VoltRon**  is a spatial omic analysis toolbox for multi-omics integration using spatial image registration. VoltRon is also capable of analyzing multiple types of spatially-aware data modalities.
   
   <ul class="maintext2">
    <li style="padding-bottom: 10px">
      <strong> Unique data structure </strong> of VoltRon allows users to seamlessly define tissue blocks, layers and multiple assay types in one R object.
    </li>
    <li style="padding-bottom: 10px">
      <strong> End-to-end downstream data analysis </strong> for distinct spatial biology technologies are supported. VoltRon visualizes and analyzes regions of interests (ROIs), spots, cells, molecules and tiles **(under development)**.
    </li>
    <li style="padding-bottom: 10px">
      <strong> Automated Image Registration </strong> incorporates <a href="https://opencv.org/">OpenCV</a> (fully embedded into the package using <a href="https://www.rcpp.org/">Rcpp</a>) to detect common features across images and achieves registration. Users may interact with built-in mini shiny apps to change alignment parameters and validate alignment accuracy.
    </li>
    <li style="padding-bottom: 10px">
      <strong> Manual Image Registration </strong> helps users to select common features across spatial datasets using reference images stored in VoltRon objects. In case automated image registration doesn't work, you can still align images by manually picking landmark points.
    </li>
    <li style="padding-bottom: 10px">
    <p style="padding-bottom: 3px"> <strong> Spatially Aware Analysis </strong> allows detecting spatial patterns across cells, spots, molecules and other entities. </p>
    <ul class="maintext3">
      <li style="padding-bottom: 10px padding-top: 12px">
      <strong>(Niche Clustering: Spots)</strong> VoltRon allows integration to single cell RNA datasets using <a href="https://satijalab.org/seurat/">Seurat</a>, <a href="https://www.bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html">SingleCellExperiment</a> and <a href="https://github.com/dmcable/spacexr">spacexr</a> for spot deconvolution. Estimated cell type abundances are then used to cluster spots into groups of cell type niches which are defined as spots with distinct composition of cell types.
      </li>
      <li style="padding-bottom: 2px">
      <strong>(Niche Clustering: Cells)</strong> VoltRon creates spatial neighborhoods around cells to cluster local cellular compositions around all cells which in turn informs users on cell types that are likely within proximity to each other.
      </li>
      <li style="padding-bottom: 10px">
      <strong>(Hot Spot Detection)</strong> VoltRon detects region of locally spatial patterns of cells/molecules/spots that are abundant in biological events and/or features.
      </li>
    </ul>  
    </li>
    <li style="padding-bottom: 10px">
    <p> <strong> Support for Big Data </strong> for VoltRon objects enables storing large feature data matrices and large microscopic images of tissues on disk without overloading memory, thus allowing analysis on large datasets with ease. VoltRon stores large images as pyramid structures to speed up visualization and data retrieval. </p>
    </li>
    <li style="padding-bottom: 10px">
    <p> <strong> Interoperability across R/Python frameworks </strong> allows users to convert VoltRon objects to a large number of objects used by other spatial omic platforms such as Seurat, Squidpy (AnnData), SpatialExperiment (BioConductor) and Giotto. </p>
    </li>
  </ul>

## Staying up-to-date

To ask questions please use VoltRon discussion forum on google groups.

- https://groups.google.com/forum/#!forum/voltron_discussion

## Installation

Install from the GitHub repository using devtools (with R version 4.3.0 or higher):

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("BIMSBbioinfo/VoltRon")
```

Depending on the number of required dependencies, installation may be completed under a minute or may take a few minutes. 

On **Windows** and **MacOS**, OpenCV will be downloaded automatically upon installation. However, [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html) may be required to be downloaded too, hence this may take some time!

On **Ubuntu** we provide a set of instructions that may help users to build OpenCV with necessary headers [here](inst/extdata/install_ubuntu.md).

On **Fedora** you may need [`opencv-devel`](https://src.fedoraproject.org/rpms/opencv):

```sh
yum install opencv-devel
```

## Tutorials

Please see the [Explore](https://artur-man.github.io/VoltRon/tutorials.html) section in the VoltRon website for tutorials, example scripts and analysis found in the [preprint](https://www.biorxiv.org/content/10.1101/2023.12.15.571667v1). Tutorials include links for accessing necessary data to run scripts across all tutorials. 

## References

Manukyan, A., Bahry, E., Wyler, E., Becher, E., Pascual-Reguant, A., Plumbom, I., ... & Akalin, A. (2023). [VoltRon: A Spatial Omics Analysis Platform for Multi-Resolution and Multi-omics Integration using Image Registration](https://www.biorxiv.org/content/10.1101/2023.12.15.571667v1). bioRxiv, 2023-12.


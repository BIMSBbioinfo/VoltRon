# VoltRon 0.2.7

-   `transferData` now allows integrating single cell data object (`Seurat` or `SingleCellExperiment`)
    for transfering features (e.g. gene expression) or metadata features (e.g. cell types, annotations).
-   `registerSpatialData` now allows registering assays with no images. In case of either one of the assays, 
    do not have images, assays can be registered with only the manual approach. 
-   An image-free alignment tutorial has been added where DBIT-Seq and a QuPath processed mIF experiment are
    aligned using manually selected landmarks.
-   Now `importImageData` and `importQuPathIF` functions only work with segments already converted to a
    list by the user, or the `generateSegments` function whose arguement is an **sf** object. 
-   The `formAssay` function now accepts only segments without the user manually generating coordinates 
    (or centroids).
    
# VoltRon 0.2.6

-   Fix anndataR conversion by removing `var_names` and `obs_names` arguments and
    fix other faulty arguments

# VoltRon 0.2.5

-   Added instructions to install VoltRon from [R-universe](https://bimsbbioinfo.r-universe.dev/VoltRon).
-   Added instructions to install SimpleITK binaries from [BIMSBbioinfo/SimpleITKRInstaller](https://github.com/BIMSBbioinfo/SimpleITKRInstaller/releases).

# VoltRon 0.2.3

-   Non-rigid fine alignment with SimpleITK combined with landmark-based manual coarse alignment. Now, `registerSpatialData` function will let users to select **BSpline (SimpleITK)** method for non-rigid alignment when **Non-Rigid** methods selected.

# VoltRon 0.2.2

-   Non-rigid fine alignment (with OpenCV coarse alignment) using SimpleITK

    <https://bioinformatics.mdc-berlin.de/VoltRon/registration.html#Non-Rigid_Alignment>

-   Added 'channels' arguement to `ImportImageData` to parse only requested channels from OME-TIFF files.

    <https://bioinformatics.mdc-berlin.de/VoltRon/importingdata.html#OME-TIFF>

-   Multi-omic clustering with multiple feature sets and tutorials (Xenium + mIF) upon registration

    <https://bioinformatics.mdc-berlin.de/VoltRon/multiomic.html#Xenium_+_IF_Tonsil_Data>

-   Added `importQuPathIF` function to generate QuPath processed images and import spatial proteomics assays as a VoltRon objects

    [https://bioinformatics.mdc-berlin.de/VoltRon/importingdata.html#Multiplex_IF\_(QuPath)](https://bioinformatics.mdc-berlin.de/VoltRon/importingdata.html#Multiplex_IF_(QuPath)){.uri}

-   Building spatial neighborhood graphs and proximity analysis across multiple entities (e.g. molecules vs cells)

    <https://bioinformatics.mdc-berlin.de/VoltRon/multiomic.html#Hot_Spot_Analysis>

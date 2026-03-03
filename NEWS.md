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

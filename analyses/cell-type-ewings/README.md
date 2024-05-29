# Cell typing of Ewing sarcoma samples

## Description

This module will include code to annotate cell types in the Ewing sarcoma samples from SCPCP000015 present on the ScPCA Portal.

The first step in cell type annotation includes classifying cells as tumor or normal cells.
Currently, we use a variety of methods to annotate tumor cells independently:

- Identify cells that express high levels of tumor marker genes.
The full list of marker genes can be found in `references/tumor-marker-genes.tsv`.
- Use the list of tumor marker genes to classify tumor cells with [`CellAssign`](https://docs.scvi-tools.org/en/stable/user_guide/models/cellassign.html).
- Identify copy number variations and annotate tumor cells using [`CopyKAT`](https://github.com/navinlabcode/copykat).
- Identify copy number variations using [`InferCNV`](https://github.com/broadinstitute/inferCNV/wiki).
This returns a proportion of each chromosome with a CNV detected.
We then calculate the mean proportion for each cell across all chromsomes and weight by the number of genes in a chromosome.
Cells with a mean proportion > 0.1 are called as tumor cells.

## Sample metadata

The `sample-metadata.tsv` file is a TSV file containing information about samples used as input for workflows or analysis in this module.
This is manually updated as new samples are used in analysis.

Each row represents a unique sample and library combination.
The following columns are present in this file:

|  |   |
| --- | ---- |
|`scpca_sample_id`| Unique sample ID. The `sample_id` corresponds to the folder name containing data files for that sample after using `download-data.py`. |
|`scpca_library_id` | Unique library ID. The `library_id` will match the prefix of all data files (`.rds` and `.h5ad`) downloaded using `download-data.py`. |
|`normal_celltypes`| A comma separated list of cell types annotated with either `SingleR` or `CellAssign` used to create a reference list of normal cells |
|`tumor_celltypes`| A comma separated list of cell typs annotated with either `SingleR` or `CellAssign` that are expected to align with tumor cells. |

**Note:** To identify the cell type annotations to use for `normal_celltypes` and `tumor_celltypes`, reference the plots found in `<library_id>_celltype-report.html`.
These can be downloaded using the `--include_reports` option in `download-data.py`.

## Usage

To annotate tumor and normal cells in the Ewing's sarcoma samples from SCPCP000015, run the `annotate-tumor-cells-workflow.sh` workflow.
**Note:** Before running this workflow be sure to run `renv::restore()`  and activate the conda environment using the following commands:

```sh
Rscript -e "renv::restore()"
conda activate openscpca-cell-type-ewings
```

The following arguments are optional and can be used to run this workflow on additional samples (default sample is `SCPCS000490`):

- `sample_id`: Unique sample ID (name of folder containing libray data)
- `normal_celltypes`: Comma separated list of cell types annotated with either `SingleR` or `CellAssign` to use as a reference list of normal cells. This should correspond to the value found in `sample-metadata.tsv` for this sample.
- `tumor_celltypes`: Comma separated list of cell typs annotated with either `SingleR` or `CellAssign` that are expected to align with tumor cells.
Any cell types used here will be used for comparing to tumor cells annotated in this workflow.
This should correspond to the value found in `sample-metadata.tsv` for this sample.

Example of running the workflow with a different sample:

```sh
sample_id="SCPCS000491" ./annotate-tumor-cells.sh
```

## Input files

This module requires the processed `SingleCellExperiment` objects (`_processed.rds`) and processed `AnnData` objects (`.hdf5` files) from SCPCP0000015.
These files were obtained using the `download-data.py` script:

```sh
# download both SCE and AnnData objects
./download-data.py --projects SCPCP000015 --format "SCE,AnnData"
```
This module also requires the following reference files:

- `references/tumor-marker-genes.tsv`: This file contains a list of marker genes for identifying Ewing sarcoma tumor cells.
- `references/all-marker-genes.tsv`: This file contains a list of marker genes for identifying both normal and Ewing sarcoma tumor cells as described in [Visser, Beligs, _et al._ (2023)](https://doi.org/10.1158/2767-9764.CRC-23-0027).

## Output files

### Final output

Running the `annotate-tumor-cells.sh` will generate the following output files in `results/annotate_tumor_cells_output` for each sample/library combination.

```
annotate_tumor_cells_output
├── <sample_id>
    ├── <library_id>_cellassign-classifications.tsv
    ├── <library_id>_cellassign-report.html
    ├── <library_id>_copykat-classifications.tsv
    ├── <library_id>_copykat-report.html
    ├── <library_id>_infercnv-classifications.tsv
    ├── <library_id>_infercnv-report.html
    ├── <library_id>_marker-gene-report.html
    ├── <library_id>_tumor-normal-classifications.tsv
    ├── annotations
    │   └── <library_id>_reference-cells.tsv
    ├── cellassign
    │   ├── <library_id>_panglao-endo-fibro_predictions.tsv
    │   ├── <library_id>_tumor-marker_predictions.tsv
    │   └── <library_id>_visser-all-marker_predictions.tsv
    ├── copykat
    │   ├── no_reference
    │   │   ├── <library_id>_copykat_heatmap.jpeg
    │   │   └── <library_id>_final-copykat.rds
    │   └── with_reference
    │       ├── <library_id>_copykat_heatmap.jpeg
    │       └── <library_id>_final-copykat.rds
    └── infercnv
        ├── <library_id>_cnv-metadata.tsv
        ├── <library_id>_cnv-obj.rds
        └── <library_id>_infercnv.png
```

All `.html` files are rendered reports summarizing use of each method (indicated in the filename) to classify tumor cells.
All `classifications.tsv` file contain the final annotation as reported by each method.
The `annotations` folder contains the reference table indicating which cells were used as "normal" or "tumor" cells in various anlaysis.
See [below](#annotation-files) for more information on this table.

### Annotation files

Additionally, a reference TSV will be generated for each library containing a table of cell types from `SingleR` and `CellAssign`.
This file contains any cell barcodes expected to line up with tumor or normal cell types specified with the `normal_celltypes` and `tumor_celltypes` arguments when running `annotate-tumor-cells.sh`.
This table will be saved in `references/cell_lists/<sample_id>/<library_id>_reference-cells.tsv` and contains the following columns:

|  |   |
| -------- | -------- |
| `barcodes` | Cell barcode |
| `reference_cell_class` | Indicates if the cell should be uses as a Normal or Tumor cell reference |
| `cellassign_celltype_annotation` | Original annotation as obtained by `CellAssign` in the processed `SingleCellExperiment` object |
| `singler_celltype_annotation` | Original annotation as obtained by `SingleR` in the processed `SingleCellExperiment` object |

## Software requirements

All notebooks and R scripts are run using the `ewings-cell-type.Rproj`.
The `renv.lock` file contains all package and version information.

All Python scripts must be run using the conda environment `openscpca-cell-type-ewings`.
To create and activate this environment from the `.yml` file use:

```sh
conda env create -f environment.yml --name openscpca-cell-type-ewings
conda activate openscpca-cell-type-ewings
```

## Computational resources

Currently, the `annotate-tumor-cells.sh` can only be run using 1 CPU if running locally.
Use `threads=1` when running the workflow to specify this.

To increase the speed of the workflow, we recommend running on a computer with at least 12 CPUs.
If this is the case, you may run the workflow using the default settings of 4 threads.
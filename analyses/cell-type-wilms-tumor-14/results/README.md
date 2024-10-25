# Results directory instructions

Files in the results directory should not be directly committed to the repository.

Instead, copy results files to an S3 bucket and add a link to the S3 location in this README file.

- On a Lightsail virtual computer, run following script to sync the results:
```bash
export OPENSCPCA_RESULTS_BUCKET=researcher-009160072044-us-east-2
cd /home/lightsail-user/git/OpenScPCA-analysis
scripts/sync-results.py cell-type-wilms-tumor-14 \
    --bucket ${OPENSCPCA_RESULTS_BUCKET}
```
#### 00. Pre-processing the provided SCE objects
No result files.

#### 01. Anchor transfer using Seurat
* Results are uploaded to `s3://researcher-009160072044-us-east-2/cell-type-wilms-tumor-14/results/01_anchor_transfer_seurat`
```
.
├── [sample_id]_celltype.csv
├── [sample_id]_celltype.pdf
├── [sample_id]_compartment.csv
└── [sample_id]_compartment.pdf
```

* The label transfer analysis was performed in two levels: `celltype` and `compartment`.
* `[sample_id]_[level].csv` label transfer result table including cell ID, predicted cell type, along with predicted scores.
* `[sample_id]_[level].pdf` label transfer result plots consisting of 3 pages:
  1. UMAP visualization colored by transferred labels and Seurat clusters, as well as a bar plot showing cell type composition of each Seurat cluster.
  2. UMAP visualization colored and split by transferred labels.
  3. Distribution for max prediction score. Note: predictions with scores < 0.5 would be labeled as "Unknown" in this analysis.

#### 02. Curating marker gene lists
TBD

#### 03. Cell type annotation with marker gene lists
TBD

#### 04. Tumor cell identification
TBD

#### 05. Sample merging and validation
TBD
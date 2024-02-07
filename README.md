# Identification of Specialized tRNA Expression in Early Human Brain Development

<!-- [![DOI](https://zenodo.org/badge/DOI_GOES_HERE.svg)](https://zenodo.org/record/YOUR_ZENODO_RECORD_ID_HERE) -->

## About

<!-- **[PLACEHOLDER FOR PAPER CITATION - Replace with the actual citation once published]** -->

In this study, we investigated the expression of tRNA-derived small RNAs (tDRs) during early human cerebral cortex development. We utilized cerebral cortical organoid models and AlkB-facilitated RNA methylation sequencing (ARM-seq) to profile tDRs across distinct developmental stages. Our analysis revealed dynamic expression patterns of diverse tDR groups originating from a wide range of tRNA isodecoders, with several distinct groups exhibiting neural-specific expression. Computational analyses of these tDRs uncovered biased sequence motifs in over-represented tRNAs, which are enriched with particular RNA modifications. This provides initial insights into the traits that change within the pool of tDRs during neural development. This expanded catalog of tDRs provides a framework for future studies on tRNA function in brain development and enhances our understanding of the complexity of tDR dynamics in neural differentiation. Read alignment and tRNA/tDR quantification were performed using [tRAX](https://github.com/UCSC-LoweLab/tRAX.git), and [tRNAgraph](https://github.com/alba1735/tRNAgraph) was used for tRNA/tDR visualization. Custom Python code was developed for further analysis and figure generation.

## Repository Contents

This repository is organized as follows:

* **`trax_analysis.ipynb`**: Jupyter Notebook for running tRAX on the ARM-seq data. Includes specific parameters for raw fastq read trimming and tRAX execution.
* **`organoid_analysis.ipynb`**: Jupyter Notebook containing code snippets and visualizations for data preprocessing, statistical analysis, and visualization of the results presented in the paper, based on the tRAX output.
* **`config/`**: Directory containing configuration files for tRAX and tRNAgraph, including sample information and metadata.
* **`supplemental/`**: Directory containing supplementary data files referenced in the analysis notebooks.
* **`requirements.yaml`**: YAML file specifying the conda environment for reproducing the analysis.

### [tRAX Analysis](trax_analysis.ipynb)

This notebook performs the initial processing of the raw sequencing data using tRAX. The key steps are:

1. **Creating the hg38 RNAdb for tRAX:**  Generating the necessary tRNA database for read alignment.
2. **Trimming raw fastq reads:** Preprocessing the sequencing reads to remove adapter sequences and low-quality bases.
3. **Running tRAX:** Executing the tRAX pipeline for read alignment and quantification of tRNA and tDR expression.

### [Organoid Analysis](organoid_analysis.ipynb)

This notebook conducts the downstream analysis and generates the figures for the publication. The steps include:

1. **Building AnnData Objects:** Creating AnnData objects from the tRAX output for subsequent analysis.
2. **Manual Addition of Metadata:**
    * **Adding Neural Categories:** Incorporating neural cell type annotations based on prior clustering and analysis.
    * **Adding tRNA Modification Information:** Integrating tRNA modification data into the AnnData object.
    * **Adding tRNA Sequence Information:** Incorporates gtRNAdb (hg38) assocaited tRNA sequence information.
3. **Clustering and HDBSCAN Annotation:** Performing clustering and annotation of the ARM-seq data using HDBSCAN.
4. **Custom Visualizations:** Generating bespoke plots and figures to highlight key findings from the ARM-seq data.
5. **Bulk Visualizations:** Producing standard tRNAgraph visualizations to explore global tRNA and tDR patterns.
6. **Supplemental Data Analysis:** Integrating data from external studies for comparative analysis and validation.
7. **CSV Exportation:** Exporting processed data to CSV files for further use and inclusion as supplementary tables in the publication.

## Data Availability

All raw sequencing data generated in this study have been submitted to the [NCBI Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) under accession number [GSE259250](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE259250).

## Installation and Running the Analysis

To reproduce the analysis, you will need to install the required software environment and download the necessary tools.

**1. Install the conda environment:**

The analysis environment can be created using `conda` or `mamba` with the `requirements.yaml` file located in the root of this repository.

```bash
mamba env create -f requirements.yaml
```

**2. Install tRAX and tRNAgraph:**

In addition to the conda dependencies, tRAX and tRNAgraph are used in this analysis. These should be downloaded into directories at the same level as the current repository (i.e., sister directories). For reproducibility, this analysis utilizes the `py3test` branch of the tRAX repository at tag `v1.1.0-beta`. tRNAgraph is used on the `main` branch at tag `v1.0.0`.

```bash
# Navigate to the parent directory of the current repository
cd ..

# Clone tRAX on the py3test branch with the specified tag
git clone --depth 1 --branch v1.1.0-beta https://github.com/UCSC-LoweLab/tRAX.git tRAX
cd tRAX
git reset --hard 17d6bfeb2e383348556af26813e6fced8faa79f3
cd ..

# Clone tRNAgraph
git clone --depth 1 --branch v1.0.0 https://github.com/alba1735/tRNAgraph.git tRNAgraph
cd tRNAgraph
git reset --hard da5dfb7089136bc189bf3f67172fe6ad0766e8a0
cd ..

# Navigate back to the repository's root directory
cd organoid_analysis # Assuming this README is in the root, otherwise adjust the path
```

**3. Run the Jupyter Notebooks:**

The analysis can be executed by running the Jupyter notebooks.

```bash
jupyter lab
```

This command will open Jupyter Lab in your web browser. Navigate to and open `trax_analysis.ipynb` and `organoid_analysis.ipynb` to execute the analysis steps.

## tRAX Configuration

The tRAX analysis utilizes the following files for its configuration:

* `config/runfile.tsv`: Defines the input raw fastq files and the desired output directories for tRAX.
* `config/samples.arm.tsv`: Lists the samples, their corresponding groups, and the location of the merged fastq files.
* `config/pairs.arm.tsv`: Specifies the pairs of conditions for differential expression analysis within tRAX.

The contents of these files are provided below:

`config/runfile.tsv`:

```tsv
rnaseq/fastq_processed/AlkBm_arm_d0_1 rnaseq/fastq_raw/AlkBm_arm_d0_1_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d0_1_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d0_2 rnaseq/fastq_raw/AlkBm_arm_d0_2_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d0_2_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d0_3 rnaseq/fastq_raw/AlkBm_arm_d0_3_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d0_3_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d14_1 rnaseq/fastq_raw/AlkBm_arm_d14_1_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d14_1_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d14_2 rnaseq/fastq_raw/AlkBm_arm_d14_2_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d14_2_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d35_1 rnaseq/fastq_raw/AlkBm_arm_d35_1_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d35_1_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d35_2 rnaseq/fastq_raw/AlkBm_arm_d35_2_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d35_2_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d35_3 rnaseq/fastq_raw/AlkBm_arm_d35_3_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d35_3_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d35_4 rnaseq/fastq_raw/AlkBm_arm_d35_4_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d35_4_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d35_5 rnaseq/fastq_raw/AlkBm_arm_d35_5_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d35_5_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d35_6 rnaseq/fastq_raw/AlkBm_arm_d35_6_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d35_6_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d70_1 rnaseq/fastq_raw/AlkBm_arm_d70_1_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d70_1_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d70_2 rnaseq/fastq_raw/AlkBm_arm_d70_2_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d70_2_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBm_arm_d70_3 rnaseq/fastq_raw/AlkBm_arm_d70_3_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBm_arm_d70_3_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d0_1 rnaseq/fastq_raw/AlkBp_arm_d0_1_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d0_1_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d0_2 rnaseq/fastq_raw/AlkBp_arm_d0_2_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d0_2_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d0_3 rnaseq/fastq_raw/AlkBp_arm_d0_3_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d0_3_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d0_4 rnaseq/fastq_raw/AlkBp_arm_d0_4_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d0_4_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d0_5 rnaseq/fastq_raw/AlkBp_arm_d0_5_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d0_5_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d0_6 rnaseq/fastq_raw/AlkBp_arm_d0_6_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d0_6_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d0_7 rnaseq/fastq_raw/AlkBp_arm_d0_7_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d0_7_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d0_8 rnaseq/fastq_raw/AlkBp_arm_d0_8_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d0_8_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d14_1 rnaseq/fastq_raw/AlkBp_arm_d14_1_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d14_1_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d14_2 rnaseq/fastq_raw/AlkBp_arm_d14_2_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d14_2_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d35_1 rnaseq/fastq_raw/AlkBp_arm_d35_1_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d35_1_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d35_2 rnaseq/fastq_raw/AlkBp_arm_d35_2_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d35_2_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d35_3 rnaseq/fastq_raw/AlkBp_arm_d35_3_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d35_3_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d35_4 rnaseq/fastq_raw/AlkBp_arm_d35_4_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d35_4_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d35_5 rnaseq/fastq_raw/AlkBp_arm_d35_5_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d35_5_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d70_1 rnaseq/fastq_raw/AlkBp_arm_d70_1_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d70_1_L001_R2.fastq.gz
rnaseq/fastq_processed/AlkBp_arm_d70_2 rnaseq/fastq_raw/AlkBp_arm_d70_2_L001_R1.fastq.gz rnaseq/fastq_raw/AlkBp_arm_d70_2_L001_R2.fastq.gz
```

`config/samples.arm.tsv`:

```tsv
AlkBm_arm_d0_1	AlkBm_arm_d0	rnaseq/fastq_processed/AlkBm_arm_d0_1_merge.fastq.gz
AlkBm_arm_d0_2	AlkBm_arm_d0	rnaseq/fastq_processed/AlkBm_arm_d0_2_merge.fastq.gz
AlkBm_arm_d0_3	AlkBm_arm_d0	rnaseq/fastq_processed/AlkBm_arm_d0_3_merge.fastq.gz
AlkBm_arm_d14_1	AlkBm_arm_d14	rnaseq/fastq_processed/AlkBm_arm_d14_1_merge.fastq.gz
AlkBm_arm_d14_2	AlkBm_arm_d14	rnaseq/fastq_processed/AlkBm_arm_d14_2_merge.fastq.gz
AlkBm_arm_d35_1	AlkBm_arm_d35	rnaseq/fastq_processed/AlkBm_arm_d35_1_merge.fastq.gz
AlkBm_arm_d35_2	AlkBm_arm_d35	rnaseq/fastq_processed/AlkBm_arm_d35_2_merge.fastq.gz
AlkBm_arm_d35_3	AlkBm_arm_d35	rnaseq/fastq_processed/AlkBm_arm_d35_3_merge.fastq.gz
AlkBm_arm_d35_4	AlkBm_arm_d35	rnaseq/fastq_processed/AlkBm_arm_d35_4_merge.fastq.gz
AlkBm_arm_d35_5	AlkBm_arm_d35	rnaseq/fastq_processed/AlkBm_arm_d35_5_merge.fastq.gz
AlkBm_arm_d35_6	AlkBm_arm_d35	rnaseq/fastq_processed/AlkBm_arm_d35_6_merge.fastq.gz
AlkBm_arm_d70_1	AlkBm_arm_d70	rnaseq/fastq_processed/AlkBm_arm_d70_1_merge.fastq.gz
AlkBm_arm_d70_2	AlkBm_arm_d70	rnaseq/fastq_processed/AlkBm_arm_d70_2_merge.fastq.gz
AlkBm_arm_d70_3	AlkBm_arm_d70	rnaseq/fastq_processed/AlkBm_arm_d70_3_merge.fastq.gz
AlkBp_arm_d0_1	AlkBp_arm_d0	rnaseq/fastq_processed/AlkBp_arm_d0_1_merge.fastq.gz
AlkBp_arm_d0_2	AlkBp_arm_d0	rnaseq/fastq_processed/AlkBp_arm_d0_2_merge.fastq.gz
AlkBp_arm_d0_3	AlkBp_arm_d0	rnaseq/fastq_processed/AlkBp_arm_d0_3_merge.fastq.gz
AlkBp_arm_d0_4	AlkBp_arm_d0	rnaseq/fastq_processed/AlkBp_arm_d0_4_merge.fastq.gz
AlkBp_arm_d0_5	AlkBp_arm_d0	rnaseq/fastq_processed/AlkBp_arm_d0_5_merge.fastq.gz
AlkBp_arm_d0_6	AlkBp_arm_d0	rnaseq/fastq_processed/AlkBp_arm_d0_6_merge.fastq.gz
AlkBp_arm_d0_7	AlkBp_arm_d0	rnaseq/fastq_processed/AlkBp_arm_d0_7_merge.fastq.gz
AlkBp_arm_d0_8	AlkBp_arm_d0	rnaseq/fastq_processed/AlkBp_arm_d0_8_merge.fastq.gz
AlkBp_arm_d14_1	AlkBp_arm_d14	rnaseq/fastq_processed/AlkBp_arm_d14_1_merge.fastq.gz
AlkBp_arm_d14_2	AlkBp_arm_d14	rnaseq/fastq_processed/AlkBp_arm_d14_2_merge.fastq.gz
AlkBp_arm_d35_1	AlkBp_arm_d35	rnaseq/fastq_processed/AlkBp_arm_d35_1_merge.fastq.gz
AlkBp_arm_d35_2	AlkBp_arm_d35	rnaseq/fastq_processed/AlkBp_arm_d35_2_merge.fastq.gz
AlkBp_arm_d35_3	AlkBp_arm_d35	rnaseq/fastq_processed/AlkBp_arm_d35_3_merge.fastq.gz
AlkBp_arm_d35_4	AlkBp_arm_d35	rnaseq/fastq_processed/AlkBp_arm_d35_4_merge.fastq.gz
AlkBp_arm_d35_5	AlkBp_arm_d35	rnaseq/fastq_processed/AlkBp_arm_d35_5_merge.fastq.gz
AlkBp_arm_d70_1	AlkBp_arm_d70	rnaseq/fastq_processed/AlkBp_arm_d70_1_merge.fastq.gz
AlkBp_arm_d70_2	AlkBp_arm_d70	rnaseq/fastq_processed/AlkBp_arm_d70_2_merge.fastq.gz
```

`config/pairs.arm.tsv`:

```tsv
AlkBm_arm_d0	AlkBp_arm_d0
AlkBm_arm_d14	AlkBp_arm_d14
AlkBm_arm_d35	AlkBp_arm_d35
AlkBm_arm_d70	AlkBp_arm_d70
AlkBp_arm_d0	AlkBp_arm_d14
AlkBp_arm_d0	AlkBp_arm_d35
AlkBp_arm_d0	AlkBp_arm_d70
AlkBp_arm_d14	AlkBp_arm_d35
AlkBp_arm_d14	AlkBp_arm_d70
AlkBp_arm_d35	AlkBp_arm_d70
```

## tRNAgraph Configuration

tRNAgraph is used to generate bulk visualizations. The AnnData objects for tRNAgraph are built using the following metadata and observations:

`config/metadata.tsv`:

```tsv
sample group celltype timepoint samplenamefancy sequencingtype treatment tissuetype species
AlkBm_arm_d0_1 AlkBm_arm_d0 c305 Day 0 d0_1 ARM-seq AlkBm Cell hg38
AlkBm_arm_d0_2 AlkBm_arm_d0 c305 Day 0 d0_2 ARM-seq AlkBm Cell hg38
AlkBm_arm_d0_3 AlkBm_arm_d0 c305 Day 0 d0_3 ARM-seq AlkBm Cell hg38
AlkBm_arm_d14_1 AlkBm_arm_d14 c305 Day 14 d14_1 ARM-seq AlkBm Organoid hg38
AlkBm_arm_d14_2 AlkBm_arm_d14 c305 Day 14 d14_2 ARM-seq AlkBm Organoid hg38
AlkBm_arm_d35_1 AlkBm_arm_d35 c305 Day 35 d35_1 ARM-seq AlkBm Organoid hg38
AlkBm_arm_d35_2 AlkBm_arm_d35 c305 Day 35 d35_2 ARM-seq AlkBm Organoid hg38
AlkBm_arm_d35_3 AlkBm_arm_d35 c305 Day 35 d35_3 ARM-seq AlkBm Organoid hg38
AlkBm_arm_d35_4 AlkBm_arm_d35 c305 Day 35 d35_4 ARM-seq AlkBm Organoid hg38
AlkBm_arm_d35_5 AlkBm_arm_d35 c305 Day 35 d35_5 ARM-seq AlkBm Organoid hg38
AlkBm_arm_d35_6 AlkBm_arm_d35 c305 Day 35 d35_6 ARM-seq AlkBm Organoid hg38
AlkBm_arm_d70_1 AlkBm_arm_d70 c305 Day 70 d70_1 ARM-seq AlkBm Organoid hg38
AlkBm_arm_d70_2 AlkBm_arm_d70 c305 Day 70 d70_2 ARM-seq AlkBm Organoid hg38
AlkBm_arm_d70_3 AlkBm_arm_d70 c305 Day 70 d70_3 ARM-seq AlkBm Organoid hg38
AlkBp_arm_d0_1 AlkBp_arm_d0 c305 Day 0 d0_1 ARM-seq AlkBp Cell hg38
AlkBp_arm_d0_2 AlkBp_arm_d0 c305 Day 0 d0_2 ARM-seq AlkBp Cell hg38
AlkBp_arm_d0_3 AlkBp_arm_d0 c305 Day 0 d0_3 ARM-seq AlkBp Cell hg38
AlkBp_arm_d0_4 AlkBp_arm_d0 c305 Day 0 d0_4 ARM-seq AlkBp Cell hg38
AlkBp_arm_d0_5 AlkBp_arm_d0 c305 Day 0 d0_5 ARM-seq AlkBp Cell hg38
AlkBp_arm_d0_6 AlkBp_arm_d0 c305 Day 0 d0_6 ARM-seq AlkBp Cell hg38
AlkBp_arm_d0_7 AlkBp_arm_d0 c305 Day 0 d0_7 ARM-seq AlkBp Cell hg38
AlkBp_arm_d0_8 AlkBp_arm_d0 c305 Day 0 d0_8 ARM-seq AlkBp Cell hg38
AlkBp_arm_d14_1 AlkBp_arm_d14 c305 Day 14 d14_1 ARM-seq AlkBp Organoid hg38
AlkBp_arm_d14_2 AlkBp_arm_d14 c305 Day 14 d14_2 ARM-seq AlkBp Organoid hg38
AlkBp_arm_d35_1 AlkBp_arm_d35 c305 Day 35 d35_1 ARM-seq AlkBp Organoid hg38
AlkBp_arm_d35_2 AlkBp_arm_d35 c305 Day 35 d35_2 ARM-seq AlkBp Organoid hg38
AlkBp_arm_d35_3 AlkBp_arm_d35 c305 Day 35 d35_3 ARM-seq AlkBp Organoid hg38
AlkBp_arm_d35_4 AlkBp_arm_d35 c305 Day 35 d35_4 ARM-seq AlkBp Organoid hg38
AlkBp_arm_d35_5 AlkBp_arm_d35 c305 Day 35 d35_5 ARM-seq AlkBp Organoid hg38
AlkBp_arm_d70_1 AlkBp_arm_d70 c305 Day 70 d70_1 ARM-seq AlkBp Organoid hg38
AlkBp_arm_d70_2 AlkBp_arm_d70 c305 Day 70 d70_2 ARM-seq AlkBp Organoid hg38
```

For certain figures, the data was subsetted to include only AlkB Positive samples using the following configuration:

`config/AlkBp.json`:

```json
{
    "name": "AlkBp",
    "obs": {
        "treatment": ["AlkBp"],
        "celltype": ["c305"],
        "pseudogene": ["tRNA"]
    },
    "obs_r": {
        "amino": ["Und"] 
    }
}
```

The following color map was used for visualizations:

`config/colormap.json`:

```json
{
    "timepoint": {
        "Day 0": "#007FFF",
        "Day 14": "#00FF7F",
        "Day 35": "#FF007F",
        "Day 70": "#FFD700"
    },
    "group": {
        "AlkBm_arm_d0": "#007FFF",
        "AlkBm_arm_d14": "#00FF7F",
        "AlkBm_arm_d35": "#FF007F",
        "AlkBm_arm_d70": "#FFD700",
        "AlkBp_arm_d0": "#007FFF",
        "AlkBp_arm_d14": "#00FF7F",
        "AlkBp_arm_d35": "#FF007F",
        "AlkBp_arm_d70": "#FFD700"
    },
    "amino": {
        "Ala": "#1F77B4",
        "Arg": "#AEC7E8",
        "Asn": "#FF7F0E",
        "Asp": "#FFBB78",
        "Cys": "#6fe835",
        "Gln": "#d0f2cb",
        "Glu": "#D62728",
        "Gly": "#FF9896",
        "His": "#9258f5",
        "Ile": "#deccfc",
        "Leu": "#a65223",
        "Lys": "#ffceb3",
        "iMet": "#00d5e3",
        "Met": "#b8fbff",
        "Phe": "#edd500",
        "Pro": "#ffff99",
        "Ser": "#db56bc",
        "Thr": "#F7B6D2",
        "Trp": "#2CA02C",
        "Tyr": "#98DF8A",
        "Val": "#6a3d9a",
        "SeC": "#C5B0D5",
        "Sup": "#808080",
        "Und": "#000000"
    },
    "neural_fine": {
        "strong neural": "darkorange",
        "neural": "gold",
        "weak neural": "khaki",
        "neutral": "mediumspringgreen",
        "weak stem": "lightskyblue",
        "stem": "deepskyblue",
        "strong stem": "royalblue",
        "undetermined": "darkgray"
    },
    "neural_broad": {
        "neural": "gold",
        "neutral": "mediumspringgreen",
        "stem": "deepskyblue",
        "undetermined": "darkgray"
    },
    "sequencingtype": {
        "ARM-seq": "#9718db"
    }
}
```

## Supplementary Data

The `supplemental/` directory contains files used to enrich and validate the primary analysis. Each file is described below, including its source and how it was used in the study:

* **`ModTable_Human_labeled.txt`**: This tab-separated file contains information on known tRNA modifications in humans. It was used to add modification annotations to the AnnData objects, providing a basis for linking tRNA modifications to the observed tDR expression patterns. This data was meticulously curated and re-formatted from the information presented in:
    * Zhang et al. 2025. Human TRMT1 and TRMT1L paralogs ensure the proper modification state, stability, and function of tRNAs. *Cell Reports*, 44(1):115092. [DOI: 10.1016/j.celrep.2024.115092](https://doi.org/10.1016/j.celrep.2024.115092). The data was converted into a tab-separated format for compatibility with the analysis pipeline.

* **`mods_disease.tsv`**: This tab-separated file contains information on tRNA modifications and their association with diseases. It was used for comparative analysis of modification patterns and tDR fragment characteristics observed in our data with those linked to disease states. This file was downloaded and converted to a tab-separated format from the supplementary material of:
    * Crécy-Lagard et al. 2019. Matching tRNA modifications in humans to their known and predicted enzymes. *Nucleic Acids Research*, 47(5):2143–2159. [DOI: 10.1093/nar/gkz011](https://doi.org/10.1093/nar/gkz011).

* **`1-s2.0-S0092867418303830-mmc2.csv`**: This comma-separated value (CSV) file contains gene expression data from human cerebral organoids generated by Fiddes et al. (2018). It was used to validate gene expression patterns observed in our organoid data, providing an external point of reference. This file was directly downloaded from the supplementary data of:
    * Fiddes et al. 2018. Human-Specific NOTCH2NL Genes Affect Notch Signaling and Cortical Neurogenesis. *Cell*, 173(6):1356-1369.e22. [DOI: 10.1016/j.cell.2018.03.051](https://doi.org/10.1016/j.cell.2018.03.051).

* **`hg38_small_ncRNAs.gtf`**: This GTF (Gene Transfer Format) file provides annotations for various small non-coding RNAs (snRNAs) in the human genome assembly GRCh38/hg38. It was used to add annotations for small ncRNAs, allowing us to differentiate tDRs from other small RNA species. These annotations were compiled by Patricia Chan of the Lowe Lab.
    * [Patricia Chan's Google Scholar](https://scholar.google.com/citations?user=7FWTtUIAAAAJ&hl=en)
    * [Patricia Chan's GitHub](https://github.com/patriciaplchan)

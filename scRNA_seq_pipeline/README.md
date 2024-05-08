## Overall Goals of this Pipeline  
The goal of this pipeline is to use CellRanger and Seurat to automatically label
cell types from single cell sequencing experiments, and to compare the
expression of genes of a given cell type in one experimental condition against
the expression of the same set of genes in the same cell type in a second
condition.

Because this is a snakemake program, you can modify parameters (e.g. cell
clustering, cell type labels, or 'real' cell filtering parameters) and only the
subset of steps that rely on these parameters will be rerun. For example, if you
run this program once, it might take an hour to run the 'cellranger' step, plus
a couple minutes to label cell types. If you change the maximum allowable
mitochondrial gene expression parameter in the config file, the program will
rerun in a couple minutes (skipping the cellranger step and only modifying cell
type labels), because mitochondrial gene expression thresholds don't affect
cellranger.

Advanced users can also specify desired intermediate output files at the
commandline, to test the effects of altered parameters and avoid re-generating
final output files every time. See 'advanced usage' section for more details.


## Inputs:
- Input a folder that has multiple samples, with each sample having its own folder of fastq files
- Input the names of the samples to compare
- Modify the file "markers.csv" so that it has marker genes appropriate to your cell types of interest
- Input the location of an annotated genome folder. You can obtain a correctly formatted example genome from 10X genomics with the code below:
```bash
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
```
- Input a name of a folder where the outputs will go for your experiment
- Input filtering parameters for cells that are considered 'real' cells (e.g. minimum number of expressed genes, maximum number of expressed genes, maximum percent mitochondrial gene expression)
- genes of interest whose expression you'd like to see in the dataset before clustering cell types in Seurat

## Outputs:
- violin plots showing the distribution of number of genes, expression level of genes, and percent mitochondrial gene expression in each cell
- violin plots showing expression of genes of interest pre-clustering
- violin plots showing expression of genes of interest in each cluster post-clustering
- elbow plot showing amount of variance explained by each principal component
- scatter plots showing gene expression levels vs. number of genes expressed and vs. mitochondrial concentration in each sample
- rds files for export into other projects
- "actual_cluster_labels" used to define cells of interest based on marker genes of interest - user can edit this file to redefine cluster labels as needed
- a labeled umap assigning a cell type to each cluster
- differential expression of various genes within a partiular cell type between experimental and control conditions


## Installation from terminal of necessary tools  
Prerequisites: This tutorial assumes a unix-like operating system (e.g. Mac or Linux) but with small modifications it may also work in Windows.


- install cellranger

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

- install mambaforge or a comparable version of conda

https://mamba.readthedocs.io/en/latest/installation.html

- create a conda environment for this pipeline with a name you can remember
```bash
conda create -n nameofmyenvironment
```
- activate the environment (you can deactivate at any time with "mamba deactivate")
```bash
conda activate nameofmyenvironment
```
- install snakemake into the environment
```bash
conda install -c conda-forge -c bioconda snakemake
```
- install seurat into the environment
```bash
conda install -c conda-forge r-seurat
```




## Using the pipeline
- Copy the pipeline folder into the directory where you want the outputs to show up
- Modify the markers.csv file in your copy to suit your needs
- Open your copy of the config file
- Modify the directories for fastq_folder, reference_genome, and gene_cluster_markers - make sure to put a terminating '/' at the end of folders
- Modify the samples to match the pair of samples in your fastq_folder that you will be comparing.  Each run will be comparing a pair of samples to eachother.   The first sample is treated as the control and the second sample is experimental
- Run the pipeline by activating your environment you created during the installation phase, changing directory into the folder that contains the Snakefile, then typing the command below (you can choose to use more cores or less if desired).  Cellranger will take a while to run if this is the first run for the chosen pair.  Subsequent runs will be much faster.
```bash
snakemake --cores 4
```
- Review the files in the 'intermediate' folder to see if you want to adjust any QC parameters or marker features in the config file.  Rerun the 'snakemake --cores 4' command if you make any adjustments
- Review the 'actual_cluster_labels' file in the intermediate folder and edit it to make sure the clusters are named appropriately.  Rerun the 'snakemake --cores 4' command when you are finished
- Change the cell_type_of_interest in the config file to match any of the cluster names you finalized while renaming the actual_cluster_labels, then rerun the 'snakemake --cores 4' command


## Advanced Options and things to watch out for
- If you change parameters that affect clustering, such as the QC parameters, resolution, or principal_components, the actual_cluster_names will be re-generated automatically (to match the newly formed clusters) and you will need to manually edit this file to make sure clusters are named appropriately.
- If you wish to compare results with different parameters you can adjust the run_number in the config file.  This will create a new folder labeled with your chosen run_number and put all outputs there.  The old folder will not be touched unless you change the run number back again
- Principal components and resolution can be adjusted at the bottom of the config file, but most users will probably want to leave them as is
- You can re-run the pipeline only up to a certain point (without asking the
pipeline to run to completion). For example, if a successfully executed run
compares sample T05 against T01, and sends output to a folder STRESS_T05_vs_T00_run_1/intermediate,
and a user wants to see the effects of increasing cluster resolution to 0.7 on
the cluster labels of the file STRESS_T05_vs_T00_run_1/intermediate/actual_cluster_labels_T05_vs_T00_run_1.txt
a user could change the config file to have a resolution of 0.7, and then run

```bash
snakemake STRESS_T05_vs_T00_run_1/intermediate/actual_cluster_labels_T05_vs_T00_run_1.txt --cores 4  
```

to rerun only the portion of the pipeline that generates this intermediate file.

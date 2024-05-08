FROM satijalab/seurat
COPY automated_Merfish_Seurat_part1.R /
COPY automated_Merfish_Seurat_part2.R /
COPY automated_Merfish_Seurat_driver.py /
COPY get_cluster_markers.py /
COPY Mouse_Brain_M22_210405_merged_metadata_and_partitions.csv /
COPY cell_type_labels.tsv /
COPY custom2.py /
COPY micro /bin
RUN Rscript automated_Merfish_Seurat_part1.R

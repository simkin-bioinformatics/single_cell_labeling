import subprocess

subprocess.call(['Rscript', 'automated_Merfish_Seurat_part1.R'])
subprocess.call(['python3', 'get_cluster_markers.py'])
junk=input('go ahead and push the Enter button when finished editing '
'actual_cluster_labels.txt')
subprocess.call(['Rscript', 'automated_Merfish_Seurat_part2.R'])

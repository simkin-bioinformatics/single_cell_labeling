'''
takes the marker genes that define clusters from seurat as input.
Goal is to get markers that uniquely identify each cluster from seurat
'''

cell_type_labels = snakemake.input[0]
markers = snakemake.input[1]
actual_cluster_labels = snakemake.output[0]

marker_dict=dict([line.strip().split(',') for line in open(cell_type_labels)])
#print(marker_dict)
import custom2
cluster_dict={}
adjusted_p_c, cluster_c, gene_c=5,6,7
for line_number, line in enumerate(open(markers)):
	if line_number>0:
		line=line.strip().split()
		adjusted_p=float(line[adjusted_p_c])
		cluster=line[cluster_c]
		#print(line, adjusted_p)
		if adjusted_p<0.05:
			cluster=int(line[cluster_c])
			gene=line[gene_c]
			cluster_dict[cluster]=custom2.get(cluster_dict, cluster, set([]))|set([gene])
#print(cluster_dict)
new_cluster_dict={}
cluster_list=sorted(list(cluster_dict.keys()))
for cluster in cluster_list:
	gene_set=cluster_dict[cluster]
	temp_list=cluster_list[:]
	temp_list.pop(cluster)
	removal_set=set([])
	for cluster2 in temp_list:
		removal_set=removal_set|cluster_dict[cluster2]
	new_cluster_dict[cluster]=gene_set-removal_set
#print(new_cluster_dict)

final_dict={}
for cluster in cluster_list:
	gene_set=cluster_dict[cluster]
	for marker in marker_dict:
		if marker in gene_set:
			final_dict[cluster]=custom2.get(final_dict, cluster, [])+[marker_dict[marker]]
	if cluster not in final_dict:
		final_dict[cluster]=['no_markers']
	final_dict[cluster]=sorted(list(set(final_dict[cluster])))
	print(cluster, final_dict[cluster])

final_list=[]
for cluster in cluster_list:
	final_list.append('/'.join(final_dict[cluster]))
output_file=open(actual_cluster_labels, 'w')
output_file.write('cluster_numbers\tcluster_names\n')
for cluster_number, cluster_name in enumerate(final_list):
	output_file.write(str(cluster_number)+'\t'+cluster_name+'\n')
output_file.close()

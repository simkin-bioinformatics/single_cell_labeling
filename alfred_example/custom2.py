'''
a library of modules I use frequently in computational biology. Unlike custom,
this version is modified for use with python3
'''

def count_in_base(string, base, valuestring):
	'''
	starts with a starting string. User specifies what 'base' to count in (e.g.
	base 10 for iterating through regular numbers), the values of valuestring
	are used to iterate through all possible values.
	'''
	outlist=[]
	full_length=len(string)
	for count in range(base**len(string)):
		outlist.append(string)
		if string==valuestring[-2]*full_length:
			break
		position=len(string)-1
		value=valuestring.index(string[position])
		string=string[:position]+valuestring[value+1]+string[position+1:]
		while value>=(base-1) and count<(base**len(string)):
			string=string[:position]+valuestring[0]+string[position+1:]
			position-=1
			value=valuestring.index(string[position])
			string=string[:position]+valuestring[value+1]+string[position+1:]
	return outlist

#def revcom(sequence, sequence_type='DNA'):
#	'''
#	the guts of this module are not of my creation
#	I've added some tweaks
#	'''
#	import string
#	if sequence_type=='RNA':
#		complement = string.maketrans('ACGTUN', 'UGCAAN')
#	elif sequence_type=='DNA':
#		complement = string.maketrans('ACGTUN', 'TGCAAN')
#	return sequence.upper().translate(complement)[::-1]

def revcom(seq,nuc='DNA'):
	if nuc=='DNA':
		complement={'N':'N','A':'T','C':'G','G':'C','T':'A','a':'t','t':'a','c':'g','g':'c', 'U':'A', 'u':'a'}
	else:
		complement={'N':'N','A':'U','C':'G','G':'C','U':'A','a':'u','u':'a', 'c':'g','g':'c'}
	return ''.join(reversed([complement[base] for base in seq]))

def meanvar(value_list, vartype='sample'):
	import math
	n=len(value_list)
	total=math.fsum(value_list)
	square_sum=math.fsum([value**2 for value in value_list])
	mean=total/n
	if n>1 and vartype=='sample':
		variance=(square_sum-(2*mean*total)+(mean**2)*n)/(n-1)
	else:
		variance=(square_sum/n)-(mean**2)

	return [mean, variance]

def grabcolumn(path, column):
	outlist=[]
	for line in open(path):
		split_line=line.strip().split()
		outlist.append(split_line[column])
	return outlist

def gather_columns(path, columns=None, split_thing='\t'):
	newlist=[]
	for line in open(path):
		line_split=line.strip().split(split_thing)
		if columns==None:
			newline=line_split
		else:
			newline=[line_split[column] for column in columns]
		newlist.append(newline)
	return newlist

def read_fasta(fasta_file):
	seq, name_list, seq_list, seq_dict='', [], [], {}
	for line in open(fasta_file):
		line=line.strip()
		if '>' in line:
			name_list.append(line[1:])
			if len(seq)>0:
				seq_list.append(seq)
				seq=''
		else:
			seq=seq+line
	seq_list.append(seq)
#	for seq_number, name in enumerate(name_list):
#		seq_dict[name]=seq_list[seq_number]
	return [[name, seq_list[name_number]] for name_number, name in enumerate(name_list)]

def read_fasta_generator(fasta_name):
	'''
	Works like read_fasta but works as a generator instead of as a list. This
	allows large fasta files to be read in line by line without creating a large
	memory usage. Stolen from biopython via a stack overflow forum.
	'''
	import gzip
	if fasta_name.endswith('.gz'):
		fp=gzip.open(fasta_name, mode='rt')
	else:
		fp=open(fasta_name)
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line[1:], []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))


def print_fasta(fasta_list, outfile, mode='w', line_chars=60):
	"""
	this program prints fasta paired lists to fasta format
	"""
	output=open(outfile, mode)
	for sequence in fasta_list:
		output.write(">"+sequence[0]+"\n")
		for char_number, char in enumerate(sequence[1]):
			output.write(char)
			if char_number%line_chars==line_chars-1 or char_number==(len(sequence[1]))-1:
				output.write("\n")

def read_fastq(fastq_file):
	final_list=[]
	for line_number, line in enumerate(open(fastq_file)):
		line=line.strip()
		if line_number%4==0:
			name=line[1:]
		elif line_number%4==1:
			seq=line
		elif line_number%4==2:
			final_list.append([name, seq])
	return final_list

def read_fastq_generator(fastq_file):
	import gzip
	if fastq_file.endswith('.gz'):
		file_handle=gzip.open(fastq_file, mode='rt')
	else:
		file_handle=open(fastq_file)
	for line_number, line in enumerate(file_handle):
		line=line.strip()
		if line_number%4==0:
			name=line[1:]
		elif line_number%4==1:
			seq=line
		elif line_number%4==2:
			yield [name, seq]

def print_phylip(species_names, sequence_list, phylipfile, length=500):
	#prints paired name/sequence pairs in sequence_list to an interleaved output file
	#(phylipfile) and names each with names from species_names. Unlike the paired
	#fasta_list, the paired sequence_list expects a list of all names as item 0
	#and a list of all sequences as item 1
	output=open(phylipfile, 'w')
	total_length=len(sequence_list[1][0])
	output.write(str(len(species_names))+' '+str(total_length)+'\n')
	for species_number, species in enumerate(species_names):
		output.write(species_names[species_number]+' '+sequence_list[1][species_number][:length]+'\n')
	output.write('\n')
	for letter_number, letter in enumerate(sequence_list[1][0][length:]):
		if letter_number%length==0:
			for sequence in sequence_list[1]:
				sequence=sequence[length:]
				output.write(sequence[letter_number:letter_number+length]+'\n')
			output.write('\n')

def print_phylip2(species_names, sequence_list, phylipfile, length=500):
	#prints paired name/sequence pairs in sequence_list to an interleaved output
	#file (phylipfile) and names each with names from species_names. This
	#version accepts a standard fasta list as input
	sequence_list=[[name for name, seq in sequence_list], [seq for name, seq in sequence_list]]
	output=open(phylipfile, 'w')
	total_length=len(sequence_list[1][0])
	output.write(str(len(species_names))+' '+str(total_length)+'\n')
	for species_number, species in enumerate(species_names):
		output.write(species_names[species_number]+' '+sequence_list[1][species_number][:length]+'\n')
	output.write('\n')
	for letter_number, letter in enumerate(sequence_list[1][0][length:]):
		if letter_number%length==0:
			for sequence in sequence_list[1]:
				sequence=sequence[length:]
				output.write(sequence[letter_number:letter_number+length]+'\n')
			output.write('\n')


def print_alignment(fasta_thing, block_size=60):
	#reads an aligned fasta file and prints it 
	#interleaved in blocks of size block_size
	if type(fasta_thing) is str:
		fasta_thing=read_fasta(fasta_thing)
	else:
		fasta_list=fasta_thing
	for sequence in fasta_list:
		print(sequence[0])
	for letter_number, letter in enumerate(fasta_list[0][1]):
		if letter_number%block_size==0:
			for sequence in fasta_list:
				print(sequence[1][letter_number:letter_number+block_size])
			print("\n",)

def split_to_nodes (nodes, unsplit_path, output_folder, master_path):
	import subprocess
	import math
	output_number=1
	new_master_command=open(master_path, 'w')
	unsplit_file=open(unsplit_path)
	job_list=unsplit_file.readlines()
	total_jobs=float(len(job_list))
	divisor=math.ceil(total_jobs/nodes)
	for number, line in enumerate(job_list):
		if number%divisor==0:
			outfile=open("{0}{1}.sh".format(output_folder, output_number), 'w')
			new_master_command.write("{0}{1}.sh\n".format(output_folder, output_number))
			output_number+=1
		if number%divisor==divisor-1 or number==total_jobs-1:
			subprocess.call(["chmod", "+x", "{0}{1}.sh".format(output_folder, output_number-1)])
		outfile.write(line)
	return output_number-1
	
def theta_w(individuals, segsites, sites_w_info):
	#this program takes the number of individuals in the sample,
	#the number of segregating sites in the sample,
	#and the number of total sites considered to calculate theta_w
	a_1=0
	S=float(segsites)/sites_w_info #bases theta on total sites considered rather than total sites in locus
	for value in range(1,individuals):
		a_1=a_1+1.0/value
	theta=S/a_1
	return theta

def replace_all(input_list, conversion_list):
	#this program replaces all first elements from
	#conversion_list found in input_list with all
	#second elements from conversion_list
	for converting_line in conversion_list:
		for in_number, in_line in enumerate(input_list):
			if converting_line[0] in in_line:
				input_list[in_number]=in_line.replace(converting_line[0], converting_line[1])
	return input_list

def overlapping_search(search, string):
	#this definition searches 'string' using search term 'search'
	#and returns the coordinates of all overlapping matches, allowing for gaps
	import re
	expression, starts, ends="", [], []
	gap="-*"
	for letter in search[:-1]:
		expression+=letter+gap
	expression+=search[-1]
	for m in re.finditer("(?="+expression+")", string):
		start=m.start()
		starts.append(start)
		ends.append((re.match(expression, string[start:])).end()+start)
	return [str(start)+':'+str(ends[start_number]) for start_number, start in enumerate(starts)]

def old_choose(n, r):
	#counts ways of choosing r things from a set of n
	import math
	from operator import mul
	if n<r or n<0 or r<0:
		return 0
	if n-r>r:
		factorial_divisor=r
	else:
		factorial_divisor=n-r
	num_list=range(1, n+1)[:-(factorial_divisor+1):-1]
	denominator=math.factorial(factorial_divisor)
	if r>0 and r<n:
		numerator=reduce(mul, num_list)
	else:
		numerator=denominator
	return numerator/denominator

def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke, from StackOverflow.
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0
def binomial_prob(prob, number_observed, number_draws):
	'''
	calculates the binomial probability of seeing number_observed events in 
	number_draws total events when the probability of a single event is prob
	'''
	from decimal import Decimal
	ways=Decimal(choose(number_draws, number_observed))
	return Decimal(str(prob**number_observed*(1-prob)**(number_draws-number_observed)))*Decimal(ways)

def gcf(smaller, bigger):
	while smaller>0:
		smaller, bigger=bigger%smaller, smaller
	return bigger

def binomial_prob2(prob_num, prob_den, number_observed, number_draws, fraction=False):
	'''
	a modified (hopefully more precise) means of calculating binomial probabilities
	'''
	from decimal import Decimal
	
	ways=choose(number_draws, number_observed)
	not_prob_num=prob_den-prob_num
	not_observed=number_draws-number_observed
	prob_num=prob_num**number_observed*not_prob_num**not_observed*ways
	prob_den=prob_den**number_observed*prob_den**not_observed
	if fraction:
		return [prob_num, prob_den]
	else:
		return Decimal(prob_num)/Decimal(prob_den)

def cum_binom_prob(prob_num, prob_den, number_observed, number_draws):
	'''
	returns the probability of observing number_observed or more successes in
	number_draws trials given success probability prob_num/prob_den
	'''
	from decimal import Decimal
	if number_observed==0:
		return Decimal(1)
	else:
		for number in range(number_observed):
			if number==0:
				num, den=binomial_prob2(prob_num, prob_den, number, number_draws, 1)
			else:
				new_num, new_den=binomial_prob2(prob_num, prob_den, number, number_draws, 1)
				num, den=num*new_den+new_num*den, den*new_den
				factor=gcf(num, den)
				num, den=num/factor, den/factor
		return Decimal(den-num)/Decimal(den)
	
def check_overlap(start1, end1, start2, end2):
	#checks whether the start and end from the first two entries
	#overlaps with the start and end from the second two entries
	#returns "status" of overlap and start and end coordinates of overlapping bases
	if start1<=start2:
		bigger_small=start2
		status='smaller'
	else:
		bigger_small=start1
		status='bigger'
	if end1<=end2:
		smaller_big=end1
	else:
		smaller_big=end2
	if smaller_big>bigger_small:
		status='overlapping'
	return bigger_small, smaller_big, status

def serpentine(ranked_list, divisor):
	#takes a ranked list and distributes it to divisor
	#sub lists of approximately equal average ranks
	out_list, number=[], 0
	for element in range(divisor):
		out_list.append([])
	for rank_number, value in enumerate(ranked_list):
		if rank_number%divisor==0:
			step=0
		elif (rank_number/divisor)%2==0:
			step=1
		elif (rank_number/divisor)%2==1:
			step=-1
		number=number+step
		out_list[number].append(value)
	return out_list

def fasta_to_fastq(fasta_file_name, fastq_name=None, quality_score='z'):
	'''
	converts fasta files into fastq, All fastq entries are given a quality
	score of 'z'.
	'''
	if not fastq_name:
		fastq_name=fasta_file_name+'.fastq'
	fastq_file=open(fastq_name, 'w')
	fasta_list=read_fasta_generator(fasta_file_name)
	for name, seq in fasta_list:
		fastq_file.write('@'+name+'\n')
		fastq_file.write(seq+'\n')
		fastq_file.write('+\n')
		fastq_file.write('z'*len(seq)+'\n')
	fastq_file.close()
def fastq_to_fasta(fastq_file, fasta_file=None):
	#converts fastq files to fasta
	if not fasta_file:
		fasta_file=open(fastq_file+'.fa', 'w')
	else:
		fasta_file=open(fasta_file, 'w')
	for line_number, line in enumerate(open(fastq_file)):
		if line_number%4==0 and '@' in line:
			fasta_file.write('>'+line[1:])
		elif line_number%4==0 and '@' not in line:
			print('error')
			exit()
		elif line_number%4==1:
			fasta_file.write(line)
def locate_differences(file1, file2):
	'''
	finds first line where two files differ
	'''
	file1=open(file1)
	file2=open(file2)
	file1_line=file1.readline()
	file2_line=file2.readline()
	line_number=0
	while file1_line==file2_line:
		file1_line=file1.readline()
		file2_line=file2.readline()
		line_number+=1
	print(file1_line)
	print(file2_line)
	print(line_number)
def compute_pi(sequences):
	#computes the popgen statistic pi for a list of DNA sequences
	good_nucs='ACGTacgt'
	diff_count, bad_count=0, 0
	for sequence_number, sequence_one in enumerate(sequences):
		for sequence_two in sequences[sequence_number+1:]:
			for bp_number, bp_one in enumerate(sequence_one):
				bp_two=sequence_two[bp_number]
				if bp_two!=bp_one and bp_one in good_nucs and bp_two in good_nucs:
					diff_count+=1
				if bp_one not in good_nucs or bp_two not in good_nucs:
					bad_count+=1
	total_comparisons=float(choose(len(sequences), 2))*len(sequences[0])
	if bad_count<total_comparisons:
		return diff_count/(total_comparisons-bad_count), total_comparisons-bad_count
	else:
		return -0.02, total_comparisons-bad_count
def compute_pi2(sequences):
	#computes the popgen statistic pi for a list of DNA sequences
	good_nucs='ACGTacgt'
	diff_count, bad_count=0, 0
	for sequence_number, sequence_one in enumerate(sequences):
		for sequence_two in sequences[sequence_number+1:]:
			for bp_number, bp_one in enumerate(sequence_one):
				bp_two=sequence_two[bp_number]
				if bp_two!=bp_one and bp_one in good_nucs and bp_two in good_nucs:
					diff_count+=1
				if bp_one not in good_nucs or bp_two not in good_nucs:
					bad_count+=1
	total_comparisons=float(choose(len(sequences), 2))*len(sequences[0])
	return diff_count, total_comparisons, bad_count
def translate(genetic_sequence):
	#translates some DNA or RNA sequence into peptides
	translated=''
	translations={'UUU':'F', 'CUU':'L', 'AUU':'I', 'GUU':'V', 'UUC':'F', 'CUC':'L',
	'AUC':'I', 'GUC':'V', 'UUA':'L', 'CUA':'L', 'AUA':'I', 'GUA':'V', 'UUG':'L',
	'CUG':'L', 'AUG':'M', 'GUG':'V', 'UCU':'S', 'CCU':'P', 'ACU':'T', 'GCU':'A',
	'UCC':'S', 'CCC':'P', 'ACC':'T', 'GCC':'A', 'UCA':'S', 'CCA':'P', 'ACA':'T',
	'UCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A', 'UAU':'Y', 'CAU':'H', 'AAU':'N',
	'GAU':'D', 'UAC':'Y', 'CAC':'H', 'AAC':'N', 'GAC':'D', 'UAA':'*', 'CAA':'Q',
	'AAA':'K', 'GAA':'E', 'UAG':'*', 'CAG':'Q', 'AAG':'K', 'GAG':'E', 'UGU':'C',
	'CGU':'R', 'AGU':'S', 'GGU':'G', 'UGC':'C', 'CGC':'R', 'AGC':'S', 'GGC':'G',
	'UGA':'*', 'CGA':'R', 'AGA':'R', 'GGA':'G', 'UGG':'W', 'CGG':'R', 'AGG':'R',
	'GGG':'G', 'GCA':'A'}
	if 'T' in genetic_sequence:
		genetic_sequence=genetic_sequence.replace('T', 'U')
	for start in range(len(genetic_sequence)//3):
		codon=genetic_sequence[start*3:start*3+3]
		translated+=translations[codon]
	return translated
def make_ranks(x, y):
	'''
	takes paired x and y values and returns ranks, with correct tie handling
	'''
	import scipy.stats
	x=scipy.stats.stats.rankdata(x)
	y=scipy.stats.stats.rankdata(y)
	return x, y

def round_scientific(input_number, precision=3):
	'''
	rounds extremely small or large numbers in scientific notation intelligently
	'''
	number_str=str(input_number)
	number_list=number_str.split('e')
	number=number_list[0]
	if len(number_list)>1:
		exp=number_list[1]
	number=str(round(float(number), precision))
	if len(number_list)>1:
		return float('e'.join([number, exp]))
	else:
		return float(number)

def scatter_it(first, second, first_label, second_label, lims=False, faded=True):
	'''
	takes a list of x values ('first') and y values ('second'), along with x
	labels and y labels, and plots them as a scatter plot. Unless specified,
	limits are auto generated, and data points are faded. Figure will be saved
	as x_label_vs_y_label.png. A linear regression line is also plotted and
	listed as the title of the graph.
	'''
	import scipy.stats
	import matplotlib as mpl
	mpl.use('Agg') #comment this in if using a nongraphical interface, comment out to enable graphical interface
	import matplotlib.pyplot as plt
	import numpy
	m, b, r, linear_p_value, std_err=scipy.stats.linregress(first, second)
#	r, rank_p_value=scipy.stats.spearmanr(first, second)
	y=[m*x+b for x in first]
	if lims!=False:
		plt.xlim(lims[0])
		plt.ylim(lims[1])
	if faded:
		plt.scatter(first, second, color='red', linewidth=0, alpha=0.08)
	else:
		plt.scatter(first, second, color='red', linewidth=0)		
	plt.plot(first, y, color="black")
#	plt.tick_params(labelsize=17)
	plt.xlabel(first_label)
	plt.ylabel(second_label)
#	ax = plt.gca()
#	ax.tick_params(direction='out')
#	plt.xlabel(first_label, size=15)
#	plt.ylabel(second_label, size=15)
#	linear_p_value=numpy.round(linear_p_value, )
	plt.title('m='+str(round_scientific(m))+' b='+str(round_scientific(b))+' r2='+str(round_scientific(r**2))+' p='+str(round_scientific(linear_p_value)), size=15)
#	plt.title('m='+str(round_scientific(m))+' b='+str(round_scientific(b))+' r2='+str(round_scientific(r**2))+' p='+str(round_scientific(rank_p_value)), size=15)
#	plt.grid(True, linestyle='-')
	plt.savefig(first_label+'_vs_'+second_label+'.pdf', bbox_inches='tight')
	plt.clf()

def MannU(x_values, y_values):
	'''
	takes some x_values and y_values, ranks them together, and outputs which
	has a lower sum of ranks, and a P-value
	'''
	import scipy.stats
	merged_values=x_values+y_values	
	ranks, junk=make_ranks(merged_values, y_values)
	x_ranks=ranks[:len(x_values)]
	rank_sum=sum(x_ranks)
	rank_count=len(x_ranks)
	x_U=rank_sum-((rank_count*(rank_count+1))/2)
	y_U=len(x_values)*len(y_values)-x_U
	smaller_U, p=scipy.stats.mannwhitneyu(x_values, y_values)
	if x_U<y_U:
		statement='x_smaller'
	if y_U<x_U:
		statement='y_smaller'
	if y_U==x_U:
		statement='equal'
	return statement, smaller_U, p

def hypergeometric_prob(pop_success, pop_size, number_observed, number_draws, fraction=False):
	'''
	returns the probability of observing exactly number_observed successes in
	number_draws trials given a total of pop_success successes and a total
	population of pop_size
	'''
	from decimal import Decimal
	success_ways=choose(pop_success, number_observed)
	fail_ways=choose(pop_size-pop_success, number_draws-number_observed)
	numerator=success_ways*fail_ways
	total_ways=choose(pop_size, number_draws)
	if fraction==False:
		return Decimal(numerator)/Decimal(total_ways)
	else:
		return numerator, total_ways

def cum_hypergeometric_prob(pop_success, pop_size, number_observed, number_draws, extra=False):
	from decimal import Decimal
	'''
	returns the probability of observing number_observed successes or more in
	number_draws. Results should be	identical to a one-sided Fisher's exact
	test. (or a two sided if extra=True)
	'''
	num, den=0,1
	for number in range(number_observed):
		new_num, new_den=hypergeometric_prob(pop_success, pop_size, number, number_draws, 1)
		num, den=num*new_den+new_num*den, den*new_den
		factor=gcf(num, den)
		num, den=num/factor, den/factor
	regular=Decimal(den-num)/Decimal(den)
	if extra:
		enrich_value=Decimal(den-num)/Decimal(den)
		exact_num, exact_den=hypergeometric_prob(pop_success, pop_size, number_observed, number_draws, 1)
		deplete_value=Decimal(num)/Decimal(den)+Decimal(exact_num)/Decimal(exact_den)
		if enrich_value<deplete_value:
			return enrich_value, 'enriched'
		else:
			return deplete_value, 'depleted'
	else:
		return regular

def make_float(input_list):
	'''
	converts all numbers in a list of strings (like those read in from input
	files) into floats. Like map(float, input_list) but no error if input isn't
	float
	'''
	for column_number, column in enumerate(input_list):
		try:
			column=float(column)
			input_list[column_number]=column
		except ValueError:
			pass
	return input_list

def get_median(input_list):
	middle_value=len(input_list)/2
	if len(input_list)%2==0:
		median=(input_list[middle_value]+input_list[middle_value-1])/2.0
	else:
		median=input_list[middle_value]
	return median

def parse_csv_line(input_line):
	parsed_list, protected=[],False
	input_line=input_line.strip()
	while len(input_line)>0:
		for character_number, character in enumerate(input_line):
			if character=='"':
				if protected==False:
					protected=True
				else:
					protected=False
					parsed_list.append(input_line[1:character_number])
					input_line=input_line[character_number+2:]
					break
			if character==',' and protected==False:
				parsed_list.append(input_line[:character_number])
				input_line=input_line[character_number+1:]
				break
		if ',' not in input_line and '"' not in input_line:
			parsed_list.append(input_line)
			input_line=''
	return parsed_list

def product(input_list):
	product=1
	for item in input_list:
		product*=item
	return product

def iterate(nested_input_list):
	'''
	takes an arbitrarily deeply nested input list and iterates through all of
	the items in it
	'''
	sizes, final_list=[],[[],[]]
	for small_list in nested_input_list:
		sizes.append(len(small_list))
	total_size=product(sizes)
	for number in xrange(total_size):
		index_list,value_list=[],[]
		for size_number, current_size in enumerate(sizes):
			correction=product(sizes[size_number+1:])
			index=number/correction%current_size
			index_list.append(index)
			value_list.append(nested_input_list[size_number][index])
		final_list[0].append(index_list)
		final_list[1].append(value_list)
	return final_list

def smart_shuffle(input_list):
	'''
	randomly but efficiently shuffles an input list so that output list items
	occupy positions other than their original ones
	'''
	import random
	original_set=set(range(len(input_list)))
	new_indices=[]
	for index in range(len(input_list)):
		new_choice=random.choice(list(original_set-set([index])))
		original_set=original_set-set([new_choice])
		new_indices+=[new_choice]
	new_list=[0 for item in new_indices]
	for item_number, index in enumerate(new_indices):
		new_list[index]=input_list[item_number]
	return new_list

def make_count_dict(iterable, count_dict={}):
	'''
	nothing revolutionary here, returns the number of times each item in
	iterable occurs. I got tired of retyping this type of function every time
	'''
	for thing in iterable:
		if thing not in count_dict:
			count_dict[thing]=0
		count_dict[thing]+=1
	return count_dict

def make_list_dict(key, value, list_dict=None):
	'''
	appends the value to a list in a dictionary of lists. Chooses which list to
	append to using the key
	'''
	if list_dict==None:
		list_dict={}
	if key not in list_dict:
		list_dict[key]=[]
	list_dict[key].append(value)
	return list_dict

def use_ctl_file(control_file=None):
	'''
	a function that reads variable names from a control file and returns them
	as a dictionary with strings as values. Assumes a variable called
	'output_folder' is always present for the control file to write to, and
	that the control file's name will always be identical to the script name
	plus '.ctl'
	'''
	import subprocess
	import sys
	if not control_file:
		if len(sys.argv)>1:
			control_file=sys.argv[1]
		else:
			control_file=sys.argv[0][:-3]+'.ctl'
	v_list=[line.strip().split('#')[0].split('=') for line in open(control_file)]
	v=dict([[term.strip() for term in item] for item in v_list])
	output_folder=v['output_folder']
	subprocess.call(['mkdir', '-p', output_folder])
	subprocess.call(['cp', '-u', control_file, output_folder])
	if 'extra_parameters' in v:
		subprocess.call(['cp', v['extra_parameters'], output_folder])		
	return v

def make_ctl_file(v_dict, file_name):
	'''
	The reverse of use_ctl_file. Takes a dictionary of variable names and
	stored strings and creates a ctl file called file_name that uses these
	variable names
	'''
	ctl_file=open(file_name, 'w')
	for v in v_dict:
		ctl_file.write(f'{v}={v_dict[v]}\n')

def yield_in_base(string, valuestring, separator='|'):
	'''
	designed to replace count_in_base by creating a generator function instead
	of a list. Also eliminates the superfluous 'base' argument and the
	requirement of a user supplied separator character as the last character of
	the valuestring. Note: do not use '|' as one of the characters in the
	valuestring unless you redefine the separator as something else.
	'''
	base=len(valuestring)
	valuestring=valuestring+separator
	full_length=len(string)
	for count in range(base**len(string)):
		yield string
		if string==valuestring[-2]*full_length:
			break
		position=len(string)-1
		value=valuestring.index(string[position])
		string=string[:position]+valuestring[value+1]+string[position+1:]
		while value>=(base-1) and count<(base**len(string)):
			string=string[:position]+valuestring[0]+string[position+1:]
			position-=1
			value=valuestring.index(string[position])
			string=string[:position]+valuestring[value+1]+string[position+1:]

def get(in_dict, key, null_value='none_provided'):
	'''
	works the same as 'dict.get(key, [null])' but more efficiently
	'''
	if key not in in_dict:
		return null_value
	else:
		return in_dict[key]

def slurp_seqs(file_name_list):
	'''
	the goal of this program is to take a list of sequence files (could be fasta
	files or fastq files, gzipped or not gzipped) and to yield a single nested
	list suitable for sending to fastq or fasta format or for iterating through.
	'''
	for file_name in file_name_list:
		altered_name=file_name
		print(altered_name)
		if file_name.endswith('.gz'):
			altered_name=file_name[:-3]
		if altered_name.endswith('.fasta') or altered_name.endswith('.fa') or altered_name.endswith('.fst'):
			new_list=read_fasta_generator(file_name)
		elif altered_name.endswith('.fq') or altered_name.endswith('.fastq'):
			new_list=read_fastq_generator(file_name)
		for paired_seq in new_list:
			yield paired_seq

def pie_chart(title, labels, sizes):
	'''
	makes a pie chart labeled with title out of data with labels and sizes given
	as lists. Output chart is a pdf. Template is from here:
	https://matplotlib.org/stable/gallery/pie_and_polar_charts/pie_features.html
	'''
	import matplotlib as mpl
	mpl.use('Agg') #comment this in if using a nongraphical interface, comment out to enable graphical interface
	import matplotlib.pyplot as plt
	# Pie chart, where the slices will be ordered and plotted counter-clockwise:
	colors=['blue', 'darkorange', 'green', 'yellow', 'red', 'purple', 'brown', 'lime', 'pink', 'gray', 'fuchsia', 'olive', 'cyan']
	summed=sum(sizes)
	fractions=[size/summed for size in sizes]
	percentages=[round(size*100, 1) for size in fractions]
	labels=[label+'_'+str(percentages[label_number])+'%' for label_number, label in enumerate(labels)]
	if len(sizes)%len(colors)==1:
		colors=colors[:-2]
	#explode = (0, 0, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
	fig1, ax1 = plt.subplots()
	#ax1.pie(sizes, explode=explode, labels=labels, counterclock=False, autopct='%1.1f%%',
	#        shadow=True, startangle=90)
	ax1.pie(sizes, labels=labels, colors=colors, counterclock=False, autopct='%1.1f%%',
		    startangle=90, pctdistance=0.9, labeldistance=1.1, textprops={'fontsize': 6})
	ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
	plt.title(title)
	plt.legend(loc='upper right', bbox_to_anchor=(-0.3, 1.))
	plt.savefig(title+'.png', bbox_inches='tight')
#	plt.savefig(title+'.pdf')
	plt.clf()


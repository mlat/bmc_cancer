def amplification_deletion(smfile):
	#gives time to see how long function takes to run and initalizes all arrays
	import time;
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	cosmic_list = cosmicgenes()
	sm_data = open(smfile, 'r')
	smdata = sm_data.readlines();
	chromdata = open(out_er_pr_her2_expr_type.txt, r)
	bpdata = chromdata.readlines();
	data = []
	chrom_bp = []
	chrom_d = {}
	for line in smdata: 
		cols = line.split(',')
		gene = cols[0]
		chrom = cols[4]
		start = cols[5]
		end = cols[6]
		data.append([gene, chrom, start, end])
	
	for gene in cosmic_list: 
		key = gene
		for i in data: 
			if i[0] == key: 
				chrom_d[key].add(i[1], i[2], i[3])
			else: 
				chrom_d[key] = {i[1], i[2], i[3]}
	for line in bpdata:
		tabs = line.split('\t')
		genename = tabs[15:]
		chrom_bp.append([genename, tabs[0], tabs[1], tabs[2]])
	
	#for key in chrom_d: 
		#if key == chrom_bp[0]: 
			
		
		
def gene_network(sorteddatafile): 
	#gives time to see how long function takes to run and initalizes all arrays
	import time;
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	import numpy;
	import math
	data = age_genefrequency_table(sorteddatafile)
	over70_data = data[3:]
	ages = [50,60,70,80]
	under70 = len(data[0:3])
	over70 = len(over70_data)
	gene_length = len(data[0][3])
	total_and_genes = []
	total_and_genes_70 = []
	toptengenes = []
	toptengenes_70 = []
	topgene_compare = numpy.zeros([10,10])
	topgene_compare_70 = numpy.zeros([10,10])
	cosmicset = cosmicgenes()
	somaticgene_d = get_mutations(sorteddatafile)
	j = 1
	k = 1
	x = 1
	w = 16
	#opens current excel file and writes headers for matrix  
	from xlrd import open_workbook
	from xlutils.copy import copy
	book = open_workbook('genecomparison2.xls', formatting_info=True)
	wbook = copy(book)
	sheet11 = wbook.get_sheet(10)
	sheet12 = wbook.get_sheet(11)
	rowgenename = sheet11.row(0)
	over60genename = sheet11.row(15)
	for g in range(gene_length): 
		total = 0 
		total_70 = 0 
		for l in range(under70):
			total = total + data[l][2][g]
		total_and_genes.append([total,cosmicset[g]])
		for m in range(over70): 
			total_70 = total_70 + over70_data[m][2][g]
		total_and_genes_70.append([total_70,cosmicset[g]])	
			
	total_and_genes.sort()
	total_and_genes.reverse()
	total_and_genes_70.sort()
	total_and_genes_70.reverse()
	print 'total array: ', total_and_genes, '\n', '70 total array: ', total_and_genes_70
	toptotal_w_genes = total_and_genes[:10]
	toptotal_w_genes_70 = total_and_genes_70[:10]
	print 'top totals: ', toptotal_w_genes, '\n', '70 top totals: ', toptotal_w_genes_70
	for genename in toptotal_w_genes:
		toptengenes.append(genename[1])
	for genes in toptotal_w_genes_70:
		toptengenes_70.append(genes[1])
	print 'top ten genes: ', toptengenes, '\n', '70 top ten genes: ', toptengenes_70
	for item in toptengenes:
		rowgenename.write(j,item)
		sheet11.write(k,0,item)
		j = j + 1
		k = k + 1
	for point in toptengenes_70: 
		over60genename.write(x, point)
		sheet11.write(w,0, point)
		x = x + 1
		w = w + 1
	for key in somaticgene_d: 
		for gene in range(len(toptengenes)): 
			if toptengenes[gene] in somaticgene_d[key] and key[1] < 70: 
				for i in range(len(toptengenes)): 
					if toptengenes[i] in somaticgene_d[key]: 
						topgene_compare[gene][i] = topgene_compare[gene][i] + 1	
						
	for r in range(topgene_compare_70.shape[0]): 
		for c in range(topgene_compare.shape[1]): 
			sheet11.write(r+1,c+1,topgene_compare[r][c])
	
					
	print '10 by 10 matrix: ', topgene_compare
	print 'top ten genes: ' , toptengenes
	print 'toptotals:', toptotal_w_genes				
	
	for sample in somaticgene_d: 
		for genetics in range(len(toptengenes_70)): 
			if toptengenes_70[genetics] in somaticgene_d[sample] and sample[1] >= 70: 
				for f in range(len(toptengenes_70)): 
					if toptengenes_70[f] in somaticgene_d[sample]: 
						topgene_compare_70[genetics][f] = topgene_compare_70[genetics][f] + 1
						
	print '10 by 10 matrix over 70  ', topgene_compare_70
	print 'top ten genes over 70: ' , toptengenes_70
	print 'toptotals_70:', toptotal_w_genes_70				
	
	for row in range(topgene_compare_70.shape[0]): 
		for column in range(topgene_compare_70.shape[1]): 
			sheet11.write(row+16,column+1,topgene_compare_70[row][column])
	

	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	wbook.save('genecomparison2.xls')




#The total_missensemutation_sample function counts the total number of missesnse mutations per individual 
def total_missensemutation_sample(sorteddatafile):
	import time;
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	sorteddata = open(sorteddatafile, 'r')
	data = []
	from xlrd import open_workbook
	from xlutils.copy import copy
	book = open_workbook('genecomparison2.xls', formatting_info=True)
	wbook = copy(book)
	sheet10= wbook.get_sheet(9)
	#polyphen2_out = open('polyphen2-output.txt','r')
	#polyphen2 = polyphen2_out.readlines()
	for line in sorteddata:
		cols = line.split(',')
		variant = cols[3]
		mutation = variant.rstrip('\n')
		data.append([cols[0], cols[1], cols[2], mutation])
	sample_dictionary= {} #initialized dictionary 
	#runs through each sample number in the data list
	no_missense = 0 
	a = 0
	for i in data:  
		if i[3] == "Missense_Mutation" and i[1] != '[Not Available]': #only takes data with a mutated gene
			no_missense = no_missense + 1
			key = (i[0], int(i[1]))
			if (key in sample_dictionary): #determines if sample number already exists and adds gene if it does
				sample_dictionary[key] = sample_dictionary[key] + 1
			else: # creates a set(allows for no repeats) of genes associted with sample number if sample number does not exist 
				sample_dictionary[key] = 1 
	sheet10.write(0,2,'Missense Mutations')
	sheet10.write(0,1,'Age Group')
	sheet10.write(0,0, 'Sample No.')
	rowgenename = sheet10.row(0)
	for i in sample_dictionary:
		a = a + 1
		sheet10.write(a,0, i[0])
		sheet10.write(a, 1, i[1]) 
		sheet10.write(a, 2, sample_dictionary[i]) 
	print 'Total Missense Mutations:', no_missense 
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	wbook.save('genecomparison2.xls')
	return sample_dictionary;	




# the missense vs age function determines how many missense mutations there are in each age group of 10 years.  
class missense_count:
	total = 1

def missense_vs_age(datasorted):
	datasorted_LUSC = open(datasorted, 'r')
	data = datasorted_LUSC.readlines();
	age_dictionary = {}
	ages = [50,60,70,80] 
	from xlrd import open_workbook
	from xlutils.copy import copy
	book = open_workbook('genecomparison2.xls', formatting_info=True)
	wbook = copy(book)
	sheet9= wbook.get_sheet(8)
	missense_total = 0
	a = 0 
	for line in data: 
		cols = line.split(',')
		if cols[1] != '[Not Available]':
			age = int(cols[1])
		mutation = cols[3].rstrip('\n')
		if mutation == "Missense_Mutation":
			missense_total = missense_total + 1
			decade = agegroup(age) # agegroup function groups any ages under 30 and any ages over 80 in that group 
			if decade in age_dictionary: 
				age_dictionary[decade].total = age_dictionary[decade].total + 1
			else: 
				age_dictionary[decade] = missense_count()
	print age_dictionary[decade]
	print age_dictionary
	print 'missense total: ', missense_total 
	sheet9.write(0,1,'Missense Mutations')
	sheet9.write(0,0,'Age Group')
	rowgenename = sheet9.row(0)
	for group in ages:
		a = a +1
		sheet9.write(a, 0, group) 
		sheet9.write(a, 1, age_dictionary[group].total) 
	wbook.save('genecomparison2.xls')		

# the age_vs_matches function plots the total number of matches for all 121 cosmic genes vs age for each sample in the BRCA file
def age_vs_matches(): 
	#opens excel file and allows for python to plot 
	import matplotlib.pyplot as graph
	from pylab import polyfit
	from pylab import poly1d
	from xlrd import open_workbook
	from xlutils.copy import copy
	book = open_workbook('genecomparison2.xls', formatting_info=True)
	wbook = copy(book)
	sheet8= wbook.get_sheet(7)
	genecompare = open('genecomparison.txt', 'r')
	ages_matches = []
	ages = []
	matches = []
	genecompare.readline()
	k = 1
	#determines the number of matches per sample, and creates a list of the ages and number of matches for that age 
	for line in genecompare: 
		match = 0 
		genes = line[2:]
		cols = line.split('\t')
		ages.append(int(cols[1]))
		for i in genes: 
			if i == 'y': 
				match = match + 1
		ages_matches.append([int(cols[1]), match])
		matches.append(match)
	print ages, len(ages)
	print matches, len(matches)
	print ages_matches, len(matches)
	#plots  a scatter plot of each age and total matches for that age in python
	graph.scatter(ages,matches)
	graph.show()
	#creates a trendline for a scatter plot 
	fit = polyfit(ages,matches,1)
	fit_fn = poly1d(fit) # fit_fn is now a function which takes in a age and returns an estimate for matches

	graph.plot(ages,matches, 'yo', ages, fit_fn(ages), '--k')
	graph.show()
	data = (ages,matches)
	graph.boxplot(ages,matches,vert = False)
	graph.ylabel('Number of Cancer Gene Mutations')
	graph.xlabel('Age')
	graph.title('Total number of Cancer Gene Mutations per age ')
	graph.plot(ages, fit_fn(ages), '--k')
	graph.show()
	sheet8.write(0,0, 'age')
	sheet8.write(0,1,'total number of cancer gene mutations')
	#writes a list of ages and matches so can also plot in excel and have tabulated data 
	for i in ages_matches:
		sheet8.write(k,0,i[0])
		sheet8.write(k,1,i[1])
		k = k + 1
	wbook.save('genecomparison2.xls')

# the age_vs_frequency function plots each the frequency of each gene as a function of age in multiple ways 
#*****THE DATAFILE these functions take is the sorted data file produced from data_sort that is a list of sample, age, gene, mutation
def age_vs_frequency(datafile): 
	#gives time to see how long function takes to run and initalizes all arrays
	import time;
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	import matplotlib.pyplot as graph
	data = age_genefrequency_table(datafile)
	ages = [50,60,70,80]
	agegroups = len(data)
	gene_length = len(data[0][3])
	legend_genes = []
	legend_genes_match = []
	all_important_genes = []
	all_important_genes_match = []
	cosmicset = cosmicgenes()
	print len(cosmicset), gene_length # makes sure the 2 sets are equal 
	#opens current excel file and writes headers for each table to plot graphs 
	from xlrd import open_workbook
	from xlutils.copy import copy
	book = open_workbook('genecomparison2.xls', formatting_info=True)
	wbook = copy(book)
	sheet3 = wbook.get_sheet(2)
	sheet4 = wbook.get_sheet(3)
	sheet5 = wbook.get_sheet(4) 
	sheet6 = wbook.get_sheet(5)
	sheet7 = wbook.get_sheet(6)
	sheet3.write(0,1,'No. of Samples')
	sheet3.write(0,0,'Age Group')
	sheet4.write(0,1,'No. of Samples')
	sheet4.write(0,0,'Age Group')
	sheet5.write(0,1,'No. of Samples')
	sheet5.write(0,0,'Age Group')
	sheet6.write(0,1,'No. of Samples')
	sheet6.write(0,0,'Age Group')
	sheet7.write(0,1,'No. of Samples')
	sheet7.write(0,0,'Age Group')
	rowgenename = sheet3.row(0)
	rowgenenamematch = sheet4.row(0)
	rowgenename0_5 = sheet5.row(0)
	rowgenename0_6 = sheet6.row(0)
	rowgenename0_7 = sheet7.row(0)
	a = 0
	n = 0
	j = 2
	p = 2
	w = 2 
	c = 2 
	x = 2
	#writes the age groups and the number of samples in each sheet 
	for age in ages:
		a = a +1
		sheet3.write(a, 0, age) 
		sheet3.write(a, 1, data[n][1]) 
		sheet4.write(a, 0, age)
		sheet4.write(a, 1, data[n][1])
		sheet5.write(a, 0, age)
		sheet5.write(a, 1, data[n][1])
		sheet6.write(a, 0, age)
		sheet6.write(a, 1, data[n][1])
		sheet7.write(a, 0, age)
		sheet7.write(a, 1, data[n][1])
		n = n + 1
		if a == 6:
			print "finished writing header of excel sheets 1-7"
	#creates a matrix with genes and their frequencies if the total frequency for that gene is greater than 7 and plots the results in matlab
	for g in range(gene_length): 
		total = 0 
		for l in range(agegroups):
			total = total + data[l][2][g]
		if total >= 15:
			important_gene = [i[2][g] for i in data]
			all_important_genes.append([cosmicset[g],important_gene])
			legend_genes.append(cosmicset[g])
			graph.plot(ages,important_gene)
	print 'graph is of important genes with frequency >= 15'
	graph.legend(legend_genes)
	graph.show()
	#creates a matrix with genes and the frequencies if total somatic matches is  greater than 7 for each cosmic gene and plots the results in matlab
	for g in range(gene_length):
		match_total = 0 
		for l in range(agegroups):
			match_total = match_total + data[l][3][g]
			print match_total 
		if match_total > 7:
			important_gene_match = [i[2][g] for i in data]
			all_important_genes_match.append([cosmicset[g],important_gene_match])
			legend_genes_match.append(cosmicset[g])
			graph.plot(ages,important_gene_match)
	print 'graph is of important genes with match > 10'
	graph.legend(legend_genes_match)
	graph.show()
	# prints the matches matrix and frequency matrix above in excel
	for item in all_important_genes:
		k = 1
		rowgenename.write(j,item[0])
		for f in range(len(ages)):
			sheet3.write(k,j,item[1][f])
			k = k + 1
		j = j + 1
	for item in all_important_genes_match:
		k = 1
		rowgenenamematch.write(p,item[0])
		for q in range(len(ages)): 
			sheet4.write(k,p,item[1][q])
			k = k + 1
		p = p + 1
	print 'finsished important_gene and important_gene_match'
	PIK3CA_TP53 = []
	PIK3CA_TP53_match = []
	legend_P = []
	# removes the PIK3CA and TP53 genes from the match matrix and frequency matrix and plots resulting genes in python 
	for item in all_important_genes: 
		#if item[0] =='PIK3CA':
			#PIK3CA_TP53.append(item)
		if item[0] =='TP53':
			PIK3CA_TP53.append(item)
	for item in PIK3CA_TP53:
		all_important_genes.remove(item)
		
	print 'now printing PIK3CA TP53'
	#print PIK3CA_TP53
	for item in all_important_genes: 
		graph.plot(ages,item[1])
	#legend_genes.remove('PIK3CA')
	legend_genes.remove('TP53')
	print 'graph is for all imporant genes without TP53 and PIK3CA'
	graph.legend(legend_genes)		
	graph.show()
	for item in all_important_genes: 
		k = 1
		rowgenename0_5.write(w,item[0])
		for f in range(len(ages)): 
			sheet5.write(k,w,item[1][f])
			k = k + 1 
		w = w + 1
	for item in all_important_genes_match: 
		#if item[0] =='PIK3CA':
			#PIK3CA_TP53_match.append(item)
		if item[0] =='TP53':
			PIK3CA_TP53_match.append(item)
	for item in PIK3CA_TP53_match: 
		all_important_genes_match.remove(item)
	print 'now printing PIK3CA TP53 match'
	#print PIK3CA_TP53_match
	for item in all_important_genes_match:	
		graph.plot(ages, item[1])
	legend_genes_match.remove('TP53')
	#legend_genes_match.remove('PIK3CA')
	print 'graph is for all important_genes_match with out TP53 and PIK3CA'
	graph.legend(legend_genes_match)
	graph.show()
	# writes a table of match matrix and frequency matrix without TP53 and PIK3CA genes in excel 
	for item in all_important_genes_match: 
		k = 1
		rowgenename0_6.write(c,item[0])
		for f in range(len(ages)): 
			sheet6.write(k,c,item[1][f])
			k = k + 1 
		c = c + 1
	#plots just the TP53 and PIK3CA gene in python and creates a table in excel of just TP53 and PIK3CA frequencies 
	for gene in PIK3CA_TP53: 
		legend_P.append(gene[0])
		graph.plot(ages,gene[1])
	graph.legend(legend_P)
	graph.show()
	for item in PIK3CA_TP53:
		k = 1
		rowgenename0_7.write(x,item[0])
		for f in range(len(ages)):
			sheet7.write(k,x,item[1][f])
			k = k + 1
		x = x + 1
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	wbook.save('genecomparison2.xls')

#the age_genefrequency_table function prints the age group in ranges of 10 years, the number of samples per age group and the frequency for each for each age group
def age_genefrequency_table(datafile):
	import time;
	localtime = time.asctime( time.localtime(time.time()) )
	#initializes age groups and arrays
	ages = [50,60,70,80]
	agedata = []
	#creates sheet in excel file to write table to 
	from xlrd import open_workbook
	from xlutils.copy import copy
	book = open_workbook('genecomparison2.xls', formatting_info=True)
	sheet = book.sheet_by_index(1)
	wbook = copy(book)
	sheet2 = wbook.get_sheet(1)
	sheet2.write(0,1,'No. of samples')
	sheet2.write(0,0,'age group')
	sheet2.write(9,2,'number of matches/gene')
	sheet2.write(10,1,'No. of samples')
	sheet2.write(10,0,'age group')
	a = 1
	j = 2 
	t = 1
	i = 1
	m = 11
	y = 11
	x = 11
	row0 = sheet2.row(0)
	row10 = sheet2.row(10)
	# prints the age groups in the first column 
	for age in ages:
		sheet2.write(a, 0, age)
		a = a +1 
		sheet2.write(y,0,age)
		y = y + 1
	cosmicset = cosmicgenes() # gets the cosmic cancer gene list 
	#writes the header for the genes in the excel file 
	for gene in cosmicset: 
		row0.write(j, str(gene))
		row10.write(j, str(gene))
		j = j + 1	
	# goes through the genedictionary by age group	
	for age in ages: 
		#s = str(age)
		gened = mutations_func_age(datafile)
		total = gened[gene][age].total # gets the total of samples for that age group 
		sheet2.write(t, 1, total)
		sheet2.write(x, 1, total)
		t = t + 1
		x = x + 1 
		g = 2
		frequency_array = []
		match_array = []
		#goes through each gene and determines the frequency for that age and prints the frequency for each gene 
		for gene in cosmicset: 
			match = gened[gene][age].match
			frequency = (float(match)/total)*100
			row_freqi = sheet2.row(i)
			row_freqi.write(g, frequency)
			frequency_array.append(frequency)
			match_array.append(match)
			row_match = sheet2.row(m)
			row_match.write(g, match)
			g = g + 1 
		i = i + 1
		m = m + 1
		#stores the frequency in an array with age and sample number 
		agedata.append([age,total, frequency_array, match_array])
	print 'creating age_genefrequency_table'
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	#wbook.save('frequencybyage.xls')
	wbook.save('genecomparison2.xls')
	return agedata
			
			
# the class match_frequency, stores the total amount of samples for each age group in each gene and the number of matached genes in that age group 		
class match_frequency:
	total = 1 #total starts at 1 so if the age group is added, that sample is counted towards the total 
	match = 0 #tallies the number of matched genes/ age group - it is initialized to 0 so that each age group starts with 0 matches 
	

# the mutations_func_age function creates a dicitonary of COSMIC genes- each gene in the dictionary has an age dictionary that stores how many samples there are for that decade and which samples have a match
#to the gene in the class match_frequency
def mutations_func_age(datafile):
	sampled = get_mutations(datafile) # gets the dictionary of sample numbers with somtic mutated genes associated
	gened = {} # initializes dictionary for each cancer gene 
	cosmicset = cosmicgenes()
	for gene in cosmicset: #runs through each gene in the cosmic set 
		age_statd = {} #creates dictionary that groups sample by ages 
		for key in sampled: #runs through each sample
			decade = agegroup(key[1]) #creates the age group by rounding down for the age in the sample 
			if decade in age_statd: #if the age group is already in the dictionary, it counts the sample for the total number of samples within in each age group 
				age_statd[decade].total = age_statd[decade].total + 1 
			else:
				age_statd[decade]= match_frequency() # if the age group is not already in the dictionary the age group is added
			if gene in sampled[key]: # if the cancer gene matches a gene in the sample with the somatic mutant gene list, counts the match 
				age_statd[decade].match = age_statd[decade].match + 1
		gened[gene] = age_statd	#stores the total samples and matches by age in each gene 
	print 'created gene dictionary with age groups, number of matches/group, and total samples/group for gene '
	return gened  
			
			
# The agegroup function takes an age and rounds it down to the decade		
def agegroup(x): 
	import math
	age = x
	decimal= age/10 
	group = math.floor(decimal) 
	age_group = group*10 	
	if age_group <= 50: 
		age_group = 50
	if age_group >= 90:
		age_group = 80
		
	return age_group
	

# the cosmic gene function imports the list of cosmic genes and stores them in an array 
def cosmicgenes(): 
	#opens file and intializes array 
	cosmicdata = open('cosmicgenes.txt', 'r')
	cosmicset = []
	cosmic = cosmicdata.readlines();
	#stores genes in comparable format in array 
	for line in cosmic:
		cols = line
		cosmicgene = cols.rstrip('\n')
		cosmicset.append(cosmicgene)
	return cosmicset

	
	#comparegenes function compares a list of genes associated with a sample number (somatic mutated genes) to a list of mutated genes known to cause cancer (cosmic genes) 
def comparegenes(datafile):
	#opens cosmic files and creates file to write 
	cosmicdata = open('cosmicgenes.txt', 'r')
	comparison = open('genecomparison.txt', 'w+')
	from xlrd import open_workbook
	from xlutils.copy import copy
	book = open_workbook('genecomparison2.xls', formatting_info=True)
	sheet = book.sheet_by_index(0)
	wbook = copy(book)
	sheet1 = wbook.get_sheet(0)
	sheet1.write(0,0,'sample')
	sheet1.write(0,1,'age')
	j = 2 
	k = 1 
	row0 = sheet1.row(0)
	comparison.write('sample\tage\t')
	#gets somatic mutation data with only mutated genes
	smdata = get_mutations(datafile)
	#initializes array of cosmic genes and fills array for comparison with somatic mutated genes
	cosmicset = []
	cosmic = cosmicdata.readlines();
	for line in cosmic:
		cols = line
		cosmicgene = cols.rstrip('\n')
		cosmicset.append(cosmicgene)
		row0.write(j, str(cosmicgene))
		j = j + 1
		comparison.write(str(cosmicgene)+ '\t')
	#compares somatic mutated genes to cosmic genes and writes in file if there is a match or not
	for key in smdata:
		rowk = sheet1.row(k)
		rowk.write(0, key[0])
		rowk.write(1, key[1])
		comparison.write('\n'+ key[0] + '\t' + str(key[1]) + '\t')
		genes = smdata[key]
		l = 2
		for i in cosmicset:
			if (i in genes):
				 rowk.write(l,"y")
				 comparison.write("y" + '\t')
			else: 
				 rowk.write(l,"n")
				 comparison.write("n" + '\t')
			l = l + 1
		k = k + 1
	wbook.save('genecomparison2.xls')
	comparison.close()
	return comparison

#get_mutations function creates a dictionary of sample numbers and the age associated with that sample number.  Each sample number has a list of only mutated genes associtaed with it 
def get_mutations(datafile):
	import time;
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	#data = sort_data(smfile, agefile) #retrieves sorted data
	sorteddata = open(datafile, 'r')
	data = []
	for line in sorteddata:
		cols = line.split(',')
		variant = cols[3]
		mutation = variant.rstrip('\n')
		data.append([cols[0], cols[1], cols[2], mutation])
	d= {} #initialized dictionary 
	#runs through each sample number in the data list
	for i in data:  
		if i[3] == "Missense_Mutation" and i[1] != '[Not Available]': #only takes data with a mutated gene
			key = (i[0], int(i[1]))
			if (key in d): #determines if sample number already exists and adds gene if it does
				d[key].add(i[2])
			else: # creates a set(allows for no repeats) of genes associted with sample number if sample number does not exist 
				d[key] = { i[2] }
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	return d;
	
#the sort_data function creates a list of data that only includes the sample number, age, gene, and variant of the gene
def sort_data(smfile, agefile, outputfilename):
	import time;
	#opens the original sample data and age data downloaded from the internet, and creates an array to store the new sorted data in.   
	alldata = open(smfile, 'r')
	allages = open(agefile, 'r')
	datasorted = open(outputfilename, 'w')
	smdata= alldata.readlines();
	agedata = allages.readlines();
	data = []
	print "Processing sort data"
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	for line in smdata:# parses through the data and chooses only the data that is desired.   
		cols = line.split('\t')
		sample_no = cols[15]
		sample = sample_no[0:12]
		variant = cols [8]
		gene = cols[0]
		for line in agedata: #parses through the age data and chooses only the age 
			colsa = line.split('\t')
			age = colsa[1]
			samplenoa = colsa[0]
			samplea = samplenoa[0:12]
			#print samplea, "age"
			if samplea == sample: #matches the correct age with the correct sample and stores the desired data into a list
				data.append([sample, age, gene, variant])
				datasorted.write(sample + ',' + age + ','+ gene + ','+ variant + '\n')
				#print data
	print ' done processing sort data' 
	localtime = time.asctime( time.localtime(time.time()) )
	print "Local current time :", localtime
	datasorted.close()
	

	
	

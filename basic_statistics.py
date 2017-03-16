import argparse
import math
import subprocess
import matplotlib.pyplot as plt
import matplotlib.collections as collections
import seaborn as sns
import pandas
import numpy as np
import statistics as stats
from matplotlib.backends.backend_pdf import PdfPages

def basic_stats(programLoc, inputFile, outputPrefix):
	filename = args.inputFile.rsplit('/', 1)[-1]
	if filename[-4:] == '.ped':
		subprocess.call([programLoc, '--file', inputFile[:-4], '--freq', '--missing', '--hardy', '--mendel', '--out', outputPrefix])
	elif filename[-4:] == '.bed':
		subprocess.call([programLoc, '--bfile', inputFile[:-4], '--freq', '--missing', '--hardy', '--mendel', '--out', outputPrefix])
	else:
		print("File is not a .ped or .bed, please convert")

def visualize_stats(outputPrefix, custom_file):
	
	basic_stats_analysis = open(outputPrefix+'-basic-analysis.txt', 'w')
	
	freq_file = pandas.read_table(outputPrefix + '.frq', delim_whitespace=True)
	all_names = list(freq_file['SNP'])
	basic_stats_analysis.write('The total number of all snps in the dataset for analysis is '+str(len(all_names))+'\n')

	# The total custom snps with same snp ID and more than one position
	custom_snp_offset = 0
	# read in custom_file if user supplies one
	if custom_file != None:
		custom_snp_dict = {}
		with open(custom_file) as custom_snp_file:
			for line in custom_snp_file:
				snp, chrm, pos = line.rstrip().split(',')
				if snp not in custom_snp_dict:
					custom_snp_dict[snp] = [chrm, pos]
				# necessary because occassionally a SNP can have multiple positions
				else:
					custom_snp_dict[snp] = custom_snp_dict[snp] + [pos]
					custom_snp_offset = custom_snp_offset + 1
		
		custom_names_only = [key for key in custom_snp_dict]
		basic_stats_analysis.write('The total number of custom SNPs is ' + str(len(custom_names_only)+custom_snp_offset)+'\n')
		basic_stats_analysis.write('There are ' + str(custom_snp_offset) + 'SNP IDs in the custom set that have multiple positions under the same ID. They are:' +'\n')
		basic_stats_analysis.write('\n'.join([snp for snp in custom_snp_dict if len(custom_snp_dict[snp]) > 2])+'\n')
		basic_stats_analysis.write('The total number of custom SNP IDs not found in the entire data set is ' + str(len(set(custom_names_only) - set(all_names)))+'\n')
		basic_stats_analysis.write('The IDs of these missing SNPs are as follows:'+'\n'+'\n'.join(list(set(custom_names_only) - set(all_names))))

	def maf_analysis():
		# distribution of SNP counts by MAF
		freq_file = pandas.read_table(outputPrefix + '.frq', delim_whitespace=True)
		plt.figure()
		maf_dist = freq_file['MAF'].plot.hist(bins=10)
		plt.xlabel("MAF scores")
		plt.ylabel("Total Numer of SNPs")
		plt.title("Distribution of MAF across all SNPs", fontsize=15)
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# graphing of MAF information by number SNPs retained at each threshold; answers the question "How many SNPs have a given MAF score within each MAF interval"
		MAF_iterables = [0.0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.20]
		total_snps = []
		total_custom_captured = []
		for MAFscore in MAF_iterables:
			total_snps.append(freq_file[freq_file['MAF'] >= float(MAFscore)].count()['MAF'])
			if custom_file != None:
				total_custom_captured.append(len(set(list(freq_file[freq_file['MAF'] >= float(MAFscore)]['SNP'])).intersection(custom_names_only)))
		str_MAF_iterables = [str(x) for x in MAF_iterables]
		str_MAF_iterables.append('NA')
		total_snps.append(freq_file['MAF'].isnull().sum())
		percentage_total_snps = [float(x)/float(len(freq_file.index)) for x in total_snps]
		
		# graph is generated only if custom file is provided
		if custom_file != None:
			percentage_custom_captured = [float(i)/float(len(custom_snp_dict)) for i in total_custom_captured]
			custom_null=len(set(list(freq_file[freq_file['MAF'].isnull()]['SNP'])).intersection(custom_names_only))
			percentage_custom_null=len(set(list(freq_file[freq_file['MAF'].isnull()]['SNP'])).intersection(custom_names_only))/float(len(custom_names_only))
			total_custom_captured.append(custom_null)
			percentage_custom_captured.append(percentage_custom_null)
			
			# formatting data for easy seaborn plotting
			custom_temp_array = [[str(str_MAF_iterables[i]), 'custom', percentage_custom_captured[i]] for i in range(0, len(str_MAF_iterables))]
			all_temp_arrary = [[str(str_MAF_iterables[i]), 'all', percentage_total_snps[i]] for i in range(0, len(str_MAF_iterables))] 

			# graph percentage retained between custom and all snps in chip		
			percentage_maf_comparison = pandas.DataFrame(np.matrix(custom_temp_array+all_temp_arrary), columns=['MAF_threshold', 'content', 'percent_retained'])
			percentage_maf_comparison[['percent_retained']] = percentage_maf_comparison[['percent_retained']].astype(float)
			plt.figure()
			maf_comparison = sns.barplot(x='MAF_threshold', y='percent_retained', hue='content', data=percentage_maf_comparison)
			maf_comparison.set(xlabel='MAF threshold filtered', ylabel='percentrage of SNPs retained')
			sns.plt.title('Total percentage of SNPs remaining at each MAF threshold')
			plt.tight_layout()
			pdf.savefig()
			plt.close()
		

		# no custom file required here; answers the question "How many SNPs remain at each MAF threshold after filtering by MAF"	
		plt.figure()
		plt.bar(np.arange(len(str_MAF_iterables)), total_snps, align='center')
		plt.xticks(np.arange(len(str_MAF_iterables)), str_MAF_iterables)
		plt.title('Total Number of SNPs Retained by MAF threshold', fontsize=20)
		plt.xlabel('MAF cut-off', fontsize=15)
		plt.ylabel('Total number of SNPs retained', fontsize=15)
		plt.tight_layout
		pdf.savefig()
		plt.close()


	def hwe():
		# distribution of SNP counts by p-value of exact test of HWE
		hwe_file = pandas.read_table(outputPrefix + '.hwe', delim_whitespace=True)
		plt.figure()
		hwe_pval_dist = hwe_file['P'].plot.hist(bins=20)
		plt.xlabel("p-values")
		plt.ylabel("Total Numer of SNPs")
		plt.title("Distribution of HWE p-values across all SNPs", fontsize=15)
		plt.tight_layout()
		pdf.savefig()
		plt.close()

		# plots total number of snps that are less than the max p-value threshold as each interval 
		total_snps_pvalues = []
		total_custom_pvalues = []
		pval_binning = [0.05, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
		for i in pval_binning:
			total_snps_pvalues.append(hwe_file[hwe_file['P'] < float(i)].count()['P'])
			if custom_file != None: # only if user supplies custom snp file
				total_custom_pvalues.append(len(set(list(hwe_file[hwe_file['P'] < float(i)]['SNP'])).intersection(custom_names_only)))
		
		# plots both all and custom snps that are less than max p-value threshold at each interval of Exact test HWE
		if len(total_custom_pvalues) > 0:
			sig_pvals = pandas.DataFrame(np.transpose(np.array([total_snps_pvalues])), index=pval_binning, columns=['all'])
			plt.figure()
			barplot = sig_pvals.plot(kind='bar', color=['Y'])
			# annotates the bar graphs by adding the total at each bar
			for p in barplot.patches:
				barplot.annotate(str(int(p.get_height())), (p.get_x()+p.get_width()/2.0, p.get_height() * 1.005), ha='center')
			plt.xlabel("HWE p-value threshold")
			plt.ylabel("Total SNPs")
			plt.title("Distribution of HWE p-values across all SNPs", fontsize=15)
			plt.tight_layout()
			pdf.savefig()
			plt.close()

			# plots the custom; not on same graph as all due to scaling of custom may be substantially different than all since counts not percents
			sig_pvals = pandas.DataFrame(np.transpose(np.array([total_custom_pvalues])), index=pval_binning, columns=['custom'])
			plt.figure()
			barplot = sig_pvals.plot(kind='bar', color=['G'])
			# annotates the bar graphs by adding the total at each bar
			for p in barplot.patches:
				barplot.annotate(str(int(p.get_height())), (p.get_x()+p.get_width()/2.0, p.get_height() * 1.005), ha='center')
			plt.xlabel("HWE p-value threshold")
			plt.ylabel("Total SNPs")
			plt.title("Distribution of HWE p-values across custom SNPs", fontsize=15)
			plt.tight_layout()
			pdf.savefig()
			plt.close()
			
		else:
			# only plots all, this is the case when user DOES NOT supply a custom snp list
			sig_pvals = pandas.DataFrame(np.transpose(np.array([total_snps_pvalues])), index=pval_binning, columns=['all'])
			plt.figure()
			sig_pvals.plot(kind='bar', color=['Y'])
			# annotates the bar graphs by adding the total at each bar
			for p in barplot.patches:
				barplot.annotate(str(int(p.get_height())), (p.get_x()+p.get_width()/2.0, p.get_height() * 1.005), ha='center')
			plt.xlabel("HWE p-value threshold")
			plt.ylabel("Total SNPs")
			plt.title("Distribution of HWE p-values across all SNPs", fontsize=15)
			plt.tight_layout()
			pdf.savefig()
			plt.close()

	
	def missing():
		# missing genotypes call percent by sample
		imiss_file = pandas.read_table(outputPrefix + '.imiss', delim_whitespace=True)
		stdev_imiss_dataset =stats.stdev(imiss_file['F_MISS'])
		plt.figure()
		outlier_list = imiss_file[imiss_file['F_MISS'] > 1.5*stdev_imiss_dataset]['IID']
		plt.boxplot(list(imiss_file['F_MISS']), 1, sym ='b.', showfliers=False)
		plt.title("Distribution of missing call rate by sample", fontsize=20)
		plt.ylabel("percentage snps missing")
		plt.xticks([1], ['all samples exluding outliers'])
		plt.figtext(0.70, 0.50, 'outliers:'+'\n'+'\n'.join(list(outlier_list)), color='black', backgroundcolor='wheat',
            weight='roman', size='small')
		plt.tight_layout()
		pdf.savefig()
		plt.close()
		
		# outlier statistic and distribution
		outlier_statistics = imiss_file[imiss_file['F_MISS'] > 1.5*stdev_imiss_dataset]
		# PUT THIS IN basic-statistics.txt NOT PDF!!!!!!
		#plt.figure()
		#plt.title("Overall missing callrate in across samples")
		#plt.figtext(0.20, 0.20, 'General Stats:'+'\n'+str(imiss_file.describe()['F_MISS'])+'\n\n'+'Outlier Stats:'+'\n'+str(outlier_statistics.describe()['F_MISS']), color='black', backgroundcolor='wheat',
            weight='roman', size='small')
		#pdf.savefig()
		#plt.close()
		
		
		total_snps_missing = []
		total_custom_missing = []
		# missing genotypes call percent by SNPs
		lmiss_file = pandas.read_table(outputPrefix + '.lmiss', delim_whitespace=True)
		group_by_percent_missing = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
		for miss_thresh in group_by_percent_missing:
			total_snps_missing.append(lmiss_file[lmiss_file['F_MISS'] >= miss_thresh].count()['F_MISS'])
			if custom_file != None:
				total_custom_missing.append(len(set(list(lmiss_file[lmiss_file['F_MISS']  >= miss_thresh]['SNP'])).intersection(custom_names_only)))

		# calculate percentages instead of frequency and then reformat lists for seaborn compatability
		reformat_all_seaborn = [[str(group_by_percent_missing[x]), 'all', str(total_snps_missing[x]/float(len(lmiss_file.index)))] for x in range(0, len(group_by_percent_missing))]
		if custom_file != None:
			# calculate percentages instead of frequency and then reformat lists for seaborn compatability for custom file
			reformat_custom_seaborn = [[str(group_by_percent_missing[i]), 'custom', str(total_custom_missing[i]/float(len(custom_names_only)))] for i in range(0, len(group_by_percent_missing))]
			merge = reformat_all_seaborn + reformat_custom_seaborn
			missing_snps_seaborn_datastruc = pandas.DataFrame(np.array(merge), columns=['missing_call_rate_threshold', 'content', "percentage_snps_remaining"])
			missing_snps_seaborn_datastruc[['percentage_snps_remaining']] = missing_snps_seaborn_datastruc[['percentage_snps_remaining']].astype(float)
			plt.figure()
			snp_callrate_comparison = sns.barplot(x='missing_call_rate_threshold', y='percentage_snps_remaining', hue='content', data=missing_snps_seaborn_datastruc)
			snp_callrate_comparison.set(xlabel='Percentage filtered missing call rate threshold', ylabel='percentage snps remaining')
			sns.plt.title('Percentage of SNPs remaining after filtering missing call rate of SNPs')
			plt.tight_layout()
			pdf.savefig()
			plt.close()

		else:
			missing_snps_seaborn_datastruc = pandas.DataFrame(np.array(reformat_all_seaborn), columns=['missing_call_rate_threshold', 'data', "percentage_snps_remaining"])
			missing_snps_seaborn_datastruc[['percentage_snps_remaining']] = missing_snps_seaborn_datastruc[['percentage_snps_remaining']].astype(float)
			snp_callrate_comparison = sns.barplot(x='missing_call_rate_threshold', y='percentage_snps_remaining', hue='content', data=missing_snps_seaborn_datastruc)
			snp_callrate_comparison.set(xlabel='Percentage filtered missing call rate threshold', ylabel='percentage snps remaining')
			sns.plt.title('Percentage of SNPs remaining after filtering missing call rate of SNPs')
			plt.tight_layout()
			pdf.savefig()
			plt.close()


	def mendel():
		pass;

	
	def maf_with_missing():
		pass;

	
	with PdfPages(str(outputPrefix)+'-basic-analysis.pdf') as pdf:
		maf_analysis()
		hwe()
		missing()








if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-plinkEx', required=True, dest='plinkPath', type=str, help='Full path to PLINK executable')
	parser.add_argument('-input', required=True, dest='inputFile', type=str, help='Full path to input file')
	parser.add_argument('-mode', default='both', dest='mode_type', type=str, help='options: both, stats_only, analysis_only')
	parser.add_argument('-out', default=None, dest='outputName', type=str, help='Full path to and preferred prefix of output files')
	parser.add_argument('--custom', default=None, dest='custom', type=str, help='Full path to custom SNP file (.csv), one SNP name per line (Name,CHR,POS)')
	args=parser.parse_args()
	
	if args.outputName == None:
		args.outputName = args.inputFile[:-4]
		print args.outputName

	#if args.mode_type == 'stats_only':
	#	print "Running stats_only mode..."
	#	basic_stats(programLoc=args.plinkPath, inputFile=args.inputFile)
	#elif args.mode_type == 'analysis_only':
	#	print "Running analysis_only mode..."
	#	visualize_stats()
	#else:
	#	print "Running both mode..."
	#basic_stats(programLoc=args.plinkPath, inputFile=args.inputFile, outputPrefix=args.outputName)
	visualize_stats(outputPrefix=args.outputName, custom_file=args.custom)

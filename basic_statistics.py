import argparse
import math
import subprocess
import matplotlib.pyplot as plt
import matplotlib.collections as collections
import seaborn as sns
import pandas
import numpy as np
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
	
	# read in custom_file if user supplies one
	if custom_file != None:
		custom_snp_dict = {}
		with open(custom_file) as custom_snp_file:
			for line in custom_snp_file:
				snp, chrm, pos = line.rstrip().split(',')
				custom_snp_dict[snp] = [chrm, pos]

		custom_names_only = [key for key in custom_snp_dict]
	
	
	def maf_analysis():
		# distribution of SNP counts by MAF
		freq_file = pandas.read_table(outputPrefix + '.frq', delim_whitespace=True)
		plt.figure()
		maf_dist = freq_file['MAF'].plot.hist(bins=10)
		plt.xlabel("MAF scores")
		plt.ylabel("Total Numer of SNPs")
		plt.title("Distribution of MAF across all SNPs", fontsize=15)
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
			pdf.savefig()
			plt.close()
		

		# no custom file required here; answers the question "How many SNPs remain at each MAF threshold after filtering by MAF"	
		plt.figure()
		plt.bar(np.arange(len(str_MAF_iterables)), total_snps, align='center')
		plt.xticks(np.arange(len(str_MAF_iterables)), str_MAF_iterables)
		plt.title('Total Number of SNPs Retained by MAF threshold', fontsize=20)
		plt.xlabel('MAF cut-off', fontsize=15)
		plt.ylabel('Total number of SNPs retained', fontsize=15)
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
			pdf.savefig()
			plt.close()

	
	def missing():
		pass;

	def mendel():
		pass;

	with PdfPages(str(outputPrefix)+'-basic-analysis.pdf') as pdf:
		maf_analysis()
		hwe()








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

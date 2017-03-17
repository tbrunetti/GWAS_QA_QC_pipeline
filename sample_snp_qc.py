from chunkypipes.components import *
import os
import pandas
import matplotlib.pyplot as plt
import matplotlib.collections as collections


class Pipeline(BasePipeline):
	def dependencies(self):
		# assuming user as pip installed
		return ['pandas', 'matplotlib']

	def description(self):
		return 'Pipeline to perform sample QC (call rate, HWE, Mendelian Error, LD-pruning)'

	def configure(self):
		return {
			'plink':{
				'path': 'Full path to PLINK executable (must be version >=1.9):'
			}
		}

	def add_pipeline_args(self, parser):
		parser.add_argument('-fileLocation', required=True, help='Full path to directory in which BINARY PLINK files are located including sample name prefix')
		parser.add_argument('--maf', default=0.05, help="maximum maf cutoff threshold, i.e. remove variants BELOW this value")
		parser.add_argument('--sample_missing_callrate', default=0.05, help='maximum percent (0.0-1.0) missing call rate in a SAMPLE in order to be removed, i.e. filter out samples with missing call rate EXCEEDING this threshold')
		parser.add_argument('--snp_missing_callrate', default=0.05, help='maximum percent (in decimal form 0.0-1.0) missing call rate in a SNP in order to be removed, i.e. filter out variants that EXCEED this value')
		parser.add_argument('--hwe_cutoff', default=0.0001, help='filter below p-value cutoff for filtering Hardy-Weinberg equilibrium exact test value. i.e. remove variants that have a p-value SMALLER than this value')
		parser.add_argument('--windowSize', default=50, help='window size in kb to be used for LD pruning')
		parser.add_argument('--varStep', default=5, help='variant step side to slide window for LD pruning')
		parser.add_argument('--r2_thresh', default=0.50, help='pairwise r-squared threhold at each step for each variant pair, if correlation is greater than set thresh, var is pruned')
	
	def run_pipeline(self, pipeline_args, pipeline_config):
		plink = Software('plink', pipeline_config['plink']['path'])

		# ('/', 1)  is number of splits starting from right, so will split one time and the first time will be the first / from the right
		fileLoc, filePrefix = pipeline_args['fileLocation'].rsplit('/', 1)

		# missing call rate, sample and snp level
		plink.run(
			Parameter('--bfile', pipeline_args['fileLocation']),
			Parameter('--geno', pipeline_args['snp_missing_callrate']),
			Parameter('--mind', pipeline_args['sample_missing_callrate']),
			Parameter('--maf', pipeline_args['maf']),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['fileLocation']+'_callrate_cleanup')
			)

		# Hardy-Weinberg Equilibrium test
		# eventually this will need to be run by ethnicity
		plink.run(
			Parameter('--bfile', pipeline_args['fileLocation']+'_callrate_cleanup'),
			Parameter('--hwe', pipeline_args['hwe_cutoff']),
			Parameter('midp'), # applies mid-p adjustment to reduce favoritism of variants with missing data
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['fileLocation']+'_callrate_cleanup_HWEfiter')
			)


		# first need to calculate the per-SNP (.lmendel) and per-individual error rates (.imendel)
		# TODO: 
		# //run statistics on mendel and HWE in determine relatedness to get an idea of what values to pass into filters???

		#Mendelian Error rate filter (need to look up general error rates to set as default)
		#plink.run(
		#		Parameter('--bfile', pipeline_args['fileLocation']+'_callrate_cleanup_HWEfiter'),
		#		Parameter('--me'),
		#		Parameter('--make-bed')
		#		)

		# LD-pruning via plink (remember to add back mendel error output to --file)
		plink.run(
			Parameter('--bfile', pipeline_args['fileLocation']+'_callrate_cleanup_HWEfiter'),
			Parameter('--indep-pairwise', pipeline_args['windowSize']),
			Parameter(pipeline_args['varStep']),
			Parameter(pipeline_args['r2_thresh']),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['fileLocation']+'_callrate_cleanup_HWEfiter_mendelErrorfilter_LDpruned')
			)


#het from plink to determine heterozygosity rates and inbreeding coefficeints (LD prune first!)s
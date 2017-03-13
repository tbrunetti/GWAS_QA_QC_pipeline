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
		return 'Pipeline to perform sample QC (call rate, pruning, HWE, Mendelian Error)'

	def configure(self):
		return {
			'plink':{
				'path': 'Full path to PLINK executable (must be version >=1.9):'
			}
		}

	def add_pipeline_args(self, parser):
		parser.add_argument('-fileLocation', required=True, help='Full path to directory in which BINARY PLINK files are located including sample name prefix')
		parser.add_argument('--sample_missing_callrate', default=0.05, help='maximum percent (0.0-1.0) missing call rate in a SAMPLE in order to be removed')
		parser.add_argument('--snp_missing_callrate', default=0.05, help='maximum percent (in decimal form 0.0-1.0) missing call rate in a SNP in order to be removed')
		parser.add_argument('--hwe_cutoff', default=0.0001, help='p-value cutoff for filtering Hardy-Weinberg equilibrium exact test value')
		parser.add_argument('--windowSize', default=50, help='window size in kb to be used for LD pruning')
		parser.add_argument('--varStep', default=5, help='variant step side to slide window for LD pruning')
		parser.add_argument('--r2_thresh', default=0.50, help='pairwise r-squared threhold at each step for each variant pair, if correlation is greather than set thresh, var is pruned')
	
	def run_pipeline(self, pipeline_args, pipeline_config):
		plink = Software('plink', pipeline_config['plink']['path'])

		# missing call rate, sample and snp level
		plink.run(
			Parameter('--bfile', pipeline_args['fileLocation']),
			Parameter('--geno', pipeline_args['sample_missing_callrate']),
			Parameter('--mind', pipeline_args['snp_missing_callrate']),
			Parameter('--make-bed'),
			Parameter('--out', pipeline_args['fileLocation']+'_callrate_cleanups')
			)

		# Hardy-Weinberg Equilibrium test
		# eventually this will need to be run by ethnicity
		plink.run(
			Parameter('--bfile', pipeline_args['fileLocation']+'_callrate_cleanup'),
			Parameter('--hwe', pipeline_args['hwe_cutoff']),
			Parameter('midp'), # applies mid-p adjustment to reduce favoritism of variants with missing data
			Parameter('--out', pipeline_args['fileLocation']+'_callrate_cleanup_HWEfiter')
			)

		#Mendelian Error rate (need to look up general error rates to set as default)
		#plink.run(
		#		Parameter('--bfile', pipeline_args['fileLocation']+'_callrate_cleanup_HWEfiter'),
		#	Parameter('--me')
		#	)

		# LD-pruning via plink
		plink.run(
			Parameter('--bfile', pipeline_args['fileLocation']+'_callrate_cleanup_HWEfiter_mendelErrorfilter'),
			Parameter('--indep-pairwise', pipeline_args['windowSize']),
			Parameter(pipeline_args['varStep']),
			Parameter(pipeline_args['r2_thresh']),
			Parameter('--out', pipeline_args['fileLocation']+'_callrate_cleanup_HWEfiter_mendelErrorfilter_LDpruned')
			)

# window size (kb), variant step size, r^2 threshold
#./../TOOLS/plink --bfile CCPM-Ex_validation_02-20-17-related-removed-snp-and-sample-callRate-filter-applied --indep-pairwise 50 5 0.5 --out LDpruned

#het from plink to determine heterozygosity rates and inbreeding coefficeints (LD prune first!)
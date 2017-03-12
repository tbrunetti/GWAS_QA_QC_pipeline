from chunkypipes.components import *
import os
import pandas
import matplotlib.pyplot as plt
import matplotlib.collections as collections

class Pipeline(BasePipeline):
	def dependencies(self):
		return ['pandas']

	def description(self):
		return 'This pipeline assesses relatedness in a set of genoptying samples'

	def configure(self):
		return {
			'plink':{
				'path': 'Full path to plink executable:'
			},
			'king':{
				'path': 'Full path to KING executable:'
			}
		}

	def add_pipeline_args(self, parser):
		parser.add_argument('-location', required=True, help='Full path to directory in which PLINK files are located')
		parser.add_argument('-input', required=True, help='sample prefix name of PLINK file to be used for analysis')
		parser.add_argument('--degree', default=2, help='KING cutoff for degree of relatedness to remove')
		parser.add_argument('--maf', default=0.05, help='MAF cutoff for binary plink files')

	def run_pipeline(self, pipeline_args, pipeline_config):
		# implement software
		plink = Software('plink', pipeline_config['plink']['path'])
		king = Software('king', pipeline_config['king']['path'])

		# check if plink binaries exist, if not convert plink .ped and .map to binary format
		os.chdir(pipeline_args[location])
		
		if pipeline_args['input']+'.bed' in os.listdir(pipeline_args[location]):
			continue;	
		else:
			plink.run(
				Parameter('--file', pipeline_args['input']),
				Parameter('--maf', pipeline_args['maf']),
				Parameter('--make-bed'),
				Parameter('--out', pipeline_args['input'])
			)
		
		# run king three times to get all information
		king.run(
			Parameter('-b', pipeline_args['input']+'.bed'),
			Parameter('--kinship'),
			Parameter('--prefix', pipeline_args['input'])
			)

		king.run(
			Parameter('-b', pipeline_args['input']+'.bed'),
			Parameter('--kinship'),
			Parameter('--ibs'),
			Parameter('--prefix', pipeline_args['input'])
			)

		king.run(
			Parameter('-b', pipeline_args['input']+'.bed'),
			Parameter('--unrelated'),
			Parameter('--degree', pipeline_args['degree']),
			Parameter('--prefix', pipeline_args['input'])
			)

		# graph output of king using pandas
		# between-family relationships
		#fig = plt.figure()
		#ax = fig.add_subplot(111)
		king_kinship_matrix = pandas.read_table(pipeline_args['location'] + '/' + pipeline_args['input']+'.kin0')
		graph = king_kinship_matrix.plot(kind='scatter', x='Kinship', y='IBS0')
		# kinship coeff > 0.354 duplicates/monozygotic twins
		# 0.177 < kinship coeff < 0.354 1st degree (parent-offspring, siblings)
		# 0.0884 < kinship coeff < 0.177 2nd degree (half-sibs, grandparent-grandchild) 
		dup_twins = collections.BrokenBarHCollection([(0.354, 0.6-0.354)], (0, 0.20), facecolor='red', alpha=0.3)
		sib_par = collections.BrokenBarHCollection([(0.177, 0.354-0.177)], (0, 0.20), facecolor='yellow', alpha=0.3)
		half_grand = collections.BrokenBarHCollection([(0.0884, 0.177-0.0884)], (0, 0.20), facecolor='green', alpha=0.3)
		graph.add_collection(dup_twins)
		graph.add_collection(sib_par)
		graph.add_collection(half_grand)
		plt.title('Kinship between families')
		plt.axvline(x=0.354, linestyle='--')
		plt.axvline(x=0.177, linestyle='--')
		plt.axvline(x=0.0884, linestyle='--')
		plt.show()

		# TODO:  remove related indivuals
		#./../TOOLS/plink --file CCPM-Ex_validation_02-20-17 --keep kingunrelated.txt --make-bed --out CCPM-Ex_validation_02-20-17-related-removed
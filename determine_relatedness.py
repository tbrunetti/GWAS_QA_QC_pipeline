from chunkypipes.components import *
import os
import pandas
import matplotlib.pyplot as plt

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
		king_kinship_matrix = pandas.read_table(pipeline_args['location'] + '/' + pipeline_args['input']+'.kin0')
		graph = king_kinship_matrix.plot(kind='scatter', x='Kinship', y='IBS0')
		# kinship coeff > 0.354 duplicates/monozygotic twins
		# 0.177 < kinship coeff < 0.354 1st degree (parent-offspring, siblings)
		# 0.0884 < kinship coeff < 0.177 2nd degree (half-sibs, grandparent-grandchild) 
		plt.axvline(x=0.354)
		plt.axvline(x=0.177)
		plt.axvline(x=0.0884)
		plt.show()
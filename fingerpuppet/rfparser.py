#!/usr/bin/env python

import fingerpuppet
import tables as tbs
import fingerpuppet.ptstructs as pts


# This submodule implements the parsing of UCSC Genome Browser refFlat.txt files into PyTables



class refFlat_parser():

	
	def __init__(self, sourcefile, birthday, zoodir, halfagiraffe, species):
	
		self.source_file = sourcefile
		self.birthday = creationdate
		self.zoo_dir = zoodir
		self.half_a_giraffe = halfagiraffe
		self.species = species
		self.gene_symbols = []
		self.trans_ids = []
		self.gene_data = {}
		self.trans_data = {}
		self.exons = {}
		
		
	def parse_refFlat(self):
	
		import sys,os
		
		refFlat = open(self.source_file, 'r')
		
		for trans in refFlat.readlines():
		
			curdata = trans.strip().split("\t")
			
			curgenesymbol = curdata[0]
			curtransid = curdata[1]
			curcontig = curdata[2]
			curstrand = curdata[3]
			curtxstart = int(curdata[4])
			curtxend = int(curdata[5])
			curcdstart = int(curdata[6])
			curcdend = int(curdata[7])
			
			numexons = int(curdata[8])
			exonstarts = curdata[9].strip(",").split(",")
			exonends = curdata[10].strip(",").split(",")
			curexons = []
			for exon in range(numexons):
				curexons.append((int(exonstarts[exon]),int(exonends[exon])))
				
			if not curgenesymbol in self.gene_symbols:
				self.gene_symbols.append(curgenesymbol)
			if not curtransid in self.trans_ids:
				self.trans_ids.append(curtransid)
				
			if not curgenesymbol in self.gene_data.keys():
				self.gene_data[curgenesymbol] = {}
				self.gene_data[curgenesymbol]["gene_symbol"] = curgenesymbol
				self.gene_data[curgenesymbol]["trans_ids"] = [curtransid]
				self.gene_data[curgenesymbol]["contig"] = curcontig
				self.gene_data[curgenesymbol]["strand"] = curstrand
				self.gene_data[curgenesymbol]["trans_starts"] = [curtxstart]
				self.gene_data[curgenesymbol]["trans_ends"] = [curtxend]
				self.gene_data[curgenesymbol]["ccds_starts"] = [curcdstart]
				self.gene_data[curgenesymbol]["ccds_ends"] = [curcdend]
				self.gene_data[curgenesymbol]["exons"] = curexons
			else:
				if not curtransid in self.gene_data[curgenesymbol]["trans_ids"]:
					self.gene_data[curgenesymbol]["trans_ids"].append(curtransid)
				if not curtxstart in self.gene_data[curgenesymbol]["trans_starts"]:
					self.gene_data[curgenesymbol]["trans_starts"].append(curtxstart)
				if not curtxend in self.gene_data[curgenesymbol]["trans_ends"]:
					self.gene_data[curgenesymbol]["trans_ends"].append(curtxend)
				if not curcdstart in self.gene_data[curgenesymbol]["ccds_starts"]:
					self.gene_data[curgenesymbol]["trans_starts"].append(curcdstart)
				if not curcdend in self.gene_data[curgenesymbol]["ccds_ends"]:
					self.gene_data[curgenesymbol]["trans_ends"].append(curcdend)
				for exon in curexons:
					if not exon in self.gene_data[curgenesymbol]["exons"]:
						self.gene_data[curgenesymbol]["exons"].append(exon)
						
			if not curtransid in self.trans_data.keys():
				self.trans_data[curtransid] = {}
				self.trans_data[curtransid]["trans_id"] = curtransid
				self.trans_data[curtransid]["trans_nums"] = [0]
				self.trans_data[curtransid]["gene_symbol"] = curgenesymbol
				self.trans_data[curtransid]["contig"] = curcontig
				self.trans_data[curtransid]["strand"] = curstrand
				self.trans_data[curtransid]["trans_structures"] = []
				self.trans_data[curtransid]["trans_structures"].append([curtxstart,curtxend,curcdstart,curcdend,curexons])
			else:
				self.trans_data[curtransid]["trans_nums"].append(len(self.trans_data[curtransid]["trans_nums"]))
				self.trans_data[curtransid]["trans_structures"].append([curtxstart,curtxend,curcdstart,curcdend,curexons])
				
			for exon in curexons:
				curexonlocus = curcontig + ":" + str(exon[0]) + "-" + str(exon[1]) + ":" + curstrand
				if not curexonlocus in self.exons.keys():
					self.exons[curexonlocus] = {}
					self.exons[curexonlocus]["gene_symbols"] = [curgenesymbol]
					self.exons[curexonlocus]["trans_ids"] = [curtransid]
				else:
					if not curgenesymbol in self.exons[curexonlocus]["gene_symbols"]:
						self.exons[curexonlocus]["gene_symbols"].append(curgenesymbol)
					if not curtransid in self.exons[curexonlocus]["trans_ids"]:
						self.exons[curexonlocus]["trans_ids"].append(curtransid)
						
						
	def initialize_half_a_giraffe(self):
	
		import sys,os
		import tables
		import fingerpuppet.ptstructs as pts
		
		if not os.path(self.zoo_dir):
			os.mkdir(self.zoo_dir)
		if not os.exists(self.zoo_dir + "/" + self.half_a_giraffe):
			giraffe = tables.openFile(self.zoo_dir + "/" + self.half_a_giraffe, mode = "a", title = "fingerpuppet giraffe, born " + self.birthday)
		else:
			giraffe = tables.openFile(self.zoo_dir + "/" + self.half_a_giraffe, mode = "a")
		
		existing_nodes = giraffe.listNodes()
		if not 'summary_data' in existing_nodes:
			summary_data = giraffe.createGroup("/","summary_data","summary_data")
			structure_summary = giraffe.createTable(summary_data,'structure_summary',pts.detail_table,'structure_summary')
			window_summary = giraffe.createTable(summary_data,'window_summary',pts.detail_table,'window_summary')
			read_depth_summary = giraffe.createTable(summary_data,'read_depth_summary',pts.read_depth_summary_table,'read_depth_summary')
		else:
			summary_data = giraffe.summary_data
			structure_summary = giraffe.summary_data.structure_summary
			window_summary = giraffe.summary_data.window_summary
			read_depth_summary = giraffe.summary_data.read_depth_summary
		if not 'locus_data' in existing_nodes:
			locus_data = giraffe.createGroup("/","locus_data","locus_data")
			window_table = giraffe.createTable(locus_data,"window_table",pts.window_table,"window_table")
		else:
			locus_data = giraffe.locus_data
			window_table = giraffe.locus_data.window_table
		if not 'relationship_data' in existing_nodes:
			relationship_data = giraffe.createGroup("/","relationship_data","relationship_data")
			relationship_table = giraffe.createTable(relationship_data,"relationship_table",pts.relationship_table,"relationship_table")
		else:
			relationship_data = giraffe.relationship_data
			relationship_table = giraffe.relationship_data.relationship_table
		if not 'read_depth_data' in existing_nodes:
			read_depth_data = giraffe.createGroup("/","read_depth_data","read_depth_data")
		else:
			read_depth_data = giraffe.read_depth_data
			
		self.giraffe_legs = {"giraffe": giraffe, "summary_data": summary_data, "structure_summary": structure_summary, "window_summary": window_summary, "read_depth_summary": read_depth_summary, "locus_data": locus_data, "window_table": window_table, "relationship_table": relationship_data, "relationship_table": relationship_table, "read_depth_data": read_depth_data}
			
			
	def generate_relationships(self):
	
		 self.relationships = []
		 
		 for gene in self.gene_data.keys():
		 	for trans in self.gene_data[gene]["trans_ids"]:
		 		self.relationships.append([gene,"gene",trans,"transcript"])
		 	for ss in self.gene_data[gene]['trans_starts']:
		 		self.relationships.append([gene,"gene",ss,"trans_start"])
		 	for end in self.gene_data[gene]['trans_ends']:
		 		self.relationships.append([gene,"gene",end,"trans_end"])
		 	for ss in self.gene_data[gene]['ccds_starts']:
		 		self.relationships.append([gene,"gene",ss,"coding_start"])
		 	for end in self.gene_data[gene]['ccds_ends']:
		 		self.relationships.append([gene,"gene",end,"coding_end"])
		 	for exon in self.gene_data[gene]["exons"]:
		 		self.relationships.append([gene,"gene",self.gene_data[gene]["contig"] + ":" + str(exon[0]) + "-" + str(exon[1]) + ":" + self.gene_data[gene]["strand"],"exon"])
		 
		 for trans in self.trans_data.keys():
		 	for struc in self.trans_data[trans]["trans_structures"]:
		 		self.relationships.append([trans,"transcript",struct[0],"trans_start"])
		 		self.relationships.append([trans,"transcript",struct[1],"trans_end"])
		 		self.relationships.append([trans,"transcript",struct[2],"coding_start"])
		 		self.relationships.append([trans,"transcript",struct[3],"coding_end"])
		 		for exon in struct[4]:
		 			self.relationships.append([trans,"transcript",self.trans_data[trans]["contig"] + ":" + str(exon[0]) + "-" + str(exon[1]) + ":" + self.trans_data[trans]["strand"],"exon"])
	
	
	def make_a_baby_giraffe(self):
	
		curtable = self.giraffe_limbs["giraffe"].createTable(self.giraffe_limbs["locus_data"],"all_structures",pts.structure_table,"all_structures")
		currows = curtable.rows
		for gene in self.gene_data.keys():
			currows['id'] = self.gene_data[gene]['gene_symbol'][
		 	
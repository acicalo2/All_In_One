from Bio import SearchIO
from Bio import SeqIO
from ete3 import NCBITaxa

import sqlite3
import sys
import os
import configparser
import argparse
import re

parser = argparse.ArgumentParser(
                    prog='VS_MD_parser',
					description='Parse blast table',
					epilog='Parse the blast table from sequence comparative algorithm (MMSeqs,blastn,blastx)')
# add options
parser.add_argument("-i", "--input")
parser.add_argument("-t", "--type")
parser.add_argument("-r", "--rank")
parser.add_argument("-o", "--outdir")


args = parser.parse_args()
print(args.input,args.type,args.rank, args.outdir)
if (args.input == None):
        print(parser.usage)
        exit(0)

bdir = os.path.dirname(args.input)

#print(bdir)
input_file = args.input
sample_id = os.path.basename(input_file)
sample_id = sample_id.replace("_blastx_diamondview.tsv","")
sample_id = re.sub(r"_mmseq_.*contig_against_NT.out","", sample_id)

"_mmseq_dragonflye_contig_against_NT.out"
outdir = args.outdir
if (args.outdir == None):
	outdir = "."

#print(sample_id)
# Get db paths from config file
config = configparser.ConfigParser()
config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'VS.cfg')
config.read(config_path)
paths = config['paths']
vhunter = paths['vhunter']
ncbidb = paths['ncbi_taxadb']
ncbi = NCBITaxa(dbfile = ncbidb)
rank = args.rank
out_rank = rank if rank else "allRanks"
# open a connection to database
try:
	connector = sqlite3.connect(vhunter)
except:
	print("Cannot connect to DB: "+vhunter)

cursor_dbhsqlite = connector.cursor()

def read_FASTA_data(fastaFile):
	fa_dict = SeqIO.index(fastaFile, "fasta")
	return fa_dict
# function to determine the taxonomy lineage for a given blast hit
def PhyloType(lineage_ref, result_ref, hit_ref, COUNT, lin_count, assignment_ref):
	assigned = 1
	description = ""
	lineage = ""
	#This for loop basically just grabs the scientific name for all the taxids in the lineage and saves it to a single variable
	for temp_node_id in lineage_ref:
		temp_name = ncbi.get_taxid_translator([temp_node_id])
		lineage += temp_name[temp_node_id]+";"
		#print(lineage)
	# Is there already a significant representative blast hit for this lineage?
	if lineage in assignment_ref.keys():
		# If yes check to see if the same accession (this changes the information for aln_apan)
		ref_list = assignment_ref[lineage].split("\t")
		#print("Hit ref accession: "+hit_ref.accession)
		#print("Dictionary value: " + assignment_ref[lineage])
		if hit_ref.accession in assignment_ref[lineage]:
			#print("TEST!")
			old_aln_length = int(ref_list[5])
			curr_aln_length = hit_ref.hsps[0].aln_span
			new_target_covered = ((curr_aln_length+old_aln_length) / hit_ref.seq_len)*100
			ref_list[8] = str(new_target_covered)
		lin_count +=1
		#print(ref_list)
		ref_list[-1] = str(lin_count)
		lineage_desc = "\t".join(ref_list)
		assignment_ref[lineage] = lineage_desc
	else:
		target_covered = (hit.hsps[0].aln_span / hit.seq_len)*100
		target_span = "["+str(hit_ref.hsps[0].hit_start)+"\t"+str(hit_ref.hsps[0].hit_end)+"]"
		lineage_desc = "\t".join([lineage,temp_name[temp_node_id],hit_ref.accession,str(hit_ref.seq_len),hit_ref.description,str(hit_ref.hsps[0].aln_span),str(hit_ref.hsps[0].ident_pct)+"%",target_span,str(target_covered)+"%",str(hit_ref.hsps[0].evalue),str(lin_count)])
		assignment_ref[lineage] = lineage_desc
	return assigned,assignment_ref, lin_count

# function that returns description for assignment dictionary based on the input parsing file
def get_long_desc(first_col):
	target_covered = (hit.hsps[0].aln_span / hit.seq_len)*100
	target_span = "["+str(hit.hsps[0].hit_start)+"\t"+str(hit.hsps[0].hit_end)+"]"
	description = "\t".join([first_col+hit.accession,str(hit.seq_len),hit.description,str(hit.hsps[0].aln_span),str(hit.hsps[0].ident_pct)+"%",target_span,str(target_covered)+"%",str(hit.hsps[0].evalue),"N/A"])
	return description

def filter_lineage(curr_lineage):
	names = ncbi.get_taxid_translator(curr_lineage)
	lineage2ranks = ncbi.get_rank(names)
	for key, value in lineage2ranks.items():
		if value == rank:
			desired_rank = key
			break
		else:
			desired_rank = "NA"
		#print(species)
	if desired_rank != "NA":
		curr_lineage = curr_lineage[:curr_lineage.index(desired_rank) + 1]
	return curr_lineage


###################################################################################
if args.type == "mmseqs":
    outFile = outdir+"/"+sample_id+"."+out_rank+"_contig.mmseqs.parsed"
    e_cutoff = 1e-10
    BX=False
    MB=True
elif args.type == "blastx":
    outFile = outdir+"/"+sample_id+"."+out_rank+".blastx.parsed"
    e_cutoff = 1e-3
    BX=True
    MB=False
# Create output file
try:
    out = open(outFile, 'w')
except IOError:
	print("can not open file "+outFile)
	sys.exit(1)	


print("parsing blast output files...\n\n")

custom_fields=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","sallgi","qlen","slen"]
report = SearchIO.parse(input_file, "blast-tab", fields=custom_fields)

out.write("sample_ID\tquery_ID\tquery_length\tlineage\tlowest_classification\ttref_acc\tref_length\tref_description\talignment_length\tpercent_id\ttarget_start\ttarget_end\t%_ref_covered\tevalue\tnum_alignments\n")

# ref_length - target_span

# Go through BLAST reports one by one
for result in report:
	#print(result.__dir__())
	haveHit = 0
	assignment = {}
	assignment_ref = {}
	assignment_NotBestE = {}
	# only take the best hits
	best_e = 100
	hit_count = 0
	determined = 0
	lineage_count = 1

	for hit in result:
		#hit_desc = hit.id
		#if args.type == "mmseqs":
		#	hit.description = hit_desc.split(" ", 1)[1]
		#else:
		#	hit.description = hit_desc
		#print(hit.__dir__())
		# from hit name get hit gi number
		hit.description = hit.gi_all
		hit_name = hit.id
		hit.accession = hit.id
		print(hit_name)
		haveHit = 1
		hit_count+=1
		eval = hit[0].evalue
		if hit_count == 1:
			best_e = hit[0].evalue
			best_bit = hit[0].bitscore
		print(hit[0].bitscore)
		bit_score = hit[0].bitscore
		print(hit.hsps)
		# check whether the hit should be kept for further analysis
		if eval <= e_cutoff: # similar to known, need Phylotyped
			#if hit[0].evalue <= best_e: # only get best hits
			#if bit_score >= best_bit # only get best hits
				#get taxonomy lineage
			sth = cursor_dbhsqlite.execute("SELECT * FROM acc_taxid_prot where accession_version = '"+hit_name+"'") if BX else cursor_dbhsqlite.execute("SELECT * FROM acc_taxid_nucl where accession_version = '"+hit_name+"'")
			ref = sth.fetchone()
			if ref: # some gi don't have record in gi_taxid_nucl
				taxID = ref[2]
				taxon_name = ncbi.get_taxid_translator([taxID])
				if not taxon_name:
					description = get_long_desc("undefined taxon\t")
					assignment["other"] = description
				else:
					lineage = ncbi.get_lineage(taxID)
					if rank:
						lineage = filter_lineage(lineage)
					lineage.pop(0)
					if lineage:
						determined = 1
						success,assignment,lineage_count = PhyloType(lineage, result, hit, hit_count, lineage_count, assignment)

			else: # for situations that gi does not have corresponding taxid
				determined = 1
				description = get_long_desc("undefined taxon\t")
				assignment["other"] = description


		# finish phylotype for given hit
		elif BX: # e value is not significant. Only do this for the blastx input files
			if determined: # skip the rest hits that are not significant
				break
			else:
				target_covered = (hit.hsps[0].aln_span / hit.seq_len)*100
				target_span = "["+str(hit.hsps[0].hit_start)+"\t"+str(hit.hsps[0].hit_end)+"]"
				nosig_desc = "\t".join(["hit not significant",hit.accession,str(hit.seq_len),hit.description,str(hit.hsps[0].aln_span),str(hit.hsps[0].ident_pct)+"%",target_span,str(target_covered)+"%",str(hit.hsps[0].evalue)])
				assignment["unassigned"] = nosig_desc
				break

	if BX and not haveHit: # only do this for the blastx input files
		assignment["unassigned"] = "no hit"

	num_assignment = assignment.keys()

	#print(assignment.keys())
	#print(assignment.values())
	for assign in assignment.keys():
		out.write(sample_id+"\t"+result.id+"\t"+str(result.seq_len)+"\t"+assignment[assign]+"\n")
	for assign1 in assignment_NotBestE.keys():
		out.write(sample_id+"\t"+result.id+"\t"+str(result.seq_len)+"\t"+assignment_NotBestE[assign1]+"\n")

out.close()

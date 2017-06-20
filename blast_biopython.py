#!/usr/bin/python

"""
Le $PATH doit contenir le chemin suivant
/usr/local/bioinfo/src/NCBI_Blast+/ncbi-blast-2.2.29+/bin/
enregister cette ligen dans le /home/<user-name>/.profile
fpiumi@genotoul2 /work/fpiumi $ cat /home/fpiumi/.profile
export PATH=$PATH:/usr/local/bioinfo/src/NCBI_Blast+/ncbi-blast-2.2.29+/bin/

commande : 
blast_biopython -q parametre_1 -d parametre_2 -i parametre_3
prend en entree 3 parametres
parametre_1 : un fichier de sequences au format fasta a blaster (query)
parametre_2 des sequences contre lesquelles on va blaster (database) soit un fichier fasta donne par l'utilisateur soit un lien vers une databse de genotoul
parametre_3 : si les seuuences de la database sont deja indexees (yes) ou a indexer (comme dans le cas d'un fichier fasta = no)


ICI on ne cherche pas seulement le meilleur BLAST mais celui qui correspond a Bos taurus

outfmt=0 (pairwise) est remplace par outfmt=5 (XML Blast output) 
pour permettre le traitement de la sortie par NCBIXML.read
"""


from sys import argv
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from optparse import OptionParser
import subprocess
import shutil
import os
import sys
import re
from os import getcwd




parser = OptionParser()
parser.add_option('-q', '--query', dest='query', help='Query file')
parser.add_option('-d', '--database', dest='database', help='Database')
parser.add_option('-i', '--index', dest='index', help='Database_indexee')
parser.add_option('-s', '--species', dest='species', help='Species')
(options, args) = parser.parse_args()

if not options.query:
	options.query = raw_input('Enter Query file:')

if options.database is None:
	options.database = raw_input('Enter Database:')
elif options.database == "nt":
	chemin_indexed_database = "/bank/blastdb/nt"

if options.index is None:
	options.index = raw_input('Index yes or no:')

if options.species is None:
	options.species = raw_input('Enter Species:')



# remove old temporary directories
rep_cour = getcwd()
temp_directory = rep_cour+"/tmp/"
os.system("rm -rf tmp")
output_directory = rep_cour+"/output/"
os.system("rm -rf output")
indexed_database_directory = rep_cour+"/DBs/"
os.system("rm -rf DBs")


# create new temporary directories
rep_cour = getcwd()
temp_directory = rep_cour+"/tmp/"
os.system("mkdir tmp")
output_directory = rep_cour+"/output/"
os.system("mkdir output")
indexed_database_directory = rep_cour+"/DBs/"
os.system("mkdir DBs")



########################
# FUNCTIONS            #
##################################################################################
# FRACTIONNEMENT DU FICHIER "QUERY" EN SOUS FICHIERS DE SEQUENCES INDIVIDUELLES  #
##################################################################################

"""
decoupage du fichier fasta queries en fichiers ne contenant qu'une sequence
les fichiers sont places dans un repertoire temporaire
blastall

"""

def split_fasta (seqFile,temp_directory):
	# Split sequence file into as many fragments as appropriate depending on the size of original_fasta
	
	if not options.query:
		print "pas de fichier fasta"
		sys.exit()

	for seq_record in SeqIO.parse(seqFile,"fasta"):

		#print seq_record.id
		# seq_record.id est la premiere chaine de caractere trouvee du signe ">" jusqu'a un espace dans l'ID de la sequence

		#il ne faut pas de | dans le seq_record.id
		id_sequence = seq_record.id.replace('|','_') 
			
		file_name = id_sequence + ".fa"

		file_path = os.path.join(temp_directory, file_name)
		current_file = open(file_path, "w")
		current_file.write(seq_record.format("fasta"))


#############################
#  DATABASE INDEXATION      #
#############################


def make_BLAST_database(fasta_file):
 
    proc = subprocess.Popen([ "makeblastdb", "-in" , fasta_file, "-dbtype",
                                'nucl' ], stdout=subprocess.PIPE)
    sys.stderr.write(proc.stdout.read())
    for file_ext in ['.nhr', '.nin', '.nsq']:
        path = fasta_file + file_ext
        shutil.move(path, os.path.join('DBs', os.path.basename(path)))
    sys.stderr.write(("Getting %s and associated database files to the DBs "
                        "location\n") % (fasta_file))
    shutil.copy2(fasta_file, os.path.join('DBs', os.path.basename(fasta_file)))
    return os.path.basename(fasta_file).split('_')[0]


#######################
#        BLAST        #
#######################

def run_BLAST(temp_directory, database, output_dir):
  
 	for path, dirs, files in os.walk(temp_directory):
		for fasta_file in files:
			seq_number_search = re.search('^([0-9A-Za-z_]+).fa',fasta_file)
			seq_number = seq_number_search.group(1)
			fich_a_blaster = "tmp/"+ fasta_file
			output_file = output_dir + "blast_" + seq_number + ".xml"
			
			blastn_cline = NcbiblastnCommandline(query=fich_a_blaster, outfmt=5, evalue=0.001, db=database, out=output_file)  
			
			# two possibilities			
			#sys.stderr.write(str(blastn_cline)+"\n")
			stdout, stderr = blastn_cline()

			blastn_cline
			



###############################
# BLAST LAUNCHING  ...........

split_fasta(options.query,temp_directory)


if options.index == "yes":
	print "banque deja indexee"
	chemin_indexed_database = options.database

elif options.index == "no":
	print "banque sera indexee par le script"
	make_BLAST_database(options.database)
	chemin_indexed_database = indexed_database_directory + options.database


else:
	os.system('rm -rf output')
	os.system('rm -rf DBs')
	os.system('rm -rf tmp')
	sys.exit()


run_BLAST(temp_directory, chemin_indexed_database, output_directory)



###############################
# Traitement resultats.........

# e_value threshold
e_value_thresh = float('1e-100')

q = re.search('^(.*).fa[s]?[t]?[a]?$',options.query)
options.query_sans_extension = q.group(1)
ofh = open(options.query_sans_extension + "_best_hits.txt", "w")  

dico_res_blast = {}

# read blast results files
for root, dirs, files in os.walk("output", topdown=False):
    for name in files:
	blast_results_file = open(os.path.join(root, name))
	
	# each xml file is parsed
	blast_records = NCBIXML.parse(blast_results_file)
	for blast_record in blast_records:
		list_descript = []
		if blast_record.alignments:
			for alignment in blast_record.alignments:
				# HSP.SCORE if best blast score wanted
				#for hsp in alignment.hsps:
					#print hsp.score, hsp.expect, alignment.hit_def
				#print alignment.hit_def
				if str(blast_record.query) not in dico_res_blast:
					dico_res_blast[str(blast_record.query)]  = []
				list_descript.append(alignment.hit_def)
				dico_res_blast[str(blast_record.query)] = list_descript 

		else:
			ofh.write('pas alignment' + '\n')



for seq_id in dico_res_blast:

	search_new_seq_id = re.search('^(.*)\|ENSBTAT[0-9]+$',seq_id)

	if search_new_seq_id:
		new_seq_id = search_new_seq_id.group(1)
		ofh.write(new_seq_id + '\t')
	else:
		ofh.write(seq_id + '\t')


	acronym = ''
	for description in dico_res_blast[seq_id]:
		# search for 'description' containing species and acronym
		search_species_in_description = re.search(r".*\b(?=\w)%s\b(?!\w).*\(([A-Z0-9]+)\).*$" % options.species,description)
		search_acronym = re.search('.*\(([A-Z0-9]*)\).*$',description)
		if search_species_in_description:
			if search_acronym: # if species and acronym
				acronym = search_species_in_description.group(1)
				break
			else: # if species but no acronym
				break
			
		else:# no species in description
			continue





	########################
	# write in output file #
	########################

	if description:

		search_predicted = re.search(r"^PREDICTED: (\b(?=\w)%s\b(?!\w).*)$" % options.species,description)
		if search_predicted:
			new_description = search_predicted.group(1)
			description = new_description


		search_synthetic_construct = re.search('^Synthetic construct (.*)$',description)
		if search_synthetic_construct:
			new_description = search_synthetic_construct.group(1)
			description = new_description

		search_TPA = re.search('^TPA: (.*)$',description)
		if search_TPA:
			new_description = search_TPA.group(1)
			description = new_description

		search_TPA_inf = re.search('^TPA_inf: (.*)$',description)
		if search_TPA_inf:
			new_description = search_TPA_inf.group(1)
			description = new_description




		if acronym:
			ofh.write(description + '\t' + acronym + '\n')
		else:
			ofh.write(description + '\n')








	else: # si pas species on prend l acronym de la premiere description trouvee
		search_acronym2 = re.search('.*\(([A-Z0-9]*)\).*$',dico_res_blast[seq_id][0])
		if search_acronym2:
			acronym2 = search_acronym2.group(1)
			ofh.write(dico_res_blast[seq_id][0] + '\t' + acronym2 + '\n')
		else:
			ofh.write(dico_res_blast[seq_id][0] + '\n')




ofh.close()

# remove temporary directories
os.system("rm -rf tmp")
os.system("rm -rf output")
os.system("rm -rf DBs")

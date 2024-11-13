from Bio import SeqIO
import glob
import os


def add_new_results(protein: str) -> None:
	fasta = ""

	#Load the original fasta file content
	with open(f"./queries/{protein}.fasta") as file:
		fasta += file.read()
	#Add new line if needed
	if fasta[-1] != "\n":
		fasta += "\n"

	#Get the results
	for file in glob.glob(f"./results/{protein}_*.res"):
		specie = file.split("/")[-1][len(protein)+1:-4]
		hits = set()
		with open(file, "r") as f:
			for l in f:
				if ">>" in l:
					hits.add(l.split()[1])
		#Append the results fasta to the original fasta
		if len(hits) != 0:
			for seq_rec in SeqIO.parse(f"./metamonada_datasets/{specie}.fasta", "fasta"):
				if seq_rec.id in hits:
					fasta += f">{seq_rec.description}\n{seq_rec.seq.strip('*')}\n"

	#Write the new enriched fasta file
	with open(f"./queries/{protein}.fasta", "w") as f:
		f.write(fasta)

def build_new_hmm(protein: str) -> None:
	#Cluster sequences with 90% sequence identity
	os.system(f"mmseqs easy-cluster ./queries/{protein}.fasta ./queries/{protein}.clust ./queries/tmp --min-seq-id 0.9 -c 0.8")
	#Replace the original fasta file and clean the tmp files
	os.system("rm -R ./queries/tmp/")
	os.system(f"mv ./queries/{protein}.clust_rep_seq.fasta ./queries/{protein}.fasta")
	os.system(f"rm ./queries/{protein}.clust_*")

	#Mafft alignment && HMM profile build
	os.system(f"mafft --maxiterate 1000 --localpair ./queries/{protein}.fasta > ./queries/{protein}.ali")
	os.system(f"hmmbuild ./queries/{protein}.hmm ./queries/{protein}.ali")

def run_hmmer(protein: str) -> None:
	for file in glob.glob("./metamonada_datasets/*.fasta"):
		specie = file.split("/")[-1][:-6]
		os.system(f"hmmsearch -E 0.001 ./queries/{protein}.hmm {file} > ./results/{protein}_{specie}.res")
		print(f"hmmsearch done for {specie}.. ")

#Set number of iterative runs you want to perform
def iterative_hmmer(protein: str, iterations: int) -> None:
	for i in range(0,iterations):
		build_new_hmm(protein=protein)
		run_hmmer(protein=protein)
		add_new_results(protein=protein)
		print(f"run {i+1} finished.. ")

iterative_hmmer(protein=YOUR_PROTEIN, iterations=3)

from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from Bio import SeqIO
import glob
import os


def prepare_fasta():
	fasta = ""

	with open("vzarsky_dataset.tsv", "r") as f:
		for l in f:
			if len(l) > 1:
				line = l.split("\t")
				if line[7] == "Tim17" or line[7] == "Tim23" or line[7] == "Tim22":
					fasta += f">{line[5]} {line[6]} Tax={line[0]}\n"
					fasta += f"{line[8].strip()}\n"

	with open("vzarsky_dataset.fasta", "w") as f:
		f.write(fasta)

def download_models():
	for seq_rec in SeqIO.parse("vzarsky_dataset.fasta", "fasta"):
		seq_ID = seq_rec.id.split("|")[1]
		os.system(f"wget https://alphafold.ebi.ac.uk/files/AF-{seq_ID}-F1-model_v4.pdb")
		os.system(f"mv AF-{seq_ID}-F1-model_v4.pdb ./models/{seq_ID}.pdb")

def DeepAlign():
	for pdb_file in glob.glob("./models/*.pdb"):
		os.system(f"DeepAlign {pdb_file} 8SCX_Tim17.pdb -o {pdb_file}_ScTim17")
		os.system(f"mv {pdb_file}_ScTim17.fasta ./alignments/")
		os.system(f"rm {pdb_file}_ScTim17.*")

def read_alignments():

	results = ""

	for file in glob.glob("./alignments/*.fasta"):

		known = {}
		unknown = {}
		canal_indeces = []
		seq_ID = ""
		canal_fasta = ""

		for seq_rec in SeqIO.parse(file, "fasta"):
			if seq_rec.id.startswith("8SCX"):
				for residue_n, residue in enumerate(seq_rec.seq):
					known[residue_n] = residue
			else:
				seq_ID = seq_rec.id
				for residue_n, residue in enumerate(seq_rec.seq):
					unknown[residue_n] = residue

		sequence = ""
		current_residue_n = -1
		for seq_rec in SeqIO.parse("8SCX_canal.fasta", "fasta"):
			for i in known:
				if known[i] == "-":
					sequence += "-"
				else:
					current_residue_n += 1
					sequence += seq_rec.seq[current_residue_n]

		for position_n, position in enumerate(sequence):
			if position != "-" and not position.islower():
				canal_indeces.append(position_n)

		for i in canal_indeces:
			canal_fasta += unknown[i]

		results += f">{seq_ID}\n{canal_fasta}\n"

	with open("canals_residues.fasta", "w") as f:
		f.write(results)

def analyze_pi():

	charges_results = ""

	for seq_rec in SeqIO.parse("canals_residues.fasta", "fasta"):
		protein = PA(seq_rec.seq)
		charge = protein.charge_at_pH(7.2)
		charges_results += f"{seq_rec.id},{charge:.2f}\n"

	with open("canals_charges.csv", "w") as f:
		f.write(charges_results)



read_alignments()
analyze_pi()


# if __name__ == "__main__":
# 	prepare_fasta()
# 	download_models()
# 	DeepAlign()
# 	read_alignments()
# 	analyze_pi()
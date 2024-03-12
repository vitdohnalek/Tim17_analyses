import os
import glob
from Bio import SeqIO
from modeller import *
from modeller.automodel import *

def Align2D(model: str, template: str, chain: str="A"):
    env = Environ()
    aln = Alignment(env)
    mdl = Model(env, file=template, model_segment=('FIRST:A','LAST:A'))
    aln.append_model(mdl, align_codes=f'{template}{chain}', atom_files=f'{template}.pdb')
    aln.append(file=f'{model}.ali', align_codes=model)
    aln.align2d(max_gap_length=50)
    aln.write(file=f'{model}-{template}{chain}.ali', alignment_format='PIR')

def RunModel(model: str, template: str, chain: str="A"):
    env = Environ()
    a = AutoModel(env, alnfile=f'{model}-{template}{chain}.ali',
                  knowns=f'{template}{chain}', sequence=model,
                  assess_methods=(assess.DOPE,
                                  #soap_protein_od.Scorer(),
                                  assess.GA341))
    a.starting_model = 1
    a.ending_model = 1
    a.make()

def clean_tmp(model: str, template: str, chain: str="A"):
	os.system(f"mv {model}.*.pdb {model}.pdb")
	for file in glob.glob(f"{model}*"):
		if file != f"{model}.pdb":
			os.remove(file)

def convert_fasta_to_pir(fasta_file: str, seq_ID: str) -> str:

	pir_file = f">P1;{seq_ID}\nsequence:{seq_ID}:::::::0.00: 0.00\n"
	for seq_rec in SeqIO.parse(fasta_file, "fasta"):
		pir_file += str(seq_rec.seq)
	pir_file += "*"

	return pir_file

def multipred():
	for file in glob.glob("./fasta/*.fasta"):
		file_name = file.split("/")[-1][:-6]
		pir_file = convert_fasta_to_pir(fasta_file=file,seq_ID=file_name)

		with open(f"./{file_name}.ali", "w") as f:
			f.write(pir_file)

		Align2D(model=file_name, template="8SCX")
		RunModel(model=file_name, template="8SCX")
		clean_tmp(model=file_name, template="8SCX")

multipred()
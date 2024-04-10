import os
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.CheckSum import crc32
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
from Bio.SeqUtils.CodonUsageIndices import SharpEcoliIndex
from Bio.SeqUtils import six_frame_translations
from Bio.Seq import Seq
import gzip
from Bio import SeqIO
from math import floor
from argparse import ArgumentParser

argparse = ArgumentParser()
argparse.add_argument(
    "-i",
    "--infile",
    help="Path to input file with S. cerevisiae sequences",
    required=True,
)

argparse.add_argument(
    "-o",
    "--outfile",
    help="Path to output csv with S. cerevisiae sequences stats",
    required=True,
)
args = argparse.parse_args()

inf = args.infile
outf = args.outfile

def calculate_cai(dna, index=SharpEcoliIndex):
    cai = CodonAdaptationIndex()
    cai.set_cai_index(index)
    if len(dna) % 3 == 0:
        a = cai.cai_for_gene(dna)
    else:
        six_translated = six_frame_translations(dna)
        n = six_translated.split("\n")
        frames = {
            "0;F": n[5],
            "1;F": n[6],
            "2;F": n[7],
            "0;R": n[12],
            "1;R": n[11],
            "2;R": n[10],
        }
        ind = 0
        for i in list(frames.keys()):
            k = frames[i].replace(" ", "")
            if "M" in k and "*" in k:
                if i.split(";")[0] == "F" and k.index("M") < k.index("*"):
                    if len(k) <= len(dna) / 3:
                        ind = int(i.split("")[0])
                        break
                elif i.split(";")[0] == "R" and k.index("M") > k.index("*"):
                    if len(k) <= len(dna) / 3:
                        ind = len(dna) - int(i.split("")[0])
                        break
        if ind == 0:
            cods = 3 * floor(len(dna) / 3)
            dna = dna[:cods]
            a = cai.cai_for_gene(dna)
        elif 1 <= ind <= 2:
            if len(dna[ind:]) % 3 == 0:
                dna = dna[ind:]
            else:
                cods = 3 * floor((len(dna) - ind) / 3)
                dna = dna[ind : cods + ind]
                a = cai.cai_for_gene(dna)
        else:
            if len(dna[:ind]) % 3 == 0:
                dna = dna[ind:]
            else:
                cods = 3 * floor((len(dna) - ind) / 3)
                dna = dna[:cods]
                a = cai.cai_for_gene(dna)
    return a


def checksum(dna):
    return crc32(dna)


def hidrophobicity(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*", "")
    hydrophobicity_score = ProteinAnalysis(protein_sequence).gravy()
    return hydrophobicity_score


def isoelectric_pt(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*", "")
    isoelectric = ProteinAnalysis(protein_sequence).isoelectric_point()
    return isoelectric


def aromatic(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*", "")
    arom = ProteinAnalysis(protein_sequence).aromaticity()
    return arom


def instable(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*", "")
    inst = ProteinAnalysis(protein_sequence).instability_index()
    return inst


def weight(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*", "")
    wgt = ProteinAnalysis(protein_sequence).molecular_weight()
    return wgt


def sec_struct(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*", "")
    second_struct = ProteinAnalysis(protein_sequence).secondary_structure_fraction()
    return ",".join([str(s) for s in second_struct])


def mol_ext(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*", "")
    molar_ext = ProteinAnalysis(protein_sequence).molar_extinction_coefficient()
    return ",".join([str(s) for s in molar_ext])

def read_multiplefasta(multiple_fasta):
    genomes = []
    with open(multiple_fasta) as file:
        single_genomes = []
        counter = 0
        for line in file:
            if ">" not in line:
                seq = []
                for char in line:
                    if char != '\n':
                        seq.append(char)
                dnar = ""
                dnadef = dnar.join(seq)
                single_genomes.append(dnadef)
                counter += 1
                comp_gen = ""
                complete_genome = comp_gen.join(single_genomes)
            if ">" in line and counter == 0:
                continue
            if ">" in line and counter != 0:
                genomes.append(complete_genome) 
                single_genomes = []
        genomes.append(complete_genome)           
    
    return genomes


def check_if_in_header(classesdict: dict, header: str) -> str:
    for key in classesdict:
        if key in header:
            return classesdict[key]
    return ""



def eliminate_blankets(filepath: str):
    fp = open(filepath, "r+")
    lines = fp.readlines()
    newlines = [line for line in lines if line != "\n"]
    fp.seek(0)
    fp.truncate()
    for newline in newlines:
        fp.write(newline)
    fp.close()
    return filepath

from random import choice
def convert_degenerate_to_random(seq):
    """Converts degenerate nucleotide symbols to a random one of the ones they represent."""
    mapping = {'B': ['C','T'], 'D': ['A','T'], 'H': ['A','C'], 'K': ['G','T'], 'M': ['A','C'], 'N': ['A','C','G','T'], 'R': ['A','G'], 'S': ['C','G'], 'V': ['A','C','G'], 'W': ['A','T'], 'Y': ['C','T']}
    converted_seq = ""
    for char in seq:
        if char in mapping:
            converted_char = choice(mapping[char])
        else:
            converted_char = char
        converted_seq += converted_char
    return converted_seq

if __name__=="__main__":
    classes = {'ORF Uncharacterized': 'Uncharacterized_ORF', 'ORF Dubious': 'Dubious_ORF', 'ORF Verified': 'Verified_ORF'}
    print("Loading data...")
    fa = open(eliminate_blankets(inf), "r")
    lines = fa.readlines()
    fa.close()
    fasta = read_multiplefasta(inf)
    headers = [header for header in lines if header.startswith(">")]
    print("Done")
    print("Writing csv...")
    with open(outf, "w") as csv:
        csv.write(
            "ORF_TYPE,CAI,CHECKSUM,HIDROPHOBICITY,ISOELECTRIC,AROMATIC,INSTABLE,MW,HELIX,TURN,SHEET,MOL_EXT_RED,MOL_EXT_OX\n"
        )
        c = 0
        for i in range(len(fasta)):
            c += 1
            orf = check_if_in_header(classes, headers[i])
            if orf != "":
                try:
                    fasta[i] = convert_degenerate_to_random(fasta[i])
                    cai = calculate_cai(fasta[i])
                    cksm = checksum(fasta[i])
                    hydr = hidrophobicity(fasta[i])
                    isl = isoelectric_pt(fasta[i])
                    arm = aromatic(fasta[i])
                    inst = instable(fasta[i])
                    mw = weight(fasta[i])
                    se_st = sec_struct(fasta[i])
                    me = mol_ext(fasta[i])
                    csv.write(
                        f"{orf},{cai},{cksm},{hydr},{isl},{arm},{inst},{mw},{se_st},{me}\n"
                    )
                except Exception as e:
                    continue
            else:
                continue
            if c % 500 == 0:
                print(f"Processed {c} reads")
    csv.close()
    print("Done")


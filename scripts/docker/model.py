import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import VotingClassifier
from sklearn.ensemble import VotingClassifier, HistGradientBoostingClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.CheckSum import crc32
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
from Bio.SeqUtils.CodonUsageIndices import SharpEcoliIndex
from Bio.SeqUtils import six_frame_translations
from Bio.Seq import Seq
from Bio import SeqIO
import gzip
from math import floor
from sklearn.metrics import accuracy_score
from orfipy_core import orfs
import sys
import matplotlib.pyplot as plt


def load_data(infile):
    """Load data from infile if it is in fasta format (after having unzipped it, if it is zipped)"""
    if infile.endswith(".gz"):  # If file is gzipped, unzip it
        y = gzip.open(infile, "rt", encoding="latin-1")
        # Read file as fasta if it is fasta
        if (
            infile.endswith(".fasta.gz")
            or infile.endswith(".fna.gz")
            or infile.endswith(".fas.gz")
            or infile.endswith(".fa.gz")
        ):
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
        else:
            y.close()
            raise ValueError("File is the wrong format")
    # Read file directly as fasta if it is a not zipped fasta: handle also more uncommon extensions :-)
    elif (
        infile.endswith(".fasta")
        or infile.endswith(".fna")
        or infile.endswith(".fas")
        or infile.endswith(".fa")
    ):
        with open(infile, "r") as y:
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
    else:
        raise ValueError("File is the wrong format")


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


def longest_orf(coding):
    keys_M_starting = [
        key
        for key in list(coding.keys())
        if str(Seq(coding[key]).translate()).startswith("M")
    ]
    M_starting = [
        seq
        for seq in list(coding.values())
        if str(Seq(seq).translate()).startswith("M")
    ]
    lengths = [len(seq) for seq in M_starting]
    max_ind = lengths.index(max(lengths))
    return {keys_M_starting[max_ind]: M_starting[max_ind]}


def predict_orf(seq, minlen=45, maxlen=18000, longest_M_starting_orf_only=True):
    ls = orfs(seq, minlen=minlen, maxlen=maxlen)
    coding = {}
    count = 0
    for start, stop, strand, description in ls:
        count += 1
        coding.update({f"ORF.{count}": seq[int(start) : int(stop)]})
    if longest_M_starting_orf_only:
        print(
            "\n---------------------------\nWarning: option longest_M_starting_orf_only is set to True and thus you will get only the longest M-starting ORF; to get all the ORFs, set it to False\n---------------------------\n",
            file=sys.stderr,
        )
        return longest_orf(coding)
    return coding


def process_dna(fasta_file):
    fas = load_data(fasta_file)
    seqs = [seq for seq in list(fas.values())]
    heads = [seq for seq in list(fas.keys())]
    data = {}
    proteins = {}
    for i in range(len(seqs)):
        coding = predict_orf(seqs[i])
        open_reading_frames = list(coding.keys())
        for key in open_reading_frames:
            head = f"{heads[i]}.{key}"
            proteins.update({head: str(Seq(coding[key]).translate())})
            cai = calculate_cai(coding[key])
            cksm = checksum(coding[key])
            hydr = hidrophobicity(coding[key])
            isl = isoelectric_pt(coding[key])
            arm = aromatic(coding[key])
            inst = instable(coding[key])
            mw = weight(coding[key])
            se_st = sec_struct(coding[key]).split(",")
            se_st1 = se_st[0]
            se_st2 = se_st[1]
            se_st3 = se_st[2]
            me = mol_ext(coding[key]).split(",")
            me1 = me[0]
            me2 = me[1]
            n = pd.DataFrame(
                {
                    "CAI": [cai],
                    "CHECKSUM": [cksm],
                    "HIDROPHOBICITY": [hydr],
                    "ISOELECTRIC": [isl],
                    "AROMATIC": [arm],
                    "INSTABLE": [inst],
                    "MW": [mw],
                    "HELIX": [se_st1],
                    "TURN": [se_st2],
                    "SHEET": [se_st3],
                    "MOL_EXT_RED": [me1],
                    "MOL_EXT_OX": [me2],
                }
            )
            data.update({head: n})
    return data, proteins

if __name__ == "__main__":
    print("Loading data...")
    # Load the data from the CSV file
    data = pd.read_csv("../../data/scerevisiae.csv")
    print("Loaded data")

    print("Generating training and test data...")
    # Features
    X = data.iloc[:, 1:]

    # Labels
    y = data["ORF_TYPE"]


    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    print("Generated training and test data")

    print("Building and training the model...")
    # Create and train the Random Forest classifier
    clf4 = DecisionTreeClassifier()
    clf7 = HistGradientBoostingClassifier()
    clf8 = ExtraTreesClassifier()
    classifier = VotingClassifier([('dt', clf4), ('hgb', clf7), ('etc', clf8)], voting='hard')

    model = classifier.fit(X, y)  # Uncomment this line if clf needs training


    # Make predictions on the test set
    y_pred = model.predict(X)

    # Evaluate the accuracy of the model
    accuracy = accuracy_score(y, y_pred)
    print(f"Accuracy: {accuracy}")

    from joblib import dump

    print("Saving model...")
    dump(model, "SacCerML.joblib")
    print("Saved")

    print("All done")

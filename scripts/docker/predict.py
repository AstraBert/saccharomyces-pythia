from joblib import load
from model import process_dna



loaded_model = load("SacCerML.joblib")

def merge_fastas(fileslist):
    finale = []
    finalfile = fileslist[-1].split(".")[0]+"_mergedfastas.fasta"
    for fl in fileslist:
        f = open(fl, "r")
        lines = f.readlines()
        f.close()
        for line in lines:
            finale.append(line)
    fnlfl = open(finalfile, "w")
    for l in finale:
        if l.endswith("\n"):
            fnlfl.write(l)  
        else:
            fnlfl.write(l+"\n") 
    fnlfl.close()
    return finalfile

def predict_genes(infile, model=loaded_model):
    X, proteins = process_dna(infile)
    headers = list(X.keys())
    predictions = []
    for x in list(X.values()):
        p = model.predict(x)
        predictions.append(p)
    msg = []
    for i in range(len(predictions)):
        msg.append(
            f"{headers[i]} protein sequence is\n{proteins[headers[i]]}\nand is predicted as {predictions[i][0]}\n"
        )
    message = "".join(msg)
    return message


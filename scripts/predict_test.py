import pandas as pd
from joblib import load
from argparse import ArgumentParser
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report

argparse = ArgumentParser()
argparse.add_argument(
    "-i",
    "--infile",
    help="Path to input csv file with S. cerevisiae sequences stats",
    required=True,
)
args = argparse.parse_args()

inf = args.infile

model = load("SacCerML.joblib")

# Load the data from the CSV file
data = pd.read_csv(inf)

# Features
X = data.iloc[:, 1:]

# Labels
y = data["ORF_TYPE"] 

# Make predictions on the test set
y_pred = model.predict(X)

y_true = y

dictionnaire = {'transposable_element_gene': 0, 'pseudogene': 1, 'Uncharacterized_ORF': 2, 'Dubious_ORF': 3, 'Verified_ORF': 4, 'blocked_reading_frame': 5}
accuracy = accuracy_score(y_true, y_pred)
f1 = f1_score([dictionnaire[yt] for yt in y_true], [dictionnaire[yp] for yp in y_pred],  labels=list(set([dictionnaire[yt] for yt in y_true])), average='weighted')
precision = precision_score([dictionnaire[yt] for yt in y_true], [dictionnaire[yp] for yp in y_pred],  labels=list(set([dictionnaire[yt] for yt in y_true])), average='weighted')
recall = recall_score([dictionnaire[yt] for yt in y_true], [dictionnaire[yp] for yp in y_pred],  labels=list(set([dictionnaire[yt] for yt in y_true])), average='weighted')
classrep = classification_report(y_true, y_pred)


print(f"TEST on {inf}\n-Accuracy: {accuracy}\n-f1: {f1}\n-Precision: {precision}\n-Recall:{recall}\nClassification report:\n{classrep}\n\n")

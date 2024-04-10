import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import VotingClassifier
from sklearn.ensemble import VotingClassifier, HistGradientBoostingClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report
import matplotlib.pyplot as plt
import warnings

# Settings the warnings to be ignored 
warnings.filterwarnings('ignore') 



accs = []
recs = []
precs = []
f1s = []

dictionnaire = {'transposable_element_gene': 0, 'pseudogene': 1, 'Uncharacterized_ORF': 2, 'Dubious_ORF': 3, 'Verified_ORF': 4, 'blocked_reading_frame': 5}

import sys
kfolds = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
for kfold in kfolds:
    print("Loading data...", file=sys.stderr)
    # Load the data from the CSV file
    data = pd.read_csv("./data/scerevisiae.csv")
    print("Loaded data", file=sys.stderr)

    print("Generating training and test data...", file=sys.stderr)
    # Features
    X = data.iloc[:, 1:]

    # Labels
    y = data["ORF_TYPE"]


    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=kfold, random_state=42
    )
    print("Generated training and test data", file=sys.stderr)

    print("Building and training the model...", file=sys.stderr)
    # Create and train the Random Forest classifier
    clf4 = DecisionTreeClassifier()
    clf7 = HistGradientBoostingClassifier()
    clf8 = ExtraTreesClassifier()
    classifier = VotingClassifier([('dt', clf4), ('hgb', clf7), ('etc', clf8)], voting='hard')

    model = classifier.fit(X_train, y_train)  # Uncomment this line if clf needs training


    # Make predictions on the test set
    y_pred = model.predict(X_test)

    # Evaluate the the model
    y_true = y_test
    accuracy = accuracy_score(y_true, y_pred)
    f1 = f1_score([dictionnaire[yt] for yt in y_true], [dictionnaire[yp] for yp in y_pred],  labels=list(set([dictionnaire[yt] for yt in y_true])), average='weighted')
    precision = precision_score([dictionnaire[yt] for yt in y_true], [dictionnaire[yp] for yp in y_pred],  labels=list(set([dictionnaire[yt] for yt in y_true])), average='weighted')
    recall = recall_score([dictionnaire[yt] for yt in y_true], [dictionnaire[yp] for yp in y_pred],  labels=list(set([dictionnaire[yt] for yt in y_true])), average='weighted')
    classrep = classification_report(y_true, y_pred)
    print(f"Classification report for kfold {kfold}\n{classrep}")
    accs.append(accuracy)
    recs.append(recall)
    precs.append(precision)
    f1s.append(f1)

df = {"Kfold": kfolds, "Accuracy": accs, "Precision": precs, "Recall": recs, "F1_score": f1s}
datf = pd.DataFrame.from_dict(df)
datf.to_csv("scripts/kfoldval/kfold_crossval.csv", index=False)

for key in df:
    if key!="Kfold":
        fig, ax = plt.subplots(figsize=(10,5))
        ax.set_title(f"{key}")
        # Define label names
        labels = list(df["Kfold"])

        # Calculate height values for each bar
        heights = df[key]

        # Plot bars with specified width, color, and positional arguments
        bars = ax.bar(range(len(labels)), height=heights, align='center', color="blue")

        # Loop through all bars and annotate their heights
        for rect in bars:
            height = round(rect.get_height(), 2)  # Get rounded height value
            ax.annotate('{}'.format(height),  # Annotate the height value
                        xy=(rect.get_x() + rect.get_width() / 2, height),  # Position annotation at center of bar
                        xytext=(0, 3),  # Vertically shift annotation by 3 points
                        textcoords="offset points",
                        ha='center', va='bottom')
        # Add x axis tick labels
        ax.set_xticks([r for r in range(len(labels))])
        ax.set_xticklabels(labels)
        ax.set_ylabel("Crossvalidation value")
        ax.set_xlabel("Kfold")
        # ax.set_ylim(max([max(plotheigts), max(heights)])+50)
        fig.savefig(f"scripts/kfoldval/Kfold-{key}.png")
        plt.show()


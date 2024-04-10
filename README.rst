SacCerML: A Machine Learning Approach to Gene Prediction in Saccharomyces cerevisiae
================================================================================

Introduction
------------

SacCerML is a Python script that leverages machine learning and bioinformatics tools to predict genes in Saccharomyces cerevisiae (baker's yeast) genomic sequences. The script combines the `orfipy_core` library for ORF prediction and the `scikit-learn` library for Voting Classification.

Features
--------

The script extracts various features from the genomic sequences, including:

* Codon Adaptation Index (CAI)
* Checksum
* Hydrophobicity
* Isoelectric point
* Aromaticity
* Instability index
* Molecular weight
* Secondary structure
* Molar extinction coefficients

These features are then used to train a Decision Tree classifier to predict the ORF type of a given sequence.

Methods
--------

The script consists of four main methods:

1. **Data Loading and Preprocessing**: This method reads genomic sequences from a FASTA file, predicts ORFs, and extracts features. It also prepares a labeled dataset (`scerevisiae.csv`) for training and testing.
2. **Model Training**: This method splits the dataset into training and test sets and utilizes a Decision Tree classifier to build and train the gene prediction model.
3. **Prediction and Evaluation**: This method predicts gene types on the test set and calculates the accuracy of the model.

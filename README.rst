.. raw:: html

   <table>
   <tr>
   <td>
   <img src="https://img.shields.io/github/languages/top/AstraBert/saccharomyces-pythia" alt="GitHub top language">
   </td>
   <td>
   <img src="https://img.shields.io/badge/saccharomyces pythia-stable-green" alt="Static Badge">
   </td>
   <td>
   <img src="https://img.shields.io/badge/Release-v1.0.0-purple" alt="Static Badge">
   </td>
   <td>
   <img src="https://img.shields.io/badge/Docker_image_size-7GB-red" alt="Static Badge">
   </td>
   <td>
   <img src="https://img.shields.io/badge/Supported_platforms-linux/amd64-brown" alt="Static Badge">
   </td> 
   </tr>
   </table>

==============================================================================
saccharomyces-pythia: an ML/AI-integrated *Saccharomyces cerevisiae* assistant
==============================================================================

Introduction
============

**saccharomyces-pythia** is the new, rebranded v1.0.0 of SacCerML. 

Initially conceived as a Python script that leveraged machine learning and bioinformatics tools to predict genes in Saccharomyces cerevisiae (baker's yeast) genomic sequences, it is now a complete and AI-integrated tool that can help researchers both as a chatbot and as a ORF-predicter. 

SacCerML: the base ML model
===========================

Training
--------

Data and preprocessing
++++++++++++++++++++++

All the annotated coding DNA sequences for *S. cerevisiae* (strain S288C) were downloaded from Saccharomyces Genome Database.

These genetic sequences were split according to their ORF classification (verified, dubious, uncharacterized, pseudogene and transposable element) and for each of them the following parameters were calculated:

- Codon Adaptation Index
- Checksum 

After that, DNA was translated into aminoacids and other descriptors were retrieved:

- Hydrophobicity
- Isoelectric point
- Aromaticity
- Instability
- Molecular weight
- Secondary structure percentage (helix, turn and sheet)
- Molar extinction (both oxidized and reduced)


All the computed data were stored in a csv file, which was used to train a supervised ML model, a Voting Classifier (implemented in scikit-learn package), made up by HistGradient Boosting Classifier, a Decision Tree Clasifier and an Extra Tree Classifier.

Validation
++++++++++

The so-obtained machine-larning model (called SacCerML) was then evaluated on the entire training set, yielding a 99.93% accuracy. 

A key component of the training was k-fold crossvalidation. 

SacCerML was trained on increasingly wider percentages of the training data and tested on the remainder: it yielded a high accuracy (>84%) in all the tests, and the same goes for recall, f1 and precision score. You can find all the details in `the dedicated folder <https://github.com/AstraBert/tree/main/scripts/kfoldval>`_. 

From the classification reports it could be already seen a slight bias towards predicting verified and dubious ORFs, with more difficulty in predicting uncharacterized ORFs.

Testing
-------

Data and preprocessing
++++++++++++++++++++++

Data were collected from ORFs of 10 *Saccharomyces cerevisiae* strains, different from the one used for training:

- AWRI1631
- BC187
- BY4741
- CBS7960
- FL100
- g833-1B
- Kyokai7
- LalvinQA23
- Vin13
- YS9

A total of 54452 transcripts were collected and processed into csv file by extracting the previously mentioned features.

Validation
++++++++++

You can find all the test results `here <https://github.com/AstraBert/tree/main/test/test_results.stats>`_.

Generally, the model performed well: it had overall accuracy, f1, precision and recall score always above 86%. Nevertheless, the slight bias towards verified and dubious ORFs was confirmed, though uncharacterized ORFs too were well detected in several tests. 

saccharomyces-pythia: gene calling and AI integration
=====================================================


SacCerML has now reached a new stage of its development (v1.0.0), where it has been rebranded as **saccharomyces-pythia**.

You can now enjoy the following upgrades, that make it user-friendly and easy to install:

- `Gradio <https://www.gradio.app/>`_ chatbot interface running completely locally on your computer
- Gene calling with automated ORF detection thanks to `orfipy <https://pypi.org/project/orfipy/>`_: no need for preprocessing your reads, just upload one or more FASTA files with *S. cerevisiae* DNA sequences to the chatbot.
- AI assistant, built upon `EleutherAI/pythia-160-deduped-v0 <https://huggingface.co/EleutherAI/pythia-160m-deduped-v0>`_ finetuned on *Saccharomyces cerevisiae and its industrial applications* (Parapouli et al., 2020): this is a text-generation model that will reply to researcher questions (stil a beta feature, we become more stable in future releases).
- Docker image to download and run the application on your computer

References
----------

* Saccharomyces Genome Database: <https://www.yeastgenome.org/>
* Biopython: <https://biopython.org/>
* Scikit-learn: <https://scikit-learn.org/stable/>
* Gradio: <https://www.gradio.app/>
* orfipy: <https://pypi.org/project/orfipy/>
* EleutherAI/pythia-160-deduped-v0: <https://huggingface.co/EleutherAI/pythia-160m-deduped-v0>
* Parapouli et al., 2020: <https://doi.org/10.3934/microbiol.2020001>

Additionally, the following libraries and packages were used in the development of the machine learning model:

* NumPy: <https://numpy.org/>
* Pandas: <https://pandas.pydata.org/>

These libraries and packages were used for data manipulation, analysis, and model training.

License
-------

The project is hereby provided under MIT license.

If you are using saccharomyces-pythia for your work, please consider citing its author, `Astra Bertelli <https://astrabert.vercel.app>`_ 

*How was this README generated? Leveraging the power of AI with reAIdme, an HuggingChat assistant based on mistralai/Mixtral-8x7B-Instruct-v0.1. Go and give it a try at this link: <https://hf.co/chat/assistant/660d9a4f590a7924eed02a32!> 🤖*

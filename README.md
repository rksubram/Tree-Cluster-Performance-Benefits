# Downstream Performance Benefits of TreeCluster 

## Introduction:

Phylogenetic Clustering of data can serve as a powerful tool for studying medical infections and diseases. Data on the percentage of Disease-Causing microbes in every patient and a machine learning model can help classify a patient as potentially sick or healthy. In this report, we evaluate the performance of Tree-Cluster, an optimal clustering algorithm, based on its efficacy in downstream applications. 

## Problem Definition:

[TreeCluster](https://github.com/niemasd/TreeCluster) identifies the minimum number of leaf-clusters based on heterogeneity constraints and is useful for HIV transmission clustering. Here, we use stool samples of patients, with details of nucleotide percentages, to obtain Tree-Cluster Clustering result. This grouping is analyzed against machine learning models like Logistic Regression and Support Vector Machines with dependency on the threshold value for clustering also evaluated. Additionally, regressions are performed with clustering based on the maximum diameter of the cluster, sum of branch lengths and single linkage and with varying degrees of dimensionality reduction. The sparsity of the dataset implores us to attempt dimensionality reduction. We hypothesize that with either statistical/phylogenetic dimensionality reduction, we can filter the dataset to most relevant features with moderate accuracy predictions from our classification models. 

## Contents:
The code has been made in Jupyter Notebooks to provide an element of interactivity. You can open the notebook and run each snippet individually to see a demonstration of the workflow. The datasets have been uploaded as well, along with the detailed results dump files for the various methods followed.

## Team Members:
* Akash Boghani - aboghani@ucsd.edu
* Inderjot Singh - isaggu@ucsd.edu
* Raghav Subramanian - rksubram@ucsd.edu
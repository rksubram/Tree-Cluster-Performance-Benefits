# Tree-Cluster-Performance-Benefits
Exploring the Downstream Performance benefits of the TreeCluster 

Project Goals(Adding so we don't miss out anything) :- 

1. Develop algorithmic techniques for your problem
2. Improve scalability or accuracy of an exisiting technique
3. Perform scientific comparisons of scalabitlity and accuracy of existing methods 


Progress:- 
So far, added a shell script to be run from the main dir, that uses bash command line to filter only required field and extract a CSV with filetered metadata 

https://github.com/niemasd/TreeCluster

Things to do by Saturday
- Divide the data into training, val and test sets, choose the same seed to split for consistent comparison of results 
80:10:10 split with 43 as seed
- Run Linear Regression, SVM, Random forest models using all the features
- Compute correlation between features
- Combine features with high correlation
- Use PCA for feature dimensionality reduction, lets try different orders 5k, 1k, 100, 10 and check model performance

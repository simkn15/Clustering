Overview of files:
mainHclust.R
Main file for running the methods hclustDivisiveGap(), from file hcGap.R.

mainMeasure.R
Code used for getting the quality measure for all thresholds of a given dataset.
Runs tclust() in a range of thresholds with the method: doTclustWithinRange() from measureTclust.R

hcGap.R
Holds the method hclustDivisiveGap().
This is the final version of Transitivity Clustering extending into Hierarchical Clustering.
Includes randomization for each node in the dendrogram, with either random approach 3 or 4.
Returns the list: list(hc, gap)
hc: The object to draw the dendrogram of the clustering.
gap: All the levels for which a splits occured, including the costs for the actual data and the random data.
Method plotCostDifference() from utilities.R is used to plot the cost difference for all the splits: plotCostDifference(gap)

hcRandomizationRange.R
Used to find biggest cost gap at each node, running a range of thresholds.

fmeasure.R
HOlds all the code for the implementation of F1-measure

makeSimFileFromBlast.R
Code used for converting the BLAST output for the big dataset into a similarity file.

measureTclust.R
Holds method doTclustWithinRange() which uses tclust on the simMatrix within a specified range of thresholds.
Returns the overall best clustering.

randomization.R
Holds all methods for the randomization of a similarity matrix.

utilities.R
Holds various different methods, where many are used in the THC algorithm.
Includes several methods used for plotting.
buildSimilarityMatrix() is used to build similarity file for the small dataset.
buildSimilarityMatrixFromBlast() is used to build similarity file for the big dataset.

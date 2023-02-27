Data Format:
* Data files are stored in Data/[name of data set]/graph.txt
  * Small data sets (allsports, captchas, cora200, gym, landmarks) are provided
  * Large data sets (amazon, dblp, livejournal, orkut, youtube) are publicly available at snap.stanford.edu/data/#communities
* Delimiter is hard-coded in driver files
* For Correlation Clustering: 
  * The first line of the file must contain the total number of nodes (anything that follows it on the line will be ignored)
  * Rest of file lists positive edges as [node1] [node2]
* Bounds used for non-uniform cluster sizes are also provided in Data/[name of data set]/[name of data set]_bounds[example number].txt
  * Format: one bound per line, with each line corresponding to a node in the graph

To Compile: javac *.java
To Run: java [DriverName] [data set folder name]
* additional heap space may be needed for some experiments; increase the max heap size with the -Xmx flag

Drivers
-------

Hard coded parameters:
* delimiter for data set ("\\s" for all examples here)
* number of Pivot rounds, uniform size constraints, etc. 
* Input: positive edge list

RunImprovedAlg.java
* Method: Ji et al. LP rounding for uniform constraints
* Requires Google OR-Tools (developers.google.com/optimization)

RunMaxKCorrelation.java
* Method: Pivot, Vote, PLS, VLS for uniform constraints

RunMaxKTimed.java
* Method: Pivot, Vote, PLs, VLS for uniform constraints with 5-minute LS time limit

RunNonUniformAlg1.java
* Method: both Ji et al. LP rounding algorithms for non-uniform constraints
* Requires Google OR-Tools (developers.google.com/optimization)

RunNonUniformCorrelation.java
* Method: Pivot, Vote, PLS, VLS for non-uniform constraints

RunNonUniformOrder.java
* Method: Pivot, Vote, PLS, VLS for non-uniform constraints ordering by constraint size

RunNonUniformTimed.java
* Method: Pivot, Vote, PLS, VLS for non-uniform constraints with 5-minute LS time limit

RunNonUniformTimedOrder.java
* Method: Pivot, Vote, PLS, VLS for non-uniform constraints ordering by constraint size and with 5-minute LS time limit

RunPmAlg.java
* Method: Puleo and Milenkovic LP rounding for uniform constraints
* Requires Google OR-Tools (developers.google.com/optimization)


Code Files
----------

DNode.java
* Implementations of Vote

Helper.java
* Helper functions for reading data sets

Pair.java
* Custom data structure used for heap implementation

PKwik.java
* Pivot algorithm implementations

WriteNonUniformSizes.java
* Code to generate bounds for non-uniform cluster size constraints


Plots
-----

visualization folder contains Jupyter notebooks for plotting algorithm results

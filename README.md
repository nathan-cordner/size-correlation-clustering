Data Format:
* Data files are stored in Data/[name of data set]/graph.txt
  * All data sets are either publicly available or provided here
  * We provide data set samples to show the format for the various experiments 
* Delimiter is hard-coded in driver files
* For Correlation Clustering: 
  * The first line of the file must contain the total number of nodes (anything that follows it on the line will be ignored)
  * Rest of file lists positive edges as [node1] [node2]

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

RunMaxKCorrelation.java
* Method: Pivot, Vote, PLS, VLS for uniform constraints

RunMaxKTimed.java
* Method: Pivot, Vote, PLs, VLS for uniform constraints with 5-minute LS time limit

RunNonUniformAlg1.java
* Method: both Ji et al. LP rounding algorithms for non-uniform constraints

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

==================================================================================
README File of PFCLUST.jar	Author: Khadija Musayeva	Date:20.06.2013
==================================================================================
This README file describes how to run the PFClust program




I REQUIREMENTS
==================================================================================

	To run this program you need Java SE Runtime Environment 1.7
	For more information about Java environment refer to http://www.oracle.com/technetwork/java/javase/downloads/index.html





II CONFIGURATION
==================================================================================

	PFClust program reads an input file (similarity matrix) from the
	input directory and writes the results to the output directory. Details of the
	input directory where all similarity files should be and output directory
	where the results are written should be specified in the config.txt file. 
	Currently this file holds default values.


	NOTE: The config.txt should be in the same directory as PFClust.jar


1. CONTENTS OF config.txt
----------------------------------------------

	EXAMPLE:
	----------------------------------------------
	input_dir=data/input/
	output_dir=data/output/
	data=cath_11.txt,Density.txt,group_10.txt
	---------------------------------------------


a) input_dir 

	The input directory where all similarity files should be.
	The value of this parameter should be specified as it is shown in the above example.
	There should not be any space before and after the equals sign.
	The value of the input_dir should end with forward slash.
	It should be correctly specified relative or absolute path.
	Currently input_dir is in the same directory as PFClust.jar


b) output_dir

	The output directory where PFClust program writes the results.
	The value of this parameter should be specified as it is shown in the above example.
	There should not be any space before and after the equals sign.
	The value of the output_dir should end with forward slash.
	It should be correctly specified relative or absolute path.
	Currently output_dir is in the same directory as PFClust.jar


c) data

	It specifies the list of input files.
	It should contain at least one file.
	There should not be any space before and after the equals sign.
	There should not be any space before and after the comma sign.
	NOTE: These files should exist in the input_dir






III RUNNING INSTRUCTIONS
==================================================================================

	1. Open command line

	2. cd path_to/PFClust.jar 
	   change path_to to the real path where you have placed PFClust.jar and related files

	3. Type java -jar PFClust.jar (you might need to change file permissions)



	While the program is running it will print out the following:

	EXAMPLE:
	-----------------------------------------------------------------------------
	File: s_2.txt
	Clusters=2	Silhouette=0.656	Dunn=3.977	Threshold=0.94
	Execution time = 1.76 sec 
	-----------------------------------------------------------------------------

	File - the name of input file
	Result of the program -  Number of clusters	Silhouette	Dunn	Threshold
	Execution time - the elapsed time of clustering




IV OUTPUT
==================================================================================

	The program writes its output to output_dir.

	The output is written to two files:



	I DATA file
	==========================================================================
 	Related information about the clustering output is saved in  name_of_the_input_file.data 	

	Contents of the file:

 	EXAMPLE:
	------------------------------------------------------------------------
 	numclust	sil	dunn	threshold
	16		0.649	1.595	0.944

	numclust - Number of clusters
	sil - Silhouette
	dunn - Dunn
	threshold - Threshold
	------------------------------------------------------------------------




	II CLUSTER file
	==========================================================================
	1. Name of the file
	   The name of the file is a concatenation of  [name of the input file, number of clusters, silhoutte, dunn, threshold]     followed by .clu
           EXAMPLE: Density16064915950944.clu 


	2. Contents of the file
	   Element index	Cluster index


	   EXAMPLE:
	   ------------------------------------------------------------------------
	   Element index	Cluster index
		1	  		9
		2			9
		3			9
		4			9
		5			9
		6			9
		7			9
		8			9
		9			9
		10			9
		11			9
		12			10
	   ------------------------------------------------------------------------
===================================================================================




/**
 * Copyright (C) 2013  Khadija Musayeva, Lazaros Mavridis
 * University of St Andrews
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package main;
import java.util.ArrayList;

import utilities.Matrix;
import utilities.Util;
import clustering.ClusterConstruction;
import clustering.Clustering;
import clustering.ClusteringCollection;
import clustering.Threshold;


/**
 * This class implements the main loop and the convergence steps of the algorithm. 
 */
public class PFClust {

	/**
	 * Represents the number of iterations of the main loop
	 */
	public static int MAIN_LOOP = 4;


	/**
	 * Represents the number of iterations of the randomization loop
	 */
	public static int RANDOMIZATION_LOOP = 1000;


	/**
	 * The cut-off percentage of the selected threshold
	 */
	public static double P = 0.85;



	/**
	 * Array of files read from the input directory
	 * specified in the configuration file.
	 */
	public static String[] FILES;



	/**
	 * Represents the threshold for the Rand index
	 */
	public static double RAND_THRESHOLD = 0.99;


	/**
	 * Represents the maximum number of iterations
	 * to execute the convergence step of the algorithm.
	 */
	public static int maxIteration = 10;



	/**
	 * For each data set specified in the configuration file
	 * executes the PFClust clustering
	 */
	public static void run() {

		for(String fileName : FILES) {

			try{
				Util.print("File: " + fileName);
				long start = System.nanoTime();

				//Starts clustering
				Clustering bestClustering = process(Util.INPUT_DIR + fileName);

				long end = System.nanoTime();

				Util.printClusteringResults(bestClustering, end - start);
				
				//Writes the results of clustering to the output
				Util.writeResults(fileName, bestClustering);
				
			}

			catch(Exception e) {
				Util.print(e.getMessage());
			}

		}

	}



	/**
	 * Executes the main iteration and convergence steps of the algorithm
	 * @param String fileName
	 * @return a Clustering
	 */
	private static Clustering process(String fileName) throws Exception {

		Matrix matrix = new Matrix(fileName);

		ClusteringCollection clusteringCollection = new ClusteringCollection();
		Clustering bestClustering = null;


		for(int i = 0; i < MAIN_LOOP; i++) {

			/**
			 * Computes thresholds
			 */
			ArrayList<Double> thresholds = Threshold.estimate(matrix);


			/**
			 * Constructs clustering
			 */

			Clustering clustering  = ClusterConstruction.perform(matrix, thresholds);

			/**
			 * Adds the clustering into collection
			 */
			clusteringCollection.add(clustering);
			
		}


		double randIndex = clusteringCollection.getRandIndex();

		if(randIndex < RAND_THRESHOLD) {
			converge(matrix, clusteringCollection);
		}

		clusteringCollection.sort();
		bestClustering = clusteringCollection.getBestClustering();

		return bestClustering;
	}


	/**
	 * If the Rand index of the clustering collection is less 
	 * than a specified value, it builds a new clustering 
	 * and replaces the worst clustering in the list
	 * with the new one, if the new clustering is better 
	 * than the worst clustering in the list. 
	 * The process continues unless the Rand index of the collection
	 * is greater than 0.99
	 * @param final Matrix matrix
	 * @param ClusteringCollection clusteringCollection
	 */
	public static void converge(final Matrix matrix, ClusteringCollection clusteringCollection) {
		double randIndex = clusteringCollection.getRandIndex();
		int counter = 0;

		while(randIndex < RAND_THRESHOLD && counter <= maxIteration) {
			clusteringCollection.sort();
			Clustering worstClustering = clusteringCollection.getWorstClustering();
			ArrayList<Double> thresholds = Threshold.estimate(matrix);
			Clustering clustering  = ClusterConstruction.perform(matrix, thresholds);
			
			if(clustering.compareTo(worstClustering) > 0) {
				clusteringCollection.remove(worstClustering);
				clusteringCollection.add(clustering);
				}

			randIndex = clusteringCollection.getRandIndex();
			counter++;
		}


	}



}

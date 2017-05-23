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

package clustering;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;

import main.PFClust;
import utilities.Matrix;
import utilities.Util;


/**
 * Implements the threshold estimation step of the algorithm.
 * Contains methods to construct the underlying distribution and select
 * appropriate threshold values from this distribution.
 */

public class Threshold {

	/**
	 * Instantiates the Random class, passes it as an argument to the
	 * getDistribution method. Selects thresholds from the distribution
	 * and returns it as an output.
	 * @param matrix
	 * @return thresholds
	 */
	public static ArrayList<Double> estimate(Matrix matrix) {
		Random random = new Random();
		ArrayList<Double> distribution = getDistribution(random, matrix);
		ArrayList<Double> thresholds = selectThresholds(distribution);
		return thresholds;
	}


	/**
	 * Constructs a distribution from the similarities of randomly generated clusters
	 * @param random
	 * @param matrix
	 * @return a distribution for threshold selection
	 */
	private static ArrayList<Double> getDistribution(Random random, Matrix matrix) {
		int size = matrix.size();
		ArrayList<Double> distribution = new ArrayList<Double>();

		for(int k = 0; k < PFClust.RANDOMIZATION_LOOP; ++k) {

			int numberOfClusters = random.nextInt(size) + 1;

			int numberOfUnclusteredElements = size;

			HashSet<Integer> clusteredElements = new HashSet<Integer>();

			for( int i = 0; i < numberOfClusters; ++i ) {

				int numberOfElementsToCluster = 0;

				if(i == numberOfClusters - 1) {
					numberOfElementsToCluster = numberOfUnclusteredElements;
				}

				else {
					int upper_bound = numberOfUnclusteredElements - (numberOfClusters - i);
					numberOfElementsToCluster = random.nextInt(upper_bound+1) + 1;
					numberOfUnclusteredElements -= numberOfElementsToCluster;
				}

				Cluster cluster = new Cluster(matrix);

				for (int j = 0 ; j < numberOfElementsToCluster ; ++j) {
					int element = random.nextInt(size);

					while(clusteredElements.contains(element)) {
						element = random.nextInt(size);
					}

					cluster.addElement(element);
					clusteredElements.add(element);
				} 

				if(cluster.size() >= 3) {
					cluster.calculateSimilarity();
					double meanSimilarity = cluster.getSimilarity();
					distribution.add(meanSimilarity);
				}

			}

		}

		return distribution;

	}


	/**
	 * Selects thresholds from the provided distribution
	 * @param distribution
	 * @return thresholds
	 */
	public static ArrayList<Double> selectThresholds(ArrayList<Double> distribution) {
		//Sorts the list in descending order.
		Collections.sort(distribution, Collections.reverseOrder());

		int size = distribution.size();

		//Holds threshold values.
		ArrayList<Double> thresholds = new ArrayList<Double>();


		//Selects two values from 95% and 95.75% levels.
		thresholds.add(Double.valueOf(Util.df.format(distribution.get((int)(size * 0.05)))));
		thresholds.add(Double.valueOf(Util.df.format(distribution.get((int)(size * 0.025)))));

		//Selects eight values between 99%(including) and 100 - ((10/size) *100)% levels. 
		int start = (int) (size * 0.01);

		double values = (start - 10)/7;

		for(int i = 0; i < 8; ++i) {
			int index = (int)(start - i * values);
			thresholds.add(Double.valueOf(Util.df.format(distribution.get(index))));
		}

		//Selects the top ten values.
		for(int i = 9; i >= 0; --i) {
			thresholds.add(Double.valueOf(Util.df.format(distribution.get(i))));
		}

		return thresholds;
		
	}



}

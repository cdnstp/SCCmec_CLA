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
import utilities.Matrix;



/**
 * The purpose of this class is to build a separate thread for the clustering computation.
 */
public class ClusteringThread extends Thread {

	/**
	 * Represents the threshold according to which the clustering is formed.
	 */
	private double threshold;


	/**
	 * Contains Clustering objects formed by each thread.
	 */
	ClusteringCollection collection;


	/**
	 * Represents the similarity matrix.
	 */
	Matrix matrix;


	/**
	 * Instantiates a ClusteringThread object
	 * @param threshold
	 * @param collection
	 * @param matrix
	 */
	public ClusteringThread(double threshold, ClusteringCollection collection, Matrix matrix) {
		this.threshold = threshold;
		this.collection = collection;
		this.matrix = matrix;
	}


	/**
	 * Starts a new thread and invokes the clustering procedure.
	 * Adds the obtained clustering to the collection.
	 */
	public void run() {
		collection.add(ClusterConstruction.buildClustering(matrix, threshold));
	}


	

}















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

import utilities.Matrix;
import utilities.Util;


/**
 * Represents a clustering defined by a list of clusters. 
 * Contains methods to calculate the Silhouette coefficient, 
 * Dunn index, Davis Bouldin index, Rand index of the clustering. 
 * Provides the functionality to refine the clustering.
 * Enables the comparison of clustering objects based on their
 * quality by implementing the Comparable interface.
 */

public class Clustering implements Comparable<Clustering>{

	/**
	 * Represents a list of clusters
	 */
	private ArrayList<Cluster> clusters = new ArrayList<Cluster>(); 


	/**
	 * Represents the Silhoutte coefficient
	 */
	private double silhouetteWidth = 0;


	/**
	 * Represents the Dunn index
	 */
	private double dunnIndex = 0;


	/**
	 * Represents the Davies Bouldin index
	 */
	private double DB = 0;


	/**
	 * Constitutes a threshold upon which the clustering is formed
	 */
	private double threshold = 0;


	/**
	 * Represents the similarity matrix
	 */
	private Matrix matrix;


	/**
	 * Number of elements of the clustering 
	 * It is equal to the number of rows or columns 
	 * of the similarity matrix (this is a symmetric matrix)
	 */
	private int numberOfElements;


	/**
	 * A singleton cluster is a cluster that contains only one element.
	 * Represents the number of singletons in the whole clustering 
	 */
	private int numberOfSingletons;


	public Clustering() {}


	/**
	 * Constructs the clustering object with the specified similarity matrix
	 * @param matrix
	 */
	public Clustering(final Matrix matrix) {

		this.matrix = matrix;
		numberOfElements = matrix.size();

	}


	/**
	 * Adds a cluster to the list of clusters
	 * @param cluster
	 */
	public void addCluster(final Cluster cluster) {
		clusters.add(cluster);
	}



	/**
	 * Returns the cluster at the specified position.
	 * @param index
	 * @return the cluster at the specified index
	 */
	public Cluster getCluster(final int index) {
		return clusters.get(index);

	}




	/**
	 * @return the list of clusters
	 */
	public ArrayList<Cluster> getClusters() {
		return clusters;
	}




	/**
	 * Updates the list of clusters.
	 * @param clusters
	 */
	public void setClusters(final ArrayList<Cluster> clusters) {

		this.clusters = clusters;

	}



	/**
	 * The number of clusters in the list.
	 * @return the number of cluster
	 */
	public int size() {

		return clusters.size();

	}


	/**
	 * Removes all clusters from the list
	 */
	public void clear() {
		clusters.clear();
	}


	/**
	 * @return the number of singleton clusters
	 */
	public int getNumberOfSingletons() {
		return numberOfSingletons;
	}


	/**
	 * Updates the number of singleton clusters
	 * @param numberOfSingletons
	 */
	public void setNumberOfSingletons(int numberOfSingletons) {
		this.numberOfSingletons = numberOfSingletons;
	}


	/**
	 * @return the number of elements in the clustering
	 */

	public int getNumberOfElements() {
		return numberOfElements;
	}


	/**
	 * Updates the number of elements in the clustering
	 * @param numberOfElements
	 */
	public void setNumberOfElements(int numberOfElements) {
		this.numberOfElements = numberOfElements;
	}	



	/**
	 * Refines the clustering.
	 * Updates each element's membership based on the similarity between an element 
	 * and its cluster. The degree of similarity of an element to a cluster
	 * is decided based the given threshold.
	 * @param threshold
	 */
	public void refine(final double threshold) {
		for(int clusterIndex = 0; clusterIndex < size(); ++clusterIndex) {
			movePoints(clusterIndex,threshold); 
		}
	}


	/**
	 * Calculates the average of similarity, alpha, between an element e and its cluster c.
	 * If alpha is less than the specified threshold, then computes the average value 
	 * of similarity, beta, between the element and other clusters d.
	 * If alpha<beta then element e is assigned to d.
	 * @param clusterIndex
	 * @param threshold
	 */
	private void movePoints(int clusterIndex, double threshold) {

		Cluster currentCluster = getCluster(clusterIndex);

		boolean isClusterChanged = false;

		for(int memberIndex = 0; memberIndex < currentCluster.size(); ++memberIndex) {

			int member = currentCluster.getElementAt(memberIndex);

			double similarityToItsCluster = currentCluster.calculateElementToClusterAverageSimilarity(member); 

			if(similarityToItsCluster < threshold) {

				double bestSimilarity = 0;

				Cluster closestCluster = null;

				for(int otherClusterIndex = 0; otherClusterIndex < size(); ++otherClusterIndex) {

					if(clusterIndex == otherClusterIndex) {continue;}

					Cluster otherCluster = getCluster(otherClusterIndex);

					double similarity = ((otherCluster.size() - 1 + 0d) / otherCluster.size())
							* otherCluster.calculateElementToClusterAverageSimilarity(member);

					if(bestSimilarity < similarity) {

						bestSimilarity = similarity;
						closestCluster = otherCluster;

					}

				}


				/**
				 * If a given element is more similar to other cluster than to its own,
				 * then moved the element to that cluster. 
				 */
				if(similarityToItsCluster < bestSimilarity) {

					currentCluster.removeMember(member);
					closestCluster.addElement(member);
					isClusterChanged = true;
					memberIndex -= 1;

				}

			}

		}


		if(isClusterChanged) {
			movePoints(clusterIndex, threshold);
		}


	}


	/**
	 * @return the Silhouette coefficient of a clustering
	 */
	public double getSilhouette() {
		return Double.valueOf(Util.df.format(silhouetteWidth));
	}


	/**
	 * @return the Dunn index of a clustering
	 */
	public double getDunnIndex() {
		return Double.valueOf(Util.df.format(dunnIndex));
	}


	/**
	 * @return the Davies Bouldin index of a clustering
	 */
	public double getDBIndex() {
		return Double.valueOf(Util.df.format(DB));
	}


	/**
	 * @return the threshold upon which the clustering is formed
	 */
	public double getThreshold() {
		return Double.valueOf(Util.df.format(threshold));
	}


	/**
	 * Updates the threshold of the clustering
	 * @param threshold
	 */
	public void setThreshold(double threshold) {
		this.threshold = threshold;
	}



	/**
	 * Constructs an array that maps an element to its cluster.
	 * Index of the array constitutes the element number, whereas the value
	 * at the specific position indicates the cluster number.
	 * @return
	 */
	public int[] getElementToClusterMap() {

		int[] map = new int[numberOfElements];

		for(int i = 0; i < size(); ++i) {

			Cluster cluster = getCluster(i);

			for(int j = 0; j < cluster.size(); ++j) {
				int elm = cluster.getElementAt(j);	
				map[elm] = i;
			}
		}

		return map;

	}


	/**
	 * Calculates the Davies-Bouldin index.
	 */
	public void calculateDB() {

		if(size()==1) {
			DB = 100;
			return;
		}

		for(int i = 0; i < size() - 1; ++i) {
			double max_Ri = -1; 
			double s_i = getCluster(i).getClusterRadius();
			int v_i = getCluster(i).getCentroid(); 

			for(int j = i+1; j < size(); ++j) {
				int v_j = getCluster(j).getCentroid(); 
				double s_j = getCluster(j).getClusterRadius();
				double d_ij = 1 - matrix.getValueAt(v_i, v_j);

				double R_ij = (s_i + s_j)/(d_ij);

				if(max_Ri < R_ij) {
					max_Ri = R_ij;
				}
			}

			DB += max_Ri;

		}

		DB = DB/size();

	}


	/**
	 * Calculates the Dunn index.
	 */
	public void calculateDunnIndex() {

		double maxClusterRadius = 0;
		double minInterClusterDistance = 1;

		for(int i = 0; i < size() - 1; ++i) {

			//Finds the cluster with the maximum radius.
			double clusterSize = getCluster(i).getClusterRadius();

			if(maxClusterRadius < clusterSize) {
				maxClusterRadius = clusterSize;
			}


			//Finds the minimum inter-cluster distances.
			int centroidI = getCluster(i).getCentroid(); 

			for(int j = i+1; j < size(); ++j) {

				int centroidJ = getCluster(j).getCentroid(); 
				double interClusterDistance = 1 - matrix.getValueAt(centroidI, centroidJ);

				if(minInterClusterDistance > interClusterDistance) {
					minInterClusterDistance = interClusterDistance;
				}

			}

		}

		//The ratio of min inter-cluster distance to cluster radius is the Dunn index.
		//Checks the denominator to avoid the pathological case.
		dunnIndex = (maxClusterRadius == 0)? 0: minInterClusterDistance/maxClusterRadius;

	}


	/**
	 * Calculates the Silhouette coefficient of the whole clustering
	 */
	public void calculateSilhouette() {

		//The Silhouette width is zero if the clustering is a singleton (contains only one cluster).
		if(size() == 1) {
			silhouetteWidth = 0; 
			return;
		}


		double sumOfSilhouettes = 0;

		for(int i = 0; i < size(); ++i) {
			double silhouetteWidthCluster = calculateSilhouette(getCluster(i), i);
			sumOfSilhouettes += silhouetteWidthCluster;
		}
		
		silhouetteWidth = (sumOfSilhouettes)/matrix.size();
	
	}


	/**
	 * Calculates the silhouette of each element in the cluster,
	 * adds up all these values and returns the average.
	 * @param cluster
	 * @param clusterIndex
	 * @return the Silhouette coefficient of the cluster
	 */
	public double calculateSilhouette(final Cluster cluster, final int clusterIndex) {
		//Silhouette width is -1 if a given cluster consists of only one element
		if(cluster.size() == 1) {
			return -1;
		} 

		double sumOfSilhouettes = 0;

		for(int i = 0; i < cluster.size(); ++i) {

			int member = cluster.getElementAt(i);
			double ai = 0;
			double minbi = 100;
			double silhouetteWidth = 0;


			/**
			 * Calculates the average dissimilarity of a given element 
			 * to each element in the cluster.
			 */
			for(int k = 0; k < cluster.size(); ++k) {

				int anotherMember = cluster.getElementAt(k);
				if(member == anotherMember) {continue;}
				ai += 1 - matrix.getValueAt(member, anotherMember);

			}


			ai = (ai+0d)/(cluster.size() - 1);



			/**
			 * Calculates the average dissimilarity of a given element 
			 * to each element in other clusters.
			 */
			for(int otherClusterIndex = 0; otherClusterIndex < size(); ++otherClusterIndex) {

				double bi = 0;
				if(clusterIndex == otherClusterIndex) {continue;}
				Cluster otherCluster = getClusters().get(otherClusterIndex);

				for(int k = 0; k < otherCluster.size(); ++k) {

					int otherClusterMember = otherCluster.getElementAt(k);
					bi += 1 - matrix.getValueAt(member, otherClusterMember);

				}

				bi = (bi + 0d)/otherCluster.size();

				if(minbi > bi) { minbi = bi; }

			}


			silhouetteWidth = (minbi - ai) / Math.max(ai, minbi);
			sumOfSilhouettes += silhouetteWidth;

		}

		return sumOfSilhouettes;
	}





	/**
	 * Calculates the Rand index agreement of two clusterings.
	 * @param otherClustering
	 * @return the Rand index of pairs of clusterings
	 */
	public double getRandIndex(Clustering otherClustering) {

		int[][] confMatrix = new int[size()][otherClustering.size()];
		int[] clusters = getElementToClusterMap();
		int[] clusters2 = otherClustering.getElementToClusterMap();


		for(int i = 0; i < clusters.length; ++i) {
			confMatrix[clusters[i]][clusters2[i]] += 1;
		}


		int row = confMatrix.length;
		int column = confMatrix[0].length;


		int d = 0;
		for(int i = 0; i < row; ++i) {
			for(int j = 0; j < column; ++j) {
				int sum = 0; int value = confMatrix[i][j];
				for(int ii = i + 1; ii < row; ++ii) {
					for(int jj = j + 1; jj < column; ++jj) {
						sum += confMatrix[ii][jj];
					}
				}
				for(int ii = i + 1; ii < row; ++ii) {
					for(int jj = 0; jj < j; ++jj) {
						sum += confMatrix[ii][jj];
					}
				}
				d += value * sum;
			}
		}


		int a = 0, b = 0, c = 0, square_matrix_size = Math.min(row, column);

		for(int i = 0; i < row; ++i) {

			if(i < square_matrix_size) {
				int value = confMatrix[i][i];
				a += (value * (value - 1))/2;
			}

			for(int j = 0; j < column; ++j) {

				if(confMatrix[i][j] == 0) {continue;}

				for(int jj = (j + 1); jj < column; ++jj) {

					b += confMatrix[i][j] * confMatrix[i][jj];

				}

				break;

			}

		}


		for(int i = 0; i < column; ++i) {

			for(int j = 0; j < row; ++j) {

				if(confMatrix[j][i] == 0) {continue;}

				for(int jj = (j + 1); jj < row; ++jj) {

					if(confMatrix[jj][i] == 0) {continue;}

					c += confMatrix[j][i] * confMatrix[jj][i];

				}

				break;

			}
		}


		return (a +d + 0d)/(a + b + c + d);
	}




	/**
	 * Compares two clusterings based on their Silhouette coefficient.
	 * Uses the Dunn index, if the difference of silhouettes of given clusterings 
	 * is within a certain threshold epsilon = 0.01.
	 */
	@Override
	public int compareTo(Clustering otherClustering) {
		if(otherClustering==null) {return 1;} 
		if(Math.abs(getSilhouette()- otherClustering.getSilhouette())<=0.01) {
			if(getDunnIndex() > otherClustering.getDunnIndex()) return 1;
			else if (getDunnIndex() < otherClustering.getDunnIndex()) return -1;
			else return 0;
			}
		else {
			if (getSilhouette() - otherClustering.getSilhouette() > 0.01) {return 1;}
			else return -1;
		}
	}



	/**
	 * Constructs the string representation of the clustering.
	 */
	@Override
	public String toString() {

		return  "Clusters = " + size() + 
				"\tSilhouette = " + Util.df.format(getSilhouette()) + 
				"\tDunn = " + Util.df.format(getDunnIndex()) + 
				"\tThreshold = " + Util.df.format(getThreshold())
				;	

	}



}
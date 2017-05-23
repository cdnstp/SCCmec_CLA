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


/**
 * Represents a cluster defined by a list of elements.
 * Provides methods to calculate cluster similarity, centroid and radius 
 * along with other elementary list operations. 
 */

public class Cluster {

	/**
	 * Contains cluster elements
	 */
	private ArrayList<Integer> members;

	
	/**
	 * Represents the similarity matrix
	 */
	private Matrix matrix;

	
	/**
	 * Represents cluster similarity
	 */
	private double similarity;


	/**
	 * Initializes the list of elements
	 */
	public Cluster() {
		members = new ArrayList<Integer>();
	}

	/**
	 * Initializes the list of elements
	 * Initializes the matrix object
	 * @param matrix
	 */
	public Cluster(final Matrix matrix) {
		this();	
		this.matrix = matrix;

	}


	/**
	 * Adds an element to the list of elements
	 * @param element
	 */
	public void addElement(final int element) {
		members.add(element);
	}


	/**
	 * Retrieves the element at the specified index
	 * @param index
	 * @return an element at the specified index.
	 */
	public int getElementAt(final int index) {
		return members.get(index);
	}



	/**
	 * Returns true if a cluster contains the specified element
	 * @param element
	 * @return true if a cluster contains the specified element
	 */
	public boolean isMember(final int element) {
		for(int i = 0; i < size(); ++i) {
			if(members.get(i) == element) {
				return true;
			}
		}
		return false;
	}



	/**
	 * Calculates cluster similarity 
	 */
	public void calculateSimilarity() {
		int n = size();

		if(n < 3) {
			similarity = -1; 
		}

		else { 
			double sumOfSimilarities = calculateSum();
			similarity = ((2 + 0d) / ((n-1) * n)) * sumOfSimilarities;
		}

	}


	/**
	 * The similarity of a given element to the cluster
	 * is calculated as the sum of similarities between the element
	 * and each member of the cluster. 
	 * @param element
	 * @return the similarity of a given element to the cluster.
	 */
	public double calculateElementToClusterSimilarity(int element) {
		double sumOfSimilarity = 0;

		for(int i = 0; i < size(); ++i) {
			int member = getElementAt(i);
			if(member == element) {continue;}
			sumOfSimilarity += matrix.getValueAt(element, member); 
		}
		return sumOfSimilarity;
	}



	/**
	 * Computes the average value of the sum of similarities between
	 * the given element and the cluster. 
	 * @param element
	 * @return the average value of element to cluster similarity.
	 */
	public double calculateElementToClusterAverageSimilarity(int element) {
		double sumOfSimilarity = calculateElementToClusterSimilarity(element);
		return sumOfSimilarity/(size() - 1);
	}


	/**
	 * Retrieves the similarity of the cluster
	 * @return
	 */
	public double getSimilarity() {
		return similarity;
	}


	/**
	 * Updates the similarity of the cluster
	 * @param similarity
	 */
	public void setSimilarity(final double similarity) {
		this.similarity = similarity;
	}



	/**
	 * @return a list of cluster members
	 */
	public ArrayList<Integer> getMembers() {
		return members;
	}


	/**
	 * @return the number of elements in the cluster
	 */
	public int size() {
		return members.size();
	}


	/**
	 * Updates elements of the cluster
	 * @param members
	 * @param threshold
	 */
	public void setMembers(ArrayList<Integer> members, double threshold) {
		this.members = members;
	}


	/**
	 * Removes the element from the list of elements at the specified position
	 * @param member
	 */
	public void removeMember(int member) {
		members.remove((Integer)member);
	}



	/**
	 * The radius of the cluster is calculated as the average of the distance between 
	 * the cluster centroid and its elements.
	 * @return averaged distance between the cluster centroid and cluster elements
	 */
	public double getClusterRadius() {
		int centroid = getCentroid();
		double radius = 0;

		for(int i = 0; i < size(); ++i) {

			int element = getElementAt(i);
			if(centroid == element) {continue;}
			radius += 1 - matrix.getValueAt(element, centroid);

		}

		radius = radius / (size() - 1);
		return radius;
	}



	/**
	 * The centroid of the cluster is the element with the highest similarity to its cluster.
	 * @return an element with the highest similarity to its cluster
	 */
	public int getCentroid() {
		double maxSimilarity = -1;
		int centroid = -1;

		for(int i = 0; i < size(); ++i) {
			int element = getElementAt(i);
			double similarity = calculateElementToClusterSimilarity(element);

			if(maxSimilarity < similarity) {
				maxSimilarity = similarity;
				centroid = element;
			}
		}

		return centroid;
	}


	/**
	 * Calculates the sum of pair-wise similarities between cluster members
	 * @return sum of pair-wise similarities between cluster members
	 */
	public double calculateSum() {
		double sumOfSimilarities = 0;

		for (int i = 1; i < size(); i++) {
			for(int k = 0; k < i ; k++) {
				sumOfSimilarities += matrix.getValueAt(members.get(i), members.get(k));
			}
		}
		return sumOfSimilarities;
	}


	/**
	 * Represents the string representation of the cluster
	 */
	public String toString() {
		String output = "";

		for(int j = 0; j < members.size(); ++j) {
			output += members.get(j);
			if(j != members.size()-1) {
				output += ",";
			}

		}

		return output;

	}




}
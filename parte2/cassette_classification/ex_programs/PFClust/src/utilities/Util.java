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


package utilities;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Properties;
import java.util.Scanner;

import main.PFClust;
import clustering.Clustering;


/**
 * Contains general purpose methods
 */
public class Util {

	/**
	 * Formats decimal numbers
	 */
	public static DecimalFormat df = new DecimalFormat("#.###");


	/**
	 * Represents directory that contains input data
	 */
	public static String INPUT_DIR;


	/**
	 * Output is written to this directory
	 */
	public static String OUTPUT_DIR;



	/**
	 * Configuration details are read from this directory
	 */
	public static String CONFIG_FILE = "config.txt";



	/**
	 * Writes the clustering data to the corresponding files
	 * @param source
	 * @param clustering
	 */
	public static void writeResults(String source, final Clustering clustering) {
		//Stores elements and their corresponding clusters.
		//Each file is named in the following format:
		//datasetName_numberOfClusters_Silhoutte_Dunn_Threshold.clu
		String append = clustering.toString().replaceAll("[\\.\t]", "");
		String clusterFileName = OUTPUT_DIR + source.substring(0, source.lastIndexOf(".")) + append + ".clu"; 

		//Contains the following information about the cluster:
		//number of clusters, silhouette width, dunn index, threshold and the corresponding cluster file.
		String dataFileName = OUTPUT_DIR + source.substring(0, source.lastIndexOf(".")) + ".data"; 

		try {
			//Initializes the corresponding file writers.
			File dataFile = new File(dataFileName);
			File clusterFile = new File(clusterFileName);
			FileWriter fileWriteData = null;
			FileWriter fileWriteCluster = null;

			//Writes clustering data to the specified file
			if(!dataFile.exists()) {
				fileWriteData = new FileWriter(dataFile, true);
				fileWriteData.write("numclust\tsil\tdunn\tthreshold\r\n");
			}
			else{
				fileWriteData = new FileWriter(dataFile, true);
			}
			fileWriteData.write(clustering + "\r\n");


			//Writes clusters to the specified file
			if(!clusterFile.exists()) {
				fileWriteCluster = new FileWriter(clusterFile, true);
				int[] clusters = clustering.getElementToClusterMap();
				for(int i = 0; i < clusters.length; ++i) {
					fileWriteCluster.write((i+1) + "\t" + (clusters[i]+1) + "\r\n");
				}
			}


			fileWriteData.close();
			fileWriteCluster.close();

		}

		catch(Exception e) {}

	}



	/**
	 * Creates a dissimilarity matrix based on the given similarity matrix.
	 * @param source
	 * @param destination
	 */
	public static void writeMatrix(final String source, final String destination) {

		try {
			int i = 0;
			FileWriter fw = new FileWriter(destination, true); 

			Scanner input = new Scanner(new File(source));

			while(input.hasNextLine()) {

				Scanner colReader = new Scanner(input.nextLine());

				while(colReader.hasNextDouble()) {

					fw.write(df.format(1 - colReader.nextDouble()) + "\t");//appends the string to the file

				}

				i += 1;
				fw.write("\r\n");
				colReader.close();
			}

			print(i);
			fw.close();
			input.close();

		}

		catch(Exception e) {
			print(e.getMessage());
		}

	}



	/**
	 * Reads the configuration file
	 */
	public static void loadConfiguration() {

		Properties prop = new Properties();

		try {
			prop.load(new FileInputStream(CONFIG_FILE));
			INPUT_DIR = prop.getProperty("input_dir");
			OUTPUT_DIR = prop.getProperty("output_dir");
			PFClust.FILES = prop.getProperty("data").split(",");
		} 

		catch (IOException e) {
			print(e.getMessage());
		}

	}


	/**
	 * Prints out the object to the output
	 * @param object
	 */
	public static void print(final Object object) {
		System.out.println(object);
	}

	/**
	 * Prints out clustering results and execution time
	 * @param clustering
	 * @param time
	 */
	public static void printClusteringResults(Clustering clustering, long time) {
		System.out.println(clustering);
		System.out.println("Execution time = " + Util.df.format((time * 1e-9)) + " sec.");
		System.out.println("---------------------------------------------------------------------");
	}
	
	
	
	/**
	 * Copyright information
	 */
	public static void printCopyright() {

		StringBuilder str = new StringBuilder();
		str.append("-----------------------------------------------------------------");
		str.append("\r\nCopyright (C) 2013  Khadija Musayeva, Lazaros Mavridis\r\n");
		str.append("University of St Andrews\r\n\r\n");
		str.append("This program is free software: you can redistribute it and/or modify\r\n");
		str.append("it under the terms of the GNU General Public License as published by\r\n");
		str.append("the Free Software Foundation, either version 3 of the License, or\r\n");
		str.append("(at your option) any later version.\r\n");
		str.append("This program is distributed in the hope that it will be useful,\r\n");
		str.append("but WITHOUT ANY WARRANTY; without even the implied warranty of\r\n");
		str.append("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\r\n");
		str.append("GNU General Public License for more details.\r\n");
		str.append("You should have received a copy of the GNU General Public License\r\n");
		str.append("along with this program. If not, see <http://www.gnu.org/licenses/>\r\n");
		str.append("-----------------------------------------------------------------\r\n");
		System.out.println(str.toString());

	}


}

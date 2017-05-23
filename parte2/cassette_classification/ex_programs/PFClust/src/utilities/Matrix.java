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

import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;


/**
 * Represents the similarity matrix backed by a two-dimensional array. 
 */

public class Matrix {

	/**
	 * Represents the two dimensional array that holds 
	 * pairwise similarities of elements in the given data-set.
	 */
	private double[][] data;


	/**
	 * Instantiates a Matrix object by reading the provided similarity matrix 
	 * into the two dimensional array.
	 * @param fileName
	 */
	public Matrix(final String fileName) throws Exception {
		this.data = readMatrix(fileName);
	}


	/**
	 * Retrieves the value at the specified row and column position.
	 * @param row
	 * @param column
	 * @return array
	 */
	public double getValueAt(final int row, final int column) {
		return data[row][column];
	}


	/**
	 * @return the number of rows in the matrix
	 */
	public int getRowSize() {
		return data.length;
	}


	/**
	 * @return the number of columns in the matrix
	 */
	public int getColumnSize() {
		return data[0].length;
	}


	/**
	 * @return the number of rows (it is a symmetric matrix, so the number of rows is equal to number of columns).
	 */
	public int size() {
		return data.length;
	}


	/**
	 * Uses the FileChannel in order to efficiently read the input from the specified 
	 * file and constructs the corresponding two dimensional array
	 * @param fileName
	 * @return two dimensional array of pairwise similarities between elements in the provided data set.
	 */
	private double[][] readMatrix(final String fileName) throws Exception {
		FileInputStream fis = null;
		double[][] matrix = new double[2][];

			fis = new FileInputStream(fileName);
			FileChannel fileChannel = fis.getChannel();
			ByteBuffer byteBuffer = ByteBuffer.allocate(131072);
			int bytes = fileChannel.read(byteBuffer);
			StringBuilder str = new StringBuilder();
			int row = 0;
			while(bytes!=-1) {
				byteBuffer.flip();
				ArrayList<Double> temp = new ArrayList<Double>();
				while (byteBuffer.hasRemaining()) {
					char ch = (char)byteBuffer.get();
					str.append(ch);

					if(ch == '\n') {
						int column = 0;
						String rowString[] = str.toString().split("\\s+");

						for(String s : rowString) {
							if(row == 0) {
								temp.add(Double.valueOf(s));
							}

							else {
								matrix[row][column] = Double.valueOf(s);
								column += 1;
							}

						}

						if(row == 0) {
							matrix = new double[temp.size()][temp.size()];
							for(int m = 0; m < matrix.length; ++m) {
								matrix[0][m] = temp.get(m);
							}
						}

						str = new StringBuilder();
						row += 1;

					}
				}
				
				byteBuffer.clear();
				bytes = fileChannel.read(byteBuffer);

			}

			int column = 0;
		
			String rowString[] = str.toString().split("\\s+");

			if(rowString.length == matrix.length) {
				for(String s : rowString) {
					matrix[row][column] = Double.valueOf(s);
					column += 1;
				}
			}

			if(fis!=null) {
				fis.close();
			}



		return matrix;
		
	}




}

import java.util.*;

/*
Name:  Magdalena Karski
MacID: karskim
Student Number:  001436728
Description:  This program multiplies n matrices together, finds the inverse of the resultant matrix, 
and then prints the elements of the inverted matrix row by row on a single line.
*/


public class HWK2_karskim {
	
	// Copies over the elements of an array into a new array. This is used both by the mult_mat and inverse functions. 
	public static double[][] deep_copy(double[][] original) {
		// new matrix will also be a 2D matrix. The first dimension of this matrix is the length of the original matrix.
		double[][] copied = new double[original.length][];
		// for every row 
	      for (int i=0; i < original.length; i++) {
	    	// Arrays.copyOf can deep copy a 1D array. Each row of the original matrix is copied into the corresponding row of the new matrix.  
	        copied[i] = Arrays.copyOf(original[i], original[i].length);
	      }
	      // The copied matrix is returned.
	      return copied;
	}
	// multiplies n matrices together
	public static double[][] mult_mat(int num_mat, int[][] matdim, double[][][] matval) {
		// If there's only one matrix
		if (num_mat == 1) {
		// return that matrix
			return matval[0];
		}
		// If there is more than one matrix, the first matrix is deep-copied to mat. 
		double[][] mat = deep_copy(matval[0]);
		
		// for every matrix in matval
		for (int i = 1; i < matval.length; i++) {
			// the number of rows of the first matrix is rowA
			int rowA = mat.length;
			// the number of columns of the first matrix is colA
			int colA = mat[0].length;
			// the number of rows of the second matrix is rowB
			int rowB = matval[i].length;
			// the number of rows of the second matrix is colB
			int colB = matval[i][0].length;
			// if the number of columns of the first matrix and the number of rows of the second matrix are not equal
			if(colA != rowB) {
				// The matrices cannot be multiplied together
				throw new IllegalArgumentException("Multiplication error");
			}
			// create a new 2D array where the product will be stored
			double[][] newmat = new double[rowA][colB];
			// for all rows in the first matrix
			for (int j = 0; j < rowA; j++) {
				// for all columns in the second matrix
				for (int k = 0; k < colB; k++) {
					// for every element in the current row of the first matrix (current column of the second)
					for (int l = 0; l < colA; l++) {
						// dot product of current row of A and current column of B
						// multiply the element in the current row of the first matrix by the corresponding element in the current column of the second matrix
						// add these together to get the value of the element in the resultant matrix
						newmat[j][k] += mat[j][l] * matval[i][l][k];
							
					}
				}
			}
			// create a deep copy of new_mat
			mat = deep_copy(newmat);
		}
		// once all of the matrices have been multipled together return the result
		return mat;
	}
	
	// create an nxn identity matrix. This is used by the inverse function. 
	public static double[][] identity(int n) {
		// create a new 2D array, nxn
        double[][] I = new double[n][n];
        // for the length of the array
        for (int i = 0; i < n; i++)
        // replace the pivot position with a 1. There are already zeroes everywhere else. 
            I[i][i] = 1;
        // return the nxn identity matrix
        return I;
    }					
// Computes the inverse of nxn matrix using Gauss-Jordan elimination
	public static double[][] inverse(double[][] mat) {
		// the number of rows is assigned to row
		int row = mat.length;
		// the number of columns is assigned to column
		int col = mat[0].length;
		// check if nxn
		if (row != col) {
			// otherwise the inverse function has been given a non-invertible function
			throw new IllegalArgumentException("Matrix not invertible");
		}
		// create an identity matrix the same dimenion as the matrix you want to invert
		double[][] I = identity(row);
		for (int j = 0; j < col; j++) {
			// if there is no pivot in the jth row of the jth column
			if (mat[j][j] == 0) {
				// search for a pivot
				int found = 0;
				for (int k = 0; k < row; k++) {
				// if you find a non-zero element in another row, the rows will be swapped, and this will be the new pivot
					if (mat[k][j] != 0) {
						// switch the found variable from 0 to 1 to indicate that you have found a possible pivot
						found = 1;
						// swap rows
						// make a deep copy of the pivot row
						double[] CopyRowA = Arrays.copyOf(mat[j], col);
						// assign the row where you found a new pivot to the original pivot row
						mat[j] = mat[k];
						// assign the original pivot row to the row where you found the new pivot
						mat[k] = CopyRowA;
						// repeat the same procedure for matrix I (originally the identity matrix)
						// copy original pivot row
						double[] CopyRowAI = Arrays.copyOf(I[j], col);
						// replace old pivot row with new pivot row
						I[j] = I[k];
						// assign the copy of the old pivot row to the position of the new pivot row
						I[k] = CopyRowAI;
					}
				}
				// if you have a column of zeros, the matrix is non-invertible
				if (found == 0) {
				// Let the user know that the matrix provided to the inverse function is not invertible.
					throw new IllegalArgumentException("Matrix not invertible");
				}
			}
			// reduce the matrix. original matrix will have numbers != to zero on the diagonal and zeros everywhere else
			for (int i = 0; i < row; i++) {
				// a is the value of the pivot
				double a = mat[j][j];
				// b is the value of the element of the row you want to reduce in the current column
				double b = mat[i][j];
				// this counts zeros produced by row reduction operations. If count equals the length of the row, you have a row of zeroes and the matrix is non-invertible.
				int count = 0;
				for (int m = 0; m < col; m++) {
					// skip the row containing the pivot
					if (i != j) {
						// replace the element in the original matrix in the ith row jth column with the previous element multiplied by the value at the pivot
						// and subtract from that the corresponding element in the row containing the pivot multiplied by b
						mat[i][m] = a*mat[i][m] - b*mat[j][m];
						// repeat the same procedure described above for the identity matrix
						I[i][m] = a*I[i][m] - b*I[j][m];
					}
					// if the operations specified above leave you with an element with value zero in the original matrix
					if (mat[i][m] == 0) {
					// increment count by 1	
						count += 1;
					}
				}
				// if you get a row of zeroes (if the number of zeros is equal to the length of the row)
				if (count == col) {
				// let the user know that the matrix given to the inverse function is a non-invertible matrix.
					throw new IllegalArgumentException("Matrix not invertible");
				}
			}
		}
		// To get the inverted matrix I, last step is to divide the elements of a row in I by the value of the element at the pivot position in that row in the original matrix 
		// for every row of matrix being inverted
		for (int i = 0; i < row; i++) {
		// get the value of mat at the pivot position in the current row
			double f = mat[i][i];
			// for the length of the row
			for (int j = 0; j < col; j++) {
				// divide out f from every element in the current row
				I[i][j] = I[i][j]/f;
			}
		}
		// return the inverted matrix
		return I;
	}
 
	// This is the main body. It extracts the number of matrices, the dimensions of the matrices, and the elements of each matrix from args.
	// First the function mult_mat is called. This multiplies n matrices. The result of this function is given to the function inverse, which computes the inverse of the matrix.
	// In the case that only one matrix is given, inverse computes the inverse of that original matrix. 
	// The elements of the inverted matrix are formatted to two decimal places and printed one one line. 
	
	public static void main(String[] args) {
		// Let the user know if no arguments were provided
		if (args == null) { throw new IllegalArgumentException("No arguments"); }
		// The first element given to args is the number of matrices
		int num_mat = Integer.parseInt(args[0]);
		// args[0] must be at least 1 
		if (num_mat == 0) { throw new IllegalArgumentException("Must have at least one matrix"); };
		
		// A 2D array containing matrix dimensions is initialized.
		// To access the number of columns of the second matrix: mat_dim[1][1]
		int[][] matdim = new int[num_mat][2];
		// c is the index of args. Matrix dimensions are extracted from args and placed into a 2D array. 
		int c = 1;
		// for the total number of matrices 
		for (int k = 0; k < matdim.length; k++) {
			// fill in the dimensions of the row first and then the column
			for (int j = 0; j < matdim[0].length; j++) {
				// an element from args is parsed to an integer and placed in matdim
				matdim[k][j] = Integer.parseInt(args[c]);
				// the index of args increments by 1
				c+=1;
			}
		}
		// Matrix values are extracted and placed into a 3D array. The each element of the first level of this array contains a matrix. 
		// The second level has the rows of each matrix. The third level has the columns. 
		// The element in the 3rd row and second column of the first matrix would be accessed this way: mat_val[0][2][1]
		// c is also used here to identify the index of args
		double[][][] matval = new double[num_mat][][];
			// elements are filled for all matrices
			for (int i = 0; i < num_mat; i++) {
				// matdim stores the dimensions of the matrix. The index of the first level of the array corresponds to the same matrix as the one that the index of the 
				// first level of matval corresponds to.
				// row stores the number of rows of the current matrix
				int row = matdim[i][0];
				// col stores the number of columns of the current matrix
				int col = matdim[i][1];
				// a new matrix is created at the ith position in the first level of matval with the dimensions of the of the current matrix
				matval[i] = new double[row][col];
				// Now elements are extracted from args and placed into this matrix. 
				// the matrix is filled one row at a time, starting from the first one
				for (int j = 0; j < matdim[i][0]; j++) {
					// elements are added to each row starting from the left-most element
					for (int k = 0; k < matdim[i][1]; k++) {
						// the value at args at c is parsed into a double and placed into the jth row and kth column of the ith matrix in matval
						matval[i][j][k] = Double.parseDouble(args[c]);
						// the value of c, used to keep track of the current index of args, is incremented by 1
						c+=1;
					}
				}
			}
		// This calls the multiply matrix function and then the inverse matrix function. 
		//	If there is only one matrix, mult_mat returns that matrix.
		double[][] Mult = mult_mat(num_mat, matdim, matval);
		// inverse takes the result of mult and inverts it. If there is only one matrix, Mult does not transform that matrix and it is inverted as it is. 
		double[][] Inv = inverse(Mult);
		// Each element of the inverted matrix is printed out one element at a time from the first to the last row
		System.out.println();
		for(int i=0; i < Inv.length; i++)
		    // for every element in the current row
			for(int j = 0; j < Inv[0].length; j++){
			   // Print all the elements on a single line.
				System.out.print(Inv[i][j] + " ");
			}


	}	
	
}


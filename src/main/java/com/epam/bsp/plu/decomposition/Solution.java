package com.epam.bsp.plu.decomposition;

public class Solution {
	
	public final static double EPS = 1e-8;
	
	/**
	* Returns a PLU decomposition of a given matrix.
	*
	* NOTE: It is guaranteed that a given matrix is square.
	* NOTE: It is guaranteed that the PLU decomposition exists.
	* 
	* NOTE: To avoid issues with precision, use EPS = 1e-8 to detect zeros:
	* |x| < EPS => We assume that x is zero.
	* 
	* @param matrix a given square matrix for which we want to get the PLU decomposition.
	* @return the resulting PLU decomposition consists of the following matrices:
	* - P - permutation matrix
	* - L - lower triangular matrix
	* - U - upper triangular matrix
	*/
	public static PluView pluDecomposition(double[][] matrix) {
		//put your code here
		return new PluView();
	}
	
	/**
	* Returns the inverse for a given matrix.
	*
	* NOTE: Use the `pluDecomposition` method from the previous coding exercise.
	* NOTE: It is guaranteed that a given matrix is square.
	* NOTE: It is guaranteed that the PLU decomposition exists.
	* 
	* NOTE: To avoid issues with precision, use EPS = 1e-8 to detect zeros:
	* |x| < EPS => We assume that x is zero.
	* 
	* @param matrix a given square matrix for which we want to get the inverse matrix.
	* @return the inverse matrix.
	*/
	public static double[][] getInverseMatrix(double[][] matrix) {
		//put your code here
		return new double[0][];
	}

}

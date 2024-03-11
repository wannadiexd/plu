package com.epam.bsp.plu.decomposition;

public class Solution {
  
  public final static double EPS = 1e-8;
  
  /**
   * Returns a PLU decomposition of a given matrix.
   *
   * @param matrix a given square matrix for which we want to get the PLU decomposition.
   * @return the resulting PLU decomposition consists of the following matrices:
   * - P - permutation matrix
   * - L - lower triangular matrix
   * - U - upper triangular matrix
   */
  public static PluView pluDecomposition(double[][] matrix) {
    int n = matrix.length;
    int[] permutation = new int[n];
    double[][] L = new double[n][n];
    double[][] U = new double[n][n];

    for (int i = 0; i < n; i++) {
      permutation[i] = i;
    }

    for (int i = 0; i < n; i++) {
      L[i][i] = 1.0;
      System.arraycopy(matrix[i], 0, U[i], 0, n);
    }

    for (int i = 0; i < n; i++) {
      int pivotRow = i;
      for (int j = i + 1; j < n; j++) {
        if (Math.abs(U[j][i]) > Math.abs(U[pivotRow][i])) {
          pivotRow = j;
        }
      }

      double[] tempRow = U[i];
      U[i] = U[pivotRow];
      U[pivotRow] = tempRow;

      int temp = permutation[i];
      permutation[i] = permutation[pivotRow];
      permutation[pivotRow] = temp;

      for (int j = i + 1; j < n; j++) {
        double factor = U[j][i] / U[i][i];
        L[j][i] = factor;
        for (int k = i; k < n; k++) {
          U[j][k] -= factor * U[i][k];
        }
      }
    }

    return new PluView(permutation, L, U);
  }
  
  /**
   * Returns the inverse for a given matrix.
   *
   * @param matrix a given square matrix for which we want to get the inverse matrix.
   * @return the inverse matrix.
   */
  public static double[][] getInverseMatrix(double[][] matrix) {
    PluView plu = pluDecomposition(matrix);
    double[][] L = plu.getL();
    double[][] U = plu.getU();
    int[] permutation = plu.getPermutation();
    int n = matrix.length;
    double[][] inverse = new double[n][n];

    for (int i = 0; i < n; i++) {
      double[] b = new double[n];
      b[i] = 1.0;

      double[] Pb = new double[n];
      for (int j = 0; j < n; j++) {
        Pb[j] = b[permutation[j]];
      }

      double[] y = new double[n];
      for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++) {
          sum += L[j][k] * y[k];
        }
        y[j] = (Pb[j] - sum) / L[j][j];
      }

      double[] x = new double[n];
      for (int j = n - 1; j >= 0; j--) {
        double sum = 0.0;
        for (int k = j + 1; k < n; k++) {
          sum += U[j][k] * x[k];
        }
        x[j] = (y[j] - sum) / U[j][j];
      }

      for (int j = 0; j < n; j++) {
        inverse[j][i] = x[j];
      }
    }

    return inverse;
  }

}

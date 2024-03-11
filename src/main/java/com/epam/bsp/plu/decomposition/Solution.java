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
    int n = matrix.length;
    double[][] perm = new double[n][n];
    double[][] lower = new double[n][n];
    double[][] upper = new double[n][n];

    // Initialize perm matrix as identity matrix
    for (int i = 0; i < n; i++) {
      perm[i][i] = 1;
    }

    // Perform PLU decomposition
    for (int i = 0; i < n; i++) {
      // Pivot row
      int pivot = i;
      for (int j = i + 1; j < n; j++) {
        if (Math.abs(matrix[j][i]) > Math.abs(matrix[pivot][i])) {
          pivot = j;
        }
      }
      // Swap rows in perm matrix
      double[] temp = perm[i];
      perm[i] = perm[pivot];
      perm[pivot] = temp;
      // Swap rows in matrix
      double[] tempRow = matrix[i];
      matrix[i] = matrix[pivot];
      matrix[pivot] = tempRow;
      
      for (int j = i; j < n; j++) {
        upper[i][j] = matrix[i][j];
        for (int k = 0; k < i; k++) {
          upper[i][j] -= lower[i][k] * upper[k][j];
        }
      }
      
      for (int j = i + 1; j < n; j++) {
        lower[j][i] = matrix[j][i];
        for (int k = 0; k < i; k++) {
          lower[j][i] -= lower[j][k] * upper[k][i];
        }
        lower[j][i] /= upper[i][i];
      }
    }

    return new PluView(perm, lower, upper);
  }
  
  /**
   * Returns the inverse for a given matrix.
   *
   * NOTE: Use the pluDecomposition method from the previous coding exercise.
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
    PluView plu = pluDecomposition(matrix);
    double[][] perm = plu.getPerm();
    double[][] lower = plu.getLower();
    double[][] upper = plu.getUpper();
    int n = matrix.length;
    double[][] inverse = new double[n][n];

    // Solve LY = P
    double[][] y = new double[n][n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double sum = perm[i][j];
        for (int k = 0; k < i; k++) {
          sum -= lower[i][k] * y[k][j];
        }
        y[i][j] = sum;
      }
    }

    // Solve UX = Y
    double[][] x = new double[n][n];
    for (int i = n - 1; i >= 0; i--) {
      for (int j = 0; j < n; j++) {
        double sum = y[i][j];
        for (int k = i + 1; k < n; k++) {
          sum -= upper[i][k] * x[k][j];
        }
        x[i][j] = sum / upper[i][i];
      }
    }

    return x;
  }

}

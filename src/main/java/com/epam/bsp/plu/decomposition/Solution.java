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
    double[][] matrixCopy = new double[n][n];

    // Copy the input matrix to avoid modifying the original matrix
    for (int i = 0; i < n; i++) {
      System.arraycopy(matrix[i], 0, matrixCopy[i], 0, matrix[i].length);
    }

    // Initialize perm matrix as identity matrix
    for (int i = 0; i < n; i++) {
      perm[i][i] = 1;
    }

    // Perform PLU decomposition
    for (int i = 0; i < n; i++) {
      // Pivot row
      int pivot = i;
      for (int j = i + 1; j < n; j++) {
        if (Math.abs(matrixCopy[j][i]) > Math.abs(matrixCopy[pivot][i])) {
          pivot = j;
        }
      }
      // Swap rows in perm matrix if a pivot is found
      if (pivot != i) {
        double[] temp = perm[i];
        perm[i] = perm[pivot];
        perm[pivot] = temp;
        // Swap rows in matrixCopy
        double[] tempRow = matrixCopy[i];
        matrixCopy[i] = matrixCopy[pivot];
        matrixCopy[pivot] = tempRow;
      }

      // Fill the upper matrix
      for (int j = i; j < n; j++) {
        upper[i][j] = matrixCopy[i][j];
        for (int k = 0; k < i; k++) {
          upper[i][j] -= lower[i][k] * upper[k][j];
        }
      }

      // Fill the lower matrix
      for (int j = i + 1; j < n; j++) {
        if (Math.abs(upper[i][i]) < EPS) { // Check for division by zero
          continue; // Skip the division if the pivot element is too close to zero
        }
        lower[j][i] = matrixCopy[j][i] / upper[i][i];
        for (int k = 0; k < i; k++) {
          lower[j][i] -= lower[j][k] * upper[k][i] / upper[i][i];
        }
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
        if (Math.abs(upper[i][i]) < EPS) { // Check for division by zero
          throw new ArithmeticException("Division by zero encountered at row " + i);
        }
        x[i][j] = sum / upper[i][i];
      }
    }

    // Copy the result to the inverse matrix
    for (int i = 0; i < n; i++) {
      System.arraycopy(x[i], 0, inverse[i], 0, x[i].length);
    }

    return inverse;
  }

}

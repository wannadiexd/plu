package com.epam.bsp.plu.decomposition;

public class Solution {

    public final static double EPS = 1e-8;

    // Helper class to store the PLU decomposition
    public static class PluView {
        public double[][] P;
        public double[][] L;
        public double[][] U;

        public PluView(double[][] P, double[][] L, double[][] U) {
            this.P = P;
            this.L = L;
            this.U = U;
        }
    }

    // Method to perform PLU decomposition
    public static PluView pluDecomposition(double[][] matrix) {
        int n = matrix.length;
        double[][] P = new double[n][n];
        double[][] L = new double[n][n];
        double[][] U = new double[n][n];

        // Initialize permutation matrix to identity
        for (int i = 0; i < n; i++) {
            P[i][i] = 1;
        }

        // Copy matrix to U
        for (int i = 0; i < n; i++) {
            System.arraycopy(matrix[i], 0, U[i], 0, n);
        }

        // Decompose matrix into PLU
        for (int i = 0; i < n; i++) {
            // Permute rows if needed
            for (int j = i; j < n; j++) {
                if (Math.abs(U[j][i]) > EPS) {
                    double[] temp = U[i];
                    U[i] = U[j];
                    U[j] = temp;

                    temp = P[i];
                    P[i] = P[j];
                    P[j] = temp;
                    break;
                }
            }

            // Compute L and U matrices
            for (int j = i + 1; j < n; j++) {
                L[j][i] = U[j][i] / U[i][i];
                for (int k = i; k < n; k++) {
                    U[j][k] -= L[j][i] * U[i][k];
                }
            }
        }

        // Fill diagonal of L with 1s
        for (int i = 0; i < n; i++) {
            L[i][i] = 1;
        }

        return new PluView(P, L, U);
    }

    // Method to calculate the inverse of a matrix
    public static double[][] getInverseMatrix(double[][] matrix) {
        PluView plu = pluDecomposition(matrix);
        int n = matrix.length;
        double[][] inverse = new double[n][n];

        // Solve for the inverse using the PLU decomposition
        for (int col = 0; col < n; col++) {
            double[] e = new double[n];
            e[col] = 1;
            // Forward substitution to solve for L * y = P * e
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = e[i];
                for (int j = 0; j < i; j++) {
                    y[i] -= plu.L[i][j] * y[j];
                }
                y[i] /= plu.L[i][i];
            }

            // Backward substitution to solve for U * x = y
            for (int i = n - 1; i >= 0; i--) {
                inverse[i][col] = y[i];
                for (int j = i + 1; j < n; j++) {
                    inverse[i][col] -= plu.U[i][j] * inverse[j][col];
                }
                inverse[i][col] /= plu.U[i][i];
            }
        }

        // Multiply by P transpose to get the final result
        double[][] result = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    result[i][j] += inverse[i][k] * plu.P[k][j];
                }
            }
        }

        return result;
    }
}

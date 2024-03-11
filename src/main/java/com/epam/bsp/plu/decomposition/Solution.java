package com.epam.bsp.plu.decomposition;

public class Solution {

    public final static double EPS = 1e-8;

    public static class PluView {
        double[][] P;
        double[][] L;
        double[][] U;

        public PluView(double[][] P, double[][] L, double[][] U) {
            this.P = P;
            this.L = L;
            this.U = U;
        }
    }

    public static PluView pluDecomposition(double[][] matrix) {
        int n = matrix.length;
        double[][] P = new double[n][n];
        double[][] L = new double[n][n];
        double[][] U = new double[n][n];

        for (int i = 0; i < n; i++) {
            P[i][i] = 1;
        }

        for (int i = 0; i < n; i++) {
            System.arraycopy(matrix[i], 0, U[i], 0, n);
        }

        for (int i = 0; i < n; i++) {
            // Find pivot
            int max = i;
            for (int j = i + 1; j < n; j++) {
                if (Math.abs(U[j][i]) > Math.abs(U[max][i])) {
                    max = j;
                }
            }

            double[] temp = U[i];
            U[i] = U[max];
            U[max] = temp;

            temp = P[i];
            P[i] = P[max];
            P[max] = temp;

            for (int j = i + 1; j < n; j++) {
                L[j][i] = U[j][i] / U[i][i];
                for (int k = i; k < n; k++) {
                    U[j][k] -= L[j][i] * U[i][k];
                }
            }
        }

        for (int i = 0; i < n; i++) {
            L[i][i] = 1;
        }

        return new PluView(P, L, U);
    }

    public static double[][] getInverseMatrix(double[][] matrix) {
        int n = matrix.length;
        double[][] inverse = new double[n][n];

        PluView plu = pluDecomposition(matrix);

        for (int col = 0; col < n; col++) {
            double[] e = new double[n];
            e[col] = 1;

            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = e[(int) plu.P[i][col]];
                for (int j = 0; j < i; j++) {
                    y[i] -= plu.L[i][j] * y[j];
                }
            }

            for (int i = n - 1; i >= 0; i--) {
                inverse[i][col] = y[i] / plu.U[i][i];
                for (int j = i + 1; j < n; j++) {
                    inverse[i][col] -= plu.U[i][j] * inverse[j][col] / plu.U[i][i];
                }
            }
        }

        return inverse;
    }
}

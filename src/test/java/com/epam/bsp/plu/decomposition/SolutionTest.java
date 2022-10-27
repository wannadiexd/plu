package com.epam.bsp.plu.decomposition;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import org.junit.jupiter.api.Test;

public class SolutionTest {
	
	@Test
	public void testPluDecompositionSample() {
		double[][] matrix = {
				{1, 2, 1}, 
				{1, 2, 2},
				{2, 1, 1}
		};
		
		PluView pluView = Solution.pluDecomposition(matrix);

		assertMatricesEqual(
			new double[][] {
				{1, 0, 0},
				{0, 0, 1},
				{0, 1, 0}
			}, pluView.getPerm()
		);
		
		assertMatricesEqual(
				new double[][] {
					{1, 0, 0},
					{2, 1, 0},
					{1, 0, 1}
				}, pluView.getLower()
			);

		assertMatricesEqual(
				new double[][] {
					{1, 2, 1},
					{0, -3, -1},
					{0, 0, 1}
				}, pluView.getUpper()
			);

	}

	@Test
	public void testGetInverseMatrixSample() {
		double[][] matrix = {
				{1, 2, 1}, 
				{1, 2, 2},
				{2, 1, 1}
		};
		
		double[][] matrixCopy = getMatrixCopy(matrix);
		
		double[][] matrixInv = Solution.getInverseMatrix(matrixCopy);
		
		double[][] identityCalc = getMultiplication(matrix, matrixInv);
		
		double[][] identity = {
				{1, 0, 0}, 
				{0, 1, 0},
				{0, 0, 1}
		};
		
		assertMatricesEqual(identity, identityCalc);
		
	}
	
	private static void assertMatricesEqual(double[][] expected, double[][] initial) {
		assertTrue(expected.length == initial.length);
		
		assertTrue(IntStream.range(0, expected.length)
				.allMatch(i -> expected[i].length == initial[i].length));
		
		assertArrayEquals(
				convertMatrixToArray(expected), 
				convertMatrixToArray(initial), 
				Solution.EPS
		);

	}
	
	private static double[] convertMatrixToArray(double[][] matrix) {
		return Arrays.stream(matrix)
				.flatMapToDouble(DoubleStream::of)
				.toArray();
	}

	private static double[][] getMatrixCopy(double[][] matrix) {
		return IntStream.range(0, matrix[0].length)
				.mapToObj(r -> IntStream.range(0, matrix.length)
						.mapToDouble(c -> matrix[c][r])
						.toArray())
				.toArray(double[][]::new);
	}

	private static double[][] getMultiplication(double[][] a, double[][] b) {
		return IntStream.range(0, a.length)
				.mapToObj(i -> IntStream.range(0, b[0].length)
						.mapToDouble(j -> IntStream.range(0, b.length)
								.mapToDouble(k -> a[i][k] * b[k][j]).sum())
						.toArray())
				.toArray(double[][]::new);
	}
}

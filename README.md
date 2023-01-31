# PLU decomposition

## Purpose

The following coding exercises are designed to test your knowledge of the following concepts:
* PLU decomposition
* The inverse of a square matrix

## Overview

The coding exercises cover the following practical problems:
* Getting the PLU decomposition of a given matrix
* Using the PLU decomposition to get the inverse of a given matrix

## Coding exercises

### Exercise 1: Find the PLU decomposition of a given matrix

Your task is to implement the following function to find the PLU decomposition of a given square matrix:

```java
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
	public static PluView pluDecomposition(double[][] matrix) 
```

**Additional materials (in case you get stuck):**

* https://www.youtube.com/watch?v=E3cCRcdFGmE&ab_channel=TheBrightSideofMathematics
* https://www.nagwa.com/en/explainers/976193728703/


**Example:**

```math
A=\begin{pmatrix}
1 & 2 & 1 \\
1 & 2 & 2 \\
2 & 1 & 1
\end{pmatrix}
```

Expected result:
```math
P = \begin{pmatrix}
1 & 0 & 0\\
0 & 0 & 1\\
0 & 1 & 0
\end{pmatrix}

L = \begin{pmatrix}
1 & 0 & 0\\
2 & 1 & 0\\
1 & 0 & 1
\end{pmatrix}

U = \begin{pmatrix}
1 & 2 & 1\\
0 & -3 & -1\\
0 & 0 & 1
\end{pmatrix}
```

### Exercise 2: Find the inverse of a given matrix using the PLU decomposition 

Your task is to implement the following function to invert matrices using the PLU decomposition:

```java
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
	public static double[][] getInverseMatrix(double[][] matrix) 
```

**Additional materials (in case you get stuck):**

* http://home.cc.umanitoba.ca/~farhadi/Math2120/Inverse%20Using%20LU%20decomposition.pdf
* You may use the fact that:
```math
A^{-1} = (P * L * U)^{-1} = U^{-1} * L^{-1} * P^{-1} = U^{-1} * L^{-1} * P
```

**Example:**

```math
A=\begin{pmatrix}
1 & 2 & 1 \\
1 & 2 & 2 \\
2 & 1 & 1
\end{pmatrix}
```

Expected result:

```math
A^{-1}=\begin{pmatrix}
0 & -\frac{1}{3} & \frac{2}{3}\\
1 & -\frac{1}{3} & -\frac{1}{3}\\
-1 & 1 & 0
\end{pmatrix}
```

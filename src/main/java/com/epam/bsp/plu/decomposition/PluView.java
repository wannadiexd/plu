package com.epam.bsp.plu.decomposition;

public class PluView {
	private double[][] perm;
	private double[][] lower;
	private double[][] upper;
	
	public PluView(double[][] perm, double[][] lower, double[][] upper) {
		super();
		this.perm = perm;
		this.lower = lower;
		this.upper = upper;
	}

	public double[][] getPerm() {
		return perm;
	}

	public double[][] getLower() {
		return lower;
	}

	public double[][] getUpper() {
		return upper;
	}
	
}

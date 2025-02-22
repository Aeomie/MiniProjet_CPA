package algorithms;

import java.awt.Point;

public class Edge {
	
	protected Point p,q;
	protected double weight;
	public Edge(Point p,Point q) {
		this.p = p;
		this.q = q;
		this.weight = calcDistance();
	}
	
	public double calcDistance() {
		return p.distance(q);
	}
	
	public double getWeight() {
		return weight;
	}
}

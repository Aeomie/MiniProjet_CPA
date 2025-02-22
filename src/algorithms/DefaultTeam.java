package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class DefaultTeam {

  public int[][] calculShortestPaths(ArrayList<Point> points, int edgeThreshold) {
	    int[][] paths=new int[points.size()][points.size()];
	    for (int i=0;i<paths.length;i++) for (int j=0;j<paths.length;j++) paths[i][j]=i;

	    double[][] dist=new double[points.size()][points.size()];

	    for (int i=0;i<paths.length;i++) {
	      for (int j=0;j<paths.length;j++) {
	        if (i==j) {dist[i][i]=0; continue;}
	        if (points.get(i).distance(points.get(j))<=edgeThreshold) dist[i][j]=points.get(i).distance(points.get(j));
	        else dist[i][j]=Double.POSITIVE_INFINITY;
	        
	        paths[i][j]=j;
	      }
	    }

	    for (int k=0;k<paths.length;k++) {
	      for (int i=0;i<paths.length;i++) {
	        for (int j=0;j<paths.length;j++) {
	          if (dist[i][j]>dist[i][k] + dist[k][j]){
	            dist[i][j]=dist[i][k] + dist[k][j];
	            paths[i][j]=paths[i][k];

	          }
	        }
	      }
	    }

	    return paths;
	  }
  public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		Kruskal k = new Kruskal();
		ArrayList<Edge> mstEdges = k.kruskal(hitPoints);

		Tree2D steinerTree= k.edgesToTree(mstEdges, hitPoints.get(0));

		int[][] paths = calculShortestPaths(points, edgeThreshold);
		
		Tree2D newSteinerTree = Steiner(paths, steinerTree, points);

		return newSteinerTree;
  }
  
	private Tree2D Steiner(int[][] paths, Tree2D steinerTree, ArrayList<Point> points) {
	    Point root = steinerTree.getRoot();
	    int rootIndex = points.indexOf(root);
	    ArrayList<Tree2D> children = new ArrayList<>();

	    for (Tree2D subTree : steinerTree.getSubTrees()) {
	        Point childRoot = subTree.getRoot();
	        int childIndex = points.indexOf(childRoot);
	        int nextNodeIndex = paths[rootIndex][childIndex];
	        Point nextNode = points.get(nextNodeIndex);

	        if (nextNode.equals(childRoot)) {
	            // Point next is the normal child node
	            children.add(Steiner(paths, subTree, points));
	        } else {
	            // Point next is an intermediate node
	            Tree2D intermediateTree = new Tree2D(nextNode, new ArrayList<>(List.of(subTree)));
	            children.add(Steiner(paths, intermediateTree, points));
	        }
	    }

	    return new Tree2D(root, children);
	}
  public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {

    int[][] paths = calculShortestPaths(points, edgeThreshold);
    int budget = 1664;
    
    ArrayList<Point> opti_hitPoints = calc_OptiHitPoints(points,hitPoints,paths,budget);
    
	Kruskal k = new Kruskal();
	ArrayList<Edge> mstEdges = k.kruskal(opti_hitPoints);

	Tree2D steinerTree= k.edgesToTree(mstEdges, hitPoints.get(0));

	Tree2D newSteinerTree = Steiner(paths, steinerTree, points);

	return newSteinerTree;
  }
  
  
  public ArrayList<Point> calc_OptiHitPoints(ArrayList<Point> points,ArrayList<Point> hitPoints, int[][] paths,int budget){
	  	System.out.println("calc Opti hit points running");
	    Point main_point = hitPoints.get(0);
	    
	    ArrayList<PointPath> hitPaths = new ArrayList<PointPath>();
	    
	    for(int i = 1 ; i < hitPoints.size(); i++) {
	    	Point startPoint = main_point;
	    	Point nextPoint = new Point(-1,-1); // point that cant exist
	    	Point endPoint = hitPoints.get(i);
	    	int end_point_Index = points.indexOf(hitPoints.get(i));
	    	double distance = 0;
	    	ArrayList<Point> path = new ArrayList<Point>();
	    	
	    	path.add(startPoint);
	    	while(!nextPoint.equals(endPoint)) {
	    		int nextPoint_Index = paths[points.indexOf(startPoint)][end_point_Index];
	    		nextPoint = points.get(nextPoint_Index);
	    		
	    		path.add(nextPoint);
	    		distance += startPoint.distance(nextPoint);
	    		
	    		startPoint = nextPoint;

	    	}
	    	hitPaths.add(new PointPath(hitPoints.get(i),path, distance));
	    }
	    
	    ArrayList<Point> opti_hitPoints = Greedy(hitPaths,budget);
	    System.out.println("calc Opti hit points finished");
	    System.out.println("opti hitPoint size : " + opti_hitPoints.size());
	    return opti_hitPoints;
  }
  
  public ArrayList<Point> Greedy(ArrayList<PointPath> paths, int budget){
	System.out.println("greedy running");
	HashMap<Double,PointPath> map = new HashMap<Double,PointPath>();
	for(PointPath path: paths) {
		Double ratioed = path.distance/budget;
		map.put(ratioed,path);
	}
    List<Double> weightList = new ArrayList<>(map.keySet());
    Collections.sort(weightList);
    
	double sum = 0;
	ArrayList<Point> resultPoints = new ArrayList<Point>();
	for(Double weight : weightList) {
		if(sum + weight <= 1) {
			sum += weight;
			resultPoints.add(map.get(weight).hitPoint);
		}
	}
	for(Double weight : map.keySet()) {
		System.out.println("hitPoint : " + map.get(weight).hitPoint);
		System.out.println("weight : " + weight);
	} 
	
	System.out.println("Result size : "+ resultPoints.size());
	System.out.println("greedy finished");
	return resultPoints;
  }
  
  class PointPath{
	  private Point hitPoint;
	  private ArrayList<Point> path;
	  double distance;
	  protected PointPath(Point hitPoint ,ArrayList<Point> path , double distance) {
		  this.hitPoint = hitPoint;
		  this.path = path;
		  this.distance = distance;
	  }
	  public ArrayList<Point> getPath() {
		  return path;
	  }
	  public double getDistance() {
		  return distance;
	  }
	  public Point getHitPoint() {
		  return hitPoint;
	  }
  }
}

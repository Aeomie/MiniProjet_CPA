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
	/**
	Functionality : Builds the tree following the shortestPaths ( paths )
	 */
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
	
	// INFORMATION
	/**
	 * Anything under 70 or above 155 will make this program loop indefintely
	 * i could add a mechanism to stop it , if distance is over budget and doesnt change in 10 loop
	 *but i cant open a popUp that says in the window that the parameter value is too low/high.
	 *This approach should require a button in the GUI , that lets you put the distance parametere between nodes that's needed. and it passes it to the code.
	 */
	
	// How it Works & Why this way
	/**
		I was inspirted by an Ai/ML Approach where we have some sorts of Clusters.
		
		We define the shortestPaths that are then used in our scoring algorithm to produce 
		hitPoints paths to other hitPoints with the number of hitpoints close to them.
		
		We loop aslong as the distance is higher than the budget , and remove the points that have the least points close to them
	 */
  public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {

    int[][] paths = calculShortestPaths(points, edgeThreshold);
    
    int DistanceBetweenNodes = 150;
    int budget = 1664;
    int distance = 100000; // instead of Initializing it as Integer MAX Value , i initialized it as a high value for the graph
    
    ArrayList<Point> hitPoints_copy = (ArrayList<Point>)hitPoints.clone();
    
    ArrayList<infoPoint> infoPoints = scoring_algorithm(points,hitPoints,paths,DistanceBetweenNodes);
    Collections.sort(infoPoints, (e1,e2) -> Integer.compare(e1.nb_closePoints, e2.nb_closePoints));
   
	Kruskal k = new Kruskal();
	ArrayList<Edge> mstEdges = new ArrayList<Edge>();
	Tree2D steinerTree= new Tree2D(new Point(0,0), new ArrayList<Tree2D>());
	
	Tree2D newSteinerTree = new Tree2D(new Point(0,0), new ArrayList<Tree2D>());
	Point main_point = hitPoints.get(0);
	
	int index = 0;
    while(distance > budget) {
    	//System.out.println("distance = " + distance);
    	Point removePoint = infoPoints.get(index).root;
    	if(removePoint.equals(main_point)) continue;
    	hitPoints_copy.remove(removePoint);
    	mstEdges = k.kruskal(hitPoints_copy);
    	steinerTree= k.edgesToTree(mstEdges, main_point);

    	newSteinerTree = Steiner(paths, steinerTree, points);
    	distance = (int)distanceCalculator(newSteinerTree);
    	index++;
    }

	//System.out.println("distance c: " + (int)distanceCalculator(newSteinerTree));
	return newSteinerTree;
  }
  
  
	/**
	Functionality : Calculates the distance of all the tree
	Use : Helps on looping while over the budget.
	 */
  public double distanceCalculator(Tree2D tree) {
	  double distance = 0;
	  Point root = tree.getRoot();
	  for(Tree2D subTree : tree.getSubTrees()) {
		  distance += root.distance(subTree.getRoot());
		  distance += distanceCalculator(subTree);
	  }
	  return distance;
  }
	/**
	Functionality : Calculates close hitPoints for each hitPoints while respecting the threshold
					Also initialises the path to that hitPoint
					
	Use : Helps us create a sort of cluster that helps us determine if that point is worth taking or not.
	 */
  public ArrayList<infoPoint> scoring_algorithm(ArrayList<Point> points, ArrayList<Point> hitPoints, int[][] shortestPaths, int threshold) {
	  ArrayList<infoPoint> infoPoints = new ArrayList<infoPoint>();
	  
	  ArrayList<PointPath> paths = new ArrayList<PointPath>();
	  int nbClose = 0;
	  for(Point start : hitPoints) {
		  paths.clear();
		  nbClose = 0; // reset
		  for(Point end : hitPoints) {
			  if(start.equals(end)) continue;
			  PointPath path = getPath(points,shortestPaths,start,end);
			  if(path.distance <= threshold) {
				  nbClose++;
			  }
			  paths.add(path);
		  }
		  infoPoint info = new infoPoint(start, paths, nbClose);
		  infoPoints.add(info);
	  }
	  
	  return infoPoints;
  }
	/**
	Functionality : gets the Path from Point  A to Point B , while following the intermediary points in paths.
					
	Use : We need it for our scoring_algorithm
	 */
  public PointPath getPath(ArrayList<Point> points, int[][] paths ,Point start, Point end) {
		Point startPoint = start;
    	Point nextPoint = new Point(-1,-1); // point that cant exist in our plan
    	Point endPoint = end;
    	int end_point_Index = points.indexOf(end);
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
    	return new PointPath(start, end, path, distance);
  }
	/**
	Functionality : Class that holds distance + paths.
	*/
  class PointPath{
	  private Point startPoint;
	  private Point endPoint;
	  private ArrayList<Point> path;
	  double distance;
	  protected PointPath(Point startPoint, Point endPoint ,ArrayList<Point> path , double distance) {
		  this.startPoint = startPoint;
		  this.endPoint = endPoint;
		  this.path = path;
		  this.distance = distance;
	  }
	  public ArrayList<Point> getPath() {
		  return path;
	  }
	  public double getDistance() {
		  return distance;
	  }
	  public Point getStartPoint() {
		  return startPoint;
	  }
	  public Point getEndPoint() {
		  return endPoint;
	  }
  }
	/**
	Functionality : Holder class
	Use: Holds data for our Hitpoint , with paths to other points and number of points that are close.
	*/
  class infoPoint{
	  private Point root;
	  private ArrayList<PointPath> paths;
	  private int nb_closePoints;
	  
	  infoPoint(Point root , ArrayList<PointPath> paths, int closePoints){
		  this.root = root;
		  this.paths = paths;
		  this.nb_closePoints = closePoints;
	  }
	  
	  public Point getRoot() {
		  return root;
	  }
	  public ArrayList<PointPath> getPaths(){
		  return paths;
	  }
	  public int getClosePointsCount() {
		  return nb_closePoints;
	  }
  }
}

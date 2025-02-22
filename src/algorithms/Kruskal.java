package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Kruskal {
	  public ArrayList<Edge> kruskal(ArrayList<Point> points) {
		  
		  ArrayList<Edge> edges = new ArrayList<Edge>();
		  for (Point p: points) {
		      for (Point q: points) {
		        if (p.equals(q) || contains(edges,p,q)){
		          continue;
		        }
		        edges.add(new Edge(p,q));
		      }
		    }
	      // Sort edges by weight
	      Collections.sort(edges, (e1, e2) -> Double.compare(e1.weight, e2.weight));

	      // Union-Find initialization
	      Map<Point, Point> parent = new HashMap<>();
	      Map<Point, Integer> size = new HashMap<>();
	      for (Point p : points) {
	          parent.put(p, p);
	          size.put(p, 1);
	      }

	      ArrayList<Edge> mstEdges = new ArrayList<>();
	      int edgesAdded = 0;
	      int n = points.size();

	      // Kruskal's algorithm
	      for (Edge edge : edges) {
	          if (edgesAdded == n - 1) break;

	          Point rootP = find(parent, edge.p);
	          Point rootQ = find(parent, edge.q);

	          if (rootP != rootQ) {
	              mstEdges.add(edge);
	              union(parent, size, rootP, rootQ);
	              edgesAdded++;
	          }
	      }

	      return mstEdges;
	  }
	  private Point find(Map<Point, Point> parent, Point p) {
	      Point root = p;
	      while (!parent.get(root).equals(root)) {
	          root = parent.get(root);
	      }
	      // Path compression
	      while (!p.equals(root)) {
	          Point next = parent.get(p);
	          parent.put(p, root);
	          p = next;
	      }
	      return root;
	  }

	  private void union(Map<Point, Point> parent, Map<Point, Integer> size, Point rootP, Point rootQ) {
	      if (size.get(rootP) < size.get(rootQ)) {
	          parent.put(rootP, rootQ);
	          size.put(rootQ, size.get(rootQ) + size.get(rootP));
	      } else {
	          parent.put(rootQ, rootP);
	          size.put(rootP, size.get(rootP) + size.get(rootQ));
	      }
	  }

	  public Tree2D edgesToTree(ArrayList<Edge> edges, Point root) {
		    ArrayList<Edge> remainder = new ArrayList<Edge>();
		    ArrayList<Point> subTreeRoots = new ArrayList<Point>();
		    Edge current;
		    while (edges.size()!=0) {
		      current = edges.remove(0);
		      if (current.p.equals(root)) {
		        subTreeRoots.add(current.q);
		      } else {
		        if (current.q.equals(root)) {
		          subTreeRoots.add(current.p);
		        } else {
		          remainder.add(current);
		        }
		      }
		    }

		    ArrayList<Tree2D> subTrees = new ArrayList<Tree2D>();
		    for (Point subTreeRoot: subTreeRoots) subTrees.add(edgesToTree((ArrayList<Edge>)remainder.clone(),subTreeRoot));

		    return new Tree2D(root, subTrees);
		  }
	  private boolean contains(ArrayList<Edge> edges,Point p,Point q){
		    for (Edge e:edges){
		      if (e.p.equals(p) && e.q.equals(q) ||
		          e.p.equals(q) && e.q.equals(p) ) return true;
		    }
		    return false;
		  }
}

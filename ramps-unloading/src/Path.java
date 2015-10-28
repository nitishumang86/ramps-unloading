import java.util.*;
import java.math.*;

public class Path {
	int pathId;
	ArrayList<Location> locationSet = new ArrayList<Location>(); 
	// an ordered list of locations that the path consists of 0-1-2-3 
	// In our problem, the first location is the pad location
	
	double pathLength; // path length = 0-1 + 1-2 + 2-3 + 3-0
	
	public double CalculatePathLength(){
		// path length excluding the first and last arcs connecting the path to the pad
		
		double length = 0.0; 
		//length+= Math.sqrt(Math.pow(this.locationSet.get(0).x_coord,2) + Math.pow(this.locationSet.get(0).y_coord, 2)); 
		for(int i=1; i<this.locationSet.size(); i++){
			double term1 = Math.pow((this.locationSet.get(i).x_coord - this.locationSet.get(i-1).x_coord),2);
			double term2 = Math.pow((this.locationSet.get(i).y_coord - this.locationSet.get(i-1).y_coord),2);
			length+= Math.sqrt(term1+term2);
		}
		
		// return to pad location from the last node
		int lastIndex = locationSet.size()-1;
		double term1 = Math.pow((this.locationSet.get(0).x_coord - this.locationSet.get(lastIndex).x_coord),2);
		double term2 = Math.pow((this.locationSet.get(0).y_coord - this.locationSet.get(lastIndex).y_coord),2);
		length+= Math.sqrt(term1+term2);
		
		return length;
	}



}

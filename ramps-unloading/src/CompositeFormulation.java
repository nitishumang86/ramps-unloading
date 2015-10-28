import gurobi.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import java.math.*;

public class CompositeFormulation {

	String read = "C:\\Users\\212478865\\Desktop\\Ramps\\Testing\\InputData.txt"; // path for the input file 
	String outputfile = "C:\\Users\\212478865\\Desktop\\Ramps\\Testing\\Output.txt"; // path for the output file 
	String pathfile = "C:\\Users\\212478865\\Desktop\\Ramps\\Testing\\Paths.txt"; // path for the output file 
	
	int numpads;
	int numVins;
	int numVinHandlers; // assume a single crew
	int numLocations;
	static int numPaths;
	static int numTrips;
	
	static ArrayList<Pad> padList = new ArrayList<Pad>();
	static ArrayList<Location> locationList = new ArrayList<Location>();
	static ArrayList<ArrayList<Path>> pathList = new ArrayList<ArrayList<Path>>();
	// each array list of paths for a given pad location
	
	public static void main(String[] args) {
		CompositeFormulation obj = new CompositeFormulation();
		obj.ReadData();
		
		int pathsGenerated = 0;
		for(int i=0; i<padList.size(); i++){
			for (int p = 0; p < pathList.size(); p++) 
				pathsGenerated = pathsGenerated + pathList.get(p).size();
			pathList.add(obj.GeneratePaths(locationList, padList.get(i),pathsGenerated));
			
		}	
		// assuming that all locations can be used by every pad
		//System.out.println("Number of paths generated:"+(pathList.size()*pathList.get(0).size()));
		
		// print all paths
		
		/*for (int i = 0; i < pathList.size(); i++) {
	    	  for (int j = 0; j < pathList.get(i).size(); j++) {
	    		  System.out.print("Path "+pathList.get(i).get(j).pathId+": \t");
	    		  ArrayList<Location> locSet = pathList.get(i).get(j).locationSet;
	    		  for(int l=0; l<locSet.size(); l++)
	    			  System.out.print(locSet.get(l).x_coord+" "+locSet.get(l).y_coord+"\t");
	    		  double lengthPath = pathList.get(i).get(j).pathLength;
	    		  System.out.print(lengthPath);
	    		  System.out.println();
	    	  }
	     }*/    
		
		// calculate the number of paths
		for (int i = 0; i < pathList.size(); i++) 
	    	 numPaths = numPaths + pathList.get(i).size();
	    
		
		obj.WritePathFile();
		obj.Solve();
	}

	public void ReadData(){
		try{
			File f = new File(read);
			Scanner sc = new Scanner(f);
			
			String line = sc.nextLine(); 
			String[] temp = line.split("\t");
			numpads = Integer.parseInt(temp[1]);
			
			line = sc.nextLine();// escape header row 
			for(int j=0; j<numpads; j++){
				line = sc.nextLine(); 
				temp = line.split("\t");
				Pad pObj = new Pad();
				pObj.padId = Integer.parseInt(temp[0]);
				pObj.x_coord = Integer.parseInt(temp[1]);
				pObj.y_coord = Integer.parseInt(temp[2]);
				String[] tempS = temp[3].split(" ");
				for(int i=0; i<temp[3].length(); i++){
					Vehicle vObj = new Vehicle();
					vObj.vinId =  Integer.parseInt(tempS[0]);
					//vObj.vinPad =  pObj;
					pObj.vinList.add(vObj);
				}
				padList.add(pObj);
			}
			
			line = sc.nextLine(); 
			temp = line.split("\t");
			numVins = Integer.parseInt(temp[1]);
			
			line = sc.nextLine(); 
			temp = line.split("\t");
			numVinHandlers = Integer.parseInt(temp[1]);
		
			line = sc.nextLine(); 
			temp = line.split("\t");
			numLocations = Integer.parseInt(temp[1]);
			
			line = sc.nextLine();// escape header row 
			for(int i=0; i<numLocations; i++){
				line = sc.nextLine(); 
				temp = line.split("\t");
				Location lObj = new Location();
				lObj.locationId = Integer.parseInt(temp[0]);
				lObj.x_coord = Integer.parseInt(temp[1]);
				lObj.y_coord = Integer.parseInt(temp[2]);
				lObj.capacity = Integer.parseInt(temp[3]);
				locationList.add(lObj);
			}
		
			numTrips = (int)Math.ceil((double)numVins/numVinHandlers);
		}
		catch (IOException e) {
			System.out.println(e);
		}
	}
	
	/** return the (straight line) distance between locations l1 and l2 */
	public double Distance(Location l1, Location l2){
		double length;
		double term1 = Math.pow((l1.x_coord - l2.x_coord),2);
		double term2 = Math.pow((l1.y_coord - l2.y_coord),2);
		length= Math.sqrt(term1+term2);
		return length;
	}
	
	public int Factorial(int n){
		int prod = 1;
		for(int i=1; i<=n; i++)
			prod = prod*i;
		return prod;
	}
	
	/**p1 and p2 are lists of indices (note:0-2-3-4 and 0-4-3-2 are the same, where 0 
	is the pad location)
	*/
	public boolean SameSequence(ArrayList<Integer> p1, ArrayList<Integer> p2){
		boolean isSame = true;
		if(p1.size()!=p2.size()) // not of the same size
			isSame = false;
		
		/*else if (p1.get(1)==p2.get(1)){ // comparing 0-1-2-3-4 with 0-1-2-3-4
			for(int i=0; i<p1.size(); i++){
				if(p1.get(i)!=p2.get(i)){
					isSame = false;
					break;
				}
			}
		}*/
		
		/*else if ((p1.get(0)==p2.get(0)) && p1.get(1)==p2.get(p2.size()-1)){ 
		// comparing 0-1-2-3-4 with 0-4-3-2-1
			for(int i=1; i<p1.size(); i++){
				if(p1.get(i)!=p2.get(p2.size()-i)){
					isSame = false;
					break;
				}
			}
		}*/
		else{
			for(int i=0; i<p1.size(); i++){
				if(p1.get(i)!=p2.get(i)){
					isSame = false;
					break;
				}
			}
		}
		return isSame;
	}
	
	
	/** given value of n, return factorial (n) sequences with first element as i
	 * for example for n=3 and i=1, return following sequences: 1-2-3, 2-3-1, 2-1-3, 1-3-2, 3-1-2, 3-2-1 
	 * */
	public ArrayList<ArrayList<Integer>> GenerateAllPossibleSequences(int n, int i){
		
		ArrayList<ArrayList<Integer>> sequences = new ArrayList<ArrayList<Integer>>();
		
		if(n==1){
			ArrayList<Integer> seq1 = new ArrayList<Integer>();
			seq1.add(i); 
			sequences.add(seq1); 
			return sequences;
		}
		
		if(n==2){
			ArrayList<Integer> seq1 = new ArrayList<Integer>();
			seq1.add(i); seq1.add(i+1);
			ArrayList<Integer> seq2 = new ArrayList<Integer>();
			seq2.add(i+1); seq2.add(i);
			sequences.add(seq1); sequences.add(seq2);
			return sequences;
		}
		else{ //if n>=3
			ArrayList<ArrayList<Integer>> sequences_new = GenerateAllPossibleSequences((n-1), (i+1));
			for(int j=0; j<sequences_new.size(); j++){
				for(int k=0; k<=sequences_new.get(j).size(); k++){
					ArrayList<Integer> seq_j = new ArrayList<Integer>();
					seq_j.addAll(sequences_new.get(j));
					seq_j.add(k,i);
					sequences.add(seq_j);
				}
			}
			return sequences;
		}
		
	}
	
	/** return minimum (closed loop) path by simple enumeration for a list of nodes 
	where the first node is the depot (pad location)
	*/
	public Path ShortestPath(ArrayList<Location> lList){
		
		/* generate the list of all possible sequences of location ids */
		//System.out.println("SequencesEvaluated.size():"+sequencesEvaluated.size());
		//System.out.println("Max number of sequences to be evaluated:"+Factorial(lList.size()-1));
		/*while(sequencesEvaluated.size()<(Factorial(lList.size()-1))){
			ArrayList<Integer> newSequence = new ArrayList<Integer>();
			newSequence.add(0); // path starts at 0 index corresponding to the pad
			while(newSequence.size()<lList.size()){
				Random r = new Random();
				int index = 1+r.nextInt(lList.size()-1); //
				//System.out.println("Index value:"+index);
				while(newSequence.contains(index))
					index = 1+r.nextInt(lList.size()-1);
				newSequence.add(index); // location corresponding to that index
			}
			
			boolean isEvaluated = false;
			for(int i=0; i<sequencesEvaluated.size(); i++){
				if(SameSequence(newSequence, sequencesEvaluated.get(i))){
					//System.out.println("EVALUATED");
					isEvaluated = true;
					break;
				}
			}
			if(!isEvaluated)
				sequencesEvaluated.add(newSequence);
			//System.out.println("Number of sequences evaluated:"+sequencesEvaluated.size());
			
		}*/
		ArrayList<ArrayList<Integer>> sequencesEvaluated = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> sequences = GenerateAllPossibleSequences((lList.size()-1), 1);
		for(int i=0; i<sequences.size(); i++){
			ArrayList<Integer> seq_i = new ArrayList<Integer>();
			seq_i.addAll(sequences.get(i));
			seq_i.add(0,0); // add pad location
			sequencesEvaluated.add(seq_i);
		}
		
		
		/*
		 * given the list of sequences compute the minimum path length by simply taking the
		 * euclidean distances
		 */
		Path returnPathObject = new Path();
		double minPathLength=100000;
		//System.out.println("The number of sequences evaluated:"+sequencesEvaluated.size());
		
		for(int i=0; i<sequencesEvaluated.size(); i++){
			
			Path pathObject = new Path();
			ArrayList<Location> lsp = new ArrayList<Location>();
			
			// add the pad location as the first node in the "location set"
			Location lPad = new Location();
			lPad.x_coord = lList.get(0).x_coord;
			lPad.y_coord = lList.get(0).y_coord;
			lsp.add(lPad);
			
			for(int j=1; j<sequencesEvaluated.get(i).size(); j++){ 
				// j=0 corresponds to the pad location
				int index = sequencesEvaluated.get(i).get(j);
				lsp.add(lList.get(index));
			}
			
			
			//System.out.println("Size of location set of generated path:"+lsp.size());			
			pathObject.locationSet = lsp;
			/*for(int l=0; l<lsp.size(); l++)
			  System.out.print(lsp.get(l).x_coord+" "+lsp.get(l).y_coord+"\t");
			System.out.println();
			*/
			pathObject.pathLength = pathObject.CalculatePathLength();
			//System.out.println("Shortest path is:"+pathObject.pathLength);
			if (pathObject.pathLength<minPathLength){
				minPathLength = pathObject.pathLength;
				returnPathObject = pathObject;
			}
			
		}
		
		return returnPathObject;
	}
	
	/** generate the set of all paths for a given set of locations (EXCLUDING the pad)
	 * and given pad (path can have at most (k+1) nodes where k is the number of vin 
	 * handlers)
	*/	
	public ArrayList<Path> GeneratePaths(ArrayList<Location> lList, Pad padLocation, int pathCounter){
		System.out.println("Number of locations......"+lList.size());
		System.out.println("Generating paths.......");
		ArrayList<Path> pathPadList = new ArrayList<Path>();
		
		for(int i=1; i<=numVinHandlers; i++){
			ArrayList<ArrayList<Location>> listCombinations = Combinations(lList,i);
			System.out.println("number of combinations for "+i+" vinhandler(s):"+listCombinations.size());
			for(int j=0; j<listCombinations.size(); j++){
				
				ArrayList<Location> newLocList = new ArrayList<Location>();
				
				/** first location is the pad location */
				Location lObj = new Location();
				lObj.x_coord = padLocation.x_coord;
				lObj.y_coord = padLocation.y_coord;
				newLocList.add(lObj);
				newLocList.addAll(listCombinations.get(j));
				
				System.out.println("Input array for shortest path function:");
				for(int k=0; k<newLocList.size(); k++){
					System.out.println(newLocList.get(k).x_coord+" "+newLocList.get(k).y_coord);
				}
				
				/* for given set of locations (including pad) generate the shortest path */
				Path pObj = ShortestPath(newLocList);
				System.out.println("Shortest path is:"+pObj.pathLength);
				pObj.pathId = pathCounter;
				pathPadList.add(pObj);
				
				pathCounter++;
				System.out.println("The number of paths generated:"+pathCounter);
			}
		}
		
		return pathPadList;	
	
	}
	
	/** given a list of locations and a number K, generate all possible combinations of 
	 * locations
	 * */
	public ArrayList<ArrayList<Location>> Combinations(ArrayList<Location> lList, int KValue){
		ArrayList<ArrayList<Location>> listCombinations = new ArrayList<ArrayList<Location>>();
		
		
		if(KValue==1){
			
			for(int i=0; i<lList.size(); i++){
				ArrayList<Location> tempLocList = new ArrayList<Location>();
				tempLocList.add(lList.get(i));
				listCombinations.add(tempLocList);
			}
			return listCombinations;	
			
		}
		
		else { // if  K > 1
			
			for(int i=0; i<lList.size()-KValue+1; i++){
				int startIndex = i;
				
				ArrayList<Location> tempLocList = new ArrayList<Location>(); 
				for(int j=i+1;j<lList.size(); j++){
					//if(j!=i)
						tempLocList.add(lList.get(j));
				}
				//KCopy=KCopy-1;
				ArrayList<ArrayList<Location>> combList = Combinations(tempLocList, KValue-1);
				// generate all combinations of size (KCopy-1)
					
				for(int j=0;j<combList.size(); j++){
					
					ArrayList<Location> newList = new ArrayList<Location>();
					newList.add(lList.get(startIndex));
					
					for(int k=0; k<combList.get(j).size(); k++)
						newList.add(combList.get(j).get(k));
					
					listCombinations.add(newList);
				}
				
			}
			
			return listCombinations;
		}
	
	}
	
	/** Does given location belong to the given path */
	public int LocationPath(Path pObj, Location lObj){
		int beta = 0;
		for(int l=0; l<pObj.locationSet.size(); l++){
			if(pObj.locationSet.get(l).x_coord == lObj.x_coord && 
					pObj.locationSet.get(l).y_coord == lObj.y_coord){
				beta = 1;
				break;
			}
		}
		
		
		return beta;
	}
	
	public void WritePathFile(){
		try {
			 FileWriter fw = new FileWriter(pathfile);
			 
			 // print all paths
			 for (int i = 0; i < pathList.size(); i++) {
		    	  for (int j = 0; j < pathList.get(i).size(); j++) {
		    		  fw.write("Path "+pathList.get(i).get(j).pathId+": \t");
		    		  ArrayList<Location> locSet = pathList.get(i).get(j).locationSet;
		    		  double lengthPath = pathList.get(i).get(j).pathLength;
		    		  for(int l=0; l<locSet.size(); l++)
		    			  fw.write(locSet.get(l).x_coord+" "+locSet.get(l).y_coord+"\t");
		    		  fw.write("Path length:"+lengthPath);
		    		  fw.write("\n");
		    	  }
		      }
		      fw.close(); 
		}
	
    	catch(IOException e) {
    		System.out.println(e);
    	}
	}
	
	public void Solve(){
		try {
		      GRBEnv    env   = new GRBEnv("mip1.log");
		      GRBModel  model = new GRBModel(env);
		      model.set(GRB.StringAttr.ModelName, "ramps_unordered");
		      
		      System.out.println("Solving the model........");
		      
		      
		      /** Create variables */

		      /* x_ilt assignment variables */
		      GRBVar[][][] x = new GRBVar[numVins][numLocations][numTrips];
		      for (int i = 0; i < numVins; i++) {
		    	  for (int l = 0; l < numLocations; l++) {
		    		  for(int t=0; t<numTrips; t++){
		    			  x[i][l][t] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, 
		    					  "assignVin_"+String.valueOf(i)+"_toLocation_"+String.valueOf(l)+"_inTrip_"+String.valueOf(t));	
		    			 
		    		  }
		    	  }
		      }
		      
		      /* z_pt trip path assignment variables  */
		      GRBVar[][] z = new GRBVar[numPaths][numTrips];
		      int p=0;
		      for (int i = 0; i < pathList.size(); i++) {
		    	  for (int j = 0; j < pathList.get(0).size(); j++) {
		    		  for (int t = 0; t < numTrips; t++) {	
		    		  
		    			  z[p][t] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, 
		    				  				"assign_path "+String.valueOf(p)+" "
		    				  				+" to trip "+String.valueOf(t));	
		    	  
		    		  }
		    		  p++;
		    	  }
		      }
		      
		      // Integrate new variables
		      model.update();
		      
		       
		      
		      /** Objective function - Set objective: minimize summation (over p and t) of d_p z_pt*/
		     
		      GRBLinExpr expr = new GRBLinExpr();
		      p=0;
		      for (int i = 0; i < pathList.size(); i++) {
		    	  for (int j = 0; j < pathList.get(0).size(); j++) {
		    		  for (int t = 0; t < numTrips; t++) {	
		    		  
		    			  expr.addTerm(pathList.get(i).get(j).pathLength, z[p][t]);
		    		  }
		    		  p++;
		    	  }
		      }
		      model.setObjective(expr, GRB.MINIMIZE);
		      
		      
		      /** Constraints */
		      
		      /* capacity constraints */ 
		      for (int l = 0; l < numLocations; l++) {	
		    	  expr = new GRBLinExpr();
		    	  for (int i = 0; i < numVins; i++) {
		    		  for(int t=0; t<numTrips; t++){
		    			  expr.addTerm(1.0, x[i][l][t]);
		    		  }
		    	  }
		    	  model.addConstr(expr, GRB.LESS_EQUAL, locationList.get(l).capacity,
	    				  "capacityConst"+String.valueOf(l));
		      }

		      
		      /* unique assignment constraints */ 
		      for (int i = 0; i < numVins; i++) {
		    	  expr = new GRBLinExpr();
		    	  for (int l = 0; l < numLocations; l++) {
		    		  for(int t=0; t<numTrips; t++){
		    			  expr.addTerm(1.0, x[i][l][t]);
		    		  }
		    	  }
		    	  model.addConstr(expr, GRB.EQUAL, 1,"uniqueAssign"+String.valueOf(i));
		      }
		       
		      
		      /* limit max number of vins to each trip (cannot be more than the number of vinhandlers) */ 
		      for(int t=0; t<numTrips; t++){
		    	  expr = new GRBLinExpr();
		    	  for (int i = 0; i < numVins; i++) {
		    	    for (int l = 0; l < numLocations; l++) {
		    			  expr.addTerm(1.0, x[i][l][t]);
		    		  }
		    	  }
		    	  model.addConstr(expr, GRB.LESS_EQUAL, numVinHandlers,"maxVinsToPath"+String.valueOf(t));
		      }
		      
		      
		     /* each trip is assigned to exactly one path */
		     for (int t = 0; t < numTrips; t++) {	
		    	  expr = new GRBLinExpr();
		    	  p=0;
		    	  for (int i = 0; i < pathList.size(); i++) {
			    	  for (int j = 0; j < pathList.get(0).size(); j++) {
			    		  
		    			  expr.addTerm(1.0, z[p][t]);
		    			  p++;
		    		  }
		    		  
		    	  }
		    	  model.addConstr(expr, GRB.EQUAL, 1,"uniquePathAssignTrip"+String.valueOf(t));
		      }
		      
		      /* linking x and z variables - 1*/
		     for (int t = 0; t < numTrips; t++) {
		    	 p=0;
		    	  for (int i = 0; i < pathList.size(); i++) {
			    	  for (int j = 0; j < pathList.get(0).size(); j++) {
			    		  
			    		  for (int l = 0; l < numLocations; l++) {
			    			  if(LocationPath(pathList.get(i).get(j), locationList.get(l))==1){
			    				  expr = new GRBLinExpr();
			    				  
			    				  for(int m=0; m<numVins; m++)
			    					  expr.addTerm(1.0,  x[m][l][t]);
			    				  expr.addTerm(-1.0,  z[p][t]);
			    				  model.addConstr(expr, GRB.GREATER_EQUAL, 0,"linkingXtoZ"+String.valueOf(p)+" "
			    						  +String.valueOf(t));
			    			  }
			    		  }
			    		  p++;
			    	  }
		    	  }
		    	 
		     }
		     
		     
		     /* linking x and z variables - 2*/
		     for (int m = 0; m < numVins; m++) {
		    	  for (int l = 0; l < numLocations; l++) {
		    		  for (int t = 0; t < numTrips; t++) {
		    			  expr = new GRBLinExpr();
		    			  expr.addTerm(-1.0, x[m][l][t]);
		    			  p=0;
		    			  for (int i = 0; i < pathList.size(); i++) {
					    	  for (int j = 0; j < pathList.get(0).size(); j++) {
					    		
					    		  expr.addTerm(LocationPath(pathList.get(i).get(j), locationList.get(l)),  z[p][t]);
					    		  
					    		  p++;  
					    	  }
		    			  
		    			  }
		    		  
		    			  model.addConstr(expr, GRB.GREATER_EQUAL, 0,"linkingXtoZ2"+String.valueOf(m)+" "
		    					  +" "+String.valueOf(l)+" "+String.valueOf(t));
		    		  
		    		  }
		    	  }
		     }
		     
		     /* path variables are less than equal to 1 if less than 2 locations in the path */
		     p=0;
		      for (int i = 0; i < pathList.size(); i++) {
		    	  for (int j = 0; j < pathList.get(i).size(); j++) {	
		    		  if(pathList.get(i).get(j).locationSet.size()>2){ // first location is the pad location
		    			  expr = new GRBLinExpr();
		    			  for (int t = 0; t < numTrips; t++) 
		    				  expr.addTerm(1.0, z[p][t]);  
		    			  model.addConstr(expr, GRB.LESS_EQUAL, 1,"binary_path_"+String.valueOf(p));
		    		  }
		    		  p++;
		    	  }
		      }
		    	
		     
		     System.out.println("max number of trips:"+Math.ceil((double)numVins/numVinHandlers));
		     /** Optimize */
		     model.getEnv().set(GRB.DoubleParam.TimeLimit, 600.0);
		     model.optimize();
		       
		     PrintSolution(model, x, z);
		      
		     int status = model.get(GRB.IntAttr.Status);
		     if (status == GRB.Status.UNBOUNDED) {
		       System.out.println("The model cannot be solved "
		     + "because it is unbounded");
		       return;
		     }
		     if (status == GRB.Status.OPTIMAL) {
		       System.out.println("The optimal objective is " +
		           model.get(GRB.DoubleAttr.ObjVal));
		       return;
		     }
		     if (status == GRB.Status.INF_OR_UNBD ||
		         status == GRB.Status.INFEASIBLE    ){
		       System.out.println("Optimization was stopped with status " + status);
		       return;
		     }

		      
		     /** Dispose of model and environment */
		     model.dispose();
		     env.dispose();

		    
			} catch (GRBException e) {
		      System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
		      e.printStackTrace();
		    }
			
	
	}

	public void PrintSolution(GRBModel model, GRBVar[][][] x, GRBVar[][] z) throws GRBException {
		/** write output file */
		
	      try{
	    	  FileWriter fw = new FileWriter(outputfile);
	    	  
	    	  // print all x variable values
	    	  for (int i = 0; i < numVins; i++) {
		    	  for (int l = 0; l < numLocations; l++) {
		    		  for(int t=0; t<numTrips; t++){
		    			  if(x[i][l][t].get(GRB.DoubleAttr.X)==1){
		    				  fw.write(x[i][l][t].get(GRB.StringAttr.VarName)+ " " +
	    					  						x[i][l][t].get(GRB.DoubleAttr.X));
		    				  fw.write("\n");
		    			  }
		    		  }
		    	  }
	    	  }
	    	  
	    	  // print all z variable values
	    	 
	    	  int p=0;
		      for (int i = 0; i < pathList.size(); i++) {
		    	  for (int j = 0; j < pathList.get(0).size(); j++) {
		    		  for (int t = 0; t < numTrips; t++) {	
		    		  
		    			  if(z[p][t].get(GRB.DoubleAttr.X)==1){
		    				  fw.write(z[p][t].get(GRB.StringAttr.VarName)+ " " +
		    					  						z[p][t].get(GRB.DoubleAttr.X));
		    				  fw.write("\n");	
		    			  }
		    		  }
		    		  p++;
		    	  }
		      }
	    	  
	    	  
	    	  // print objective function value
	    	  fw.write("Total distance traveled (Objective function value): " + 
	    	  model.get(GRB.DoubleAttr.ObjVal));
	    	
	    	  // print computation time
	    	  fw.write("\n");
	    	  fw.write("Total computation time): " + 
	    	    	  model.get(GRB.DoubleAttr.Runtime));
	    	  
	    	  // print optimality gap
	    	  fw.write("\n");
	    	  fw.write("Optimality Gap: " + 
	    	    	  model.get(GRB.DoubleAttr.MIPGap));
	    	  fw.close();
	    	  
	      }
	      catch(IOException e) {
	    	  System.out.println(e);
	      }
	      
	      
	}
}

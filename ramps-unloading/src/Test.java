import java.util.ArrayList;

public class Test {
	
	static ArrayList<Location> locationList = new ArrayList<Location>();
	
	public static void main(String[] args) {
		Test obj = new Test();
		//obj.Initialize();
		System.out.println(Math.ceil(5/(double)2));
	}

	public void Initialize(){
		/*Location l1 = new Location();
		l1.x_coord = 2;
		l1.y_coord = 3;
	
		Location l2 = new Location();
		l2.x_coord = 3;
		l2.y_coord = 4;
	
		Location l3 = new Location();
		l3.x_coord = 4;
		l3.y_coord = 5;
	
		Location l4 = new Location();
		l4.x_coord = 5;
		l4.y_coord = 6;
	
		Location l5 = new Location();
		l5.x_coord = 6;
		l5.y_coord = 7;
		
		Location l6 = new Location();
		l6.x_coord = 7;
		l6.y_coord = 8;
		
		locationList.add(l1); locationList.add(l2); locationList.add(l3); locationList.add(l4);
		locationList.add(l5); locationList.add(l6);
		ArrayList<ArrayList<Location>> listCombinations = Combinations(locationList, 3);
	
		System.out.println("HEYY");
		System.out.println("Size of listCombinations:"+listCombinations.size());
		for(int i=0; i<listCombinations.size(); i++){
			System.out.println("Printing combination "+(i+1));
			for(int j=0; j<listCombinations.get(i).size(); j++){
				System.out.println(listCombinations.get(i).get(j).x_coord + " " +
						listCombinations.get(i).get(j).y_coord);
			}
		}*/
		/*
		ArrayList<Integer> list1 = new ArrayList<Integer>();
		list1.add(0);list1.add(1);list1.add(2);
		ArrayList<Integer> list2 = new ArrayList<Integer>();
		list2.add(0);list2.add(2);list2.add(1);
		System.out.println("Sequences are same?:"+SameSequence(list1, list2));
		*/
		
		ArrayList<ArrayList<Integer>> lists = GenerateAllPossibleSequences(3,2);
		for(int i=0; i<lists.size(); i++){
			for(int j=0; j<lists.get(i).size(); j++){
				System.out.print(lists.get(i).get(j)+" ");
			}
			System.out.println();
		}
		System.out.println("Number of sequences:"+lists.size());
		
		
		
	
	}
	
	public ArrayList<ArrayList<Location>> Combinations(ArrayList<Location> lList, int KValue){
		//int KCopy = K;
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
}

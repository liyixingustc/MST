/*
CSE6140 HW1
This is an example of how your experiments should look like.
Feel free to use and modify the code below, or write your own experimental code, 
as long as it produces the desired output.
*/

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.ArrayList;

public class RunExperiments{


	public static void QuickSort(int[][] arr, int low, int high){
		if (arr == null || arr.length == 0)
			return;
 
		if (low >= high)
			return;
 
		// pick the pivot
		int middle = low + (high - low) / 2;
		int pivot = arr[middle][2];
 
		// make left < pivot and right > pivot
		int i = low, j = high;
		while (i <= j) {
			while (arr[i][2] < pivot) {
				i++;
			}
 
			while (arr[j][2] > pivot) {
				j--;
			}
 
			if (i <= j) {
				int temp1 = arr[i][0];
				int temp2 = arr[i][1];
				int temp3 = arr[i][2];
				arr[i][0] = arr[j][0];
				arr[i][1] = arr[j][1];
				arr[i][2] = arr[j][2];
				arr[j][0] = temp1;
				arr[j][1] = temp2;
				arr[j][2] = temp3;
				i++;
				j--;
			}
		}
 
		// recursively sort two sub parts
		if (low < j)
			QuickSort(arr, low, j);
 
		if (high > i)
			QuickSort(arr, i, high);
	}





	public static int computeMST(EdgeWeightedGraph G, int[][] new_weight_info){

		//sort
		int vortexnum = G.V();
		int edgenum = G.E();
		int[][] weight_info = new int [edgenum][3];
		int temp1, temp2, temp3;
		int MSTweight = 0;
		int key;

		//sort by building an array
		int num1 = 0;
		int num2 = 0;
		for (int v = 0; v < vortexnum; v++) {
			num2 = 0;
            for (Edge e : G.adj[v]) {
            	int ebgn = e.v;
        		int eend = e.w;
            	if(ebgn == v){
                weight_info[num1][0] = v;
                weight_info[num1][1] = num2;
                weight_info[num1][2] = e.weight();
                num1++;
                num2++;
            	}
            }
        }

 		QuickSort(weight_info, 0, edgenum - 1);

  		//greedy
  		UF uf = new UF(G.V());
  		num1 = 0;
  		int temp4 = 0;
  		//System.out.println(mst.E());
        while (temp4 < G.V() - 1) {
        	num2 = weight_info[num1][0];
        	//System.out.println(num2);
        	temp1 = 0;
        	for (Edge e : G.adj[num2]){
        		if(temp1 == weight_info[num1][1]){
        			//System.out.println(weight_info[num1][2]);
            		int ebgn = e.v;
        			int eend = e.w;
        			if (!uf.connected(ebgn, eend)) { // v-w does not create a cycle
                		uf.union(ebgn, eend);  // merge v and w components
                		new_weight_info[temp4][0] = e.v;
                		new_weight_info[temp4][1] = e.w;
                		new_weight_info[temp4][2] = e.weight;
                		//mst.addEdge(e.v, e.w, e.weight); // add edge e to mst
               			MSTweight += e.weight();
               			temp4++;
               			
               			//System.out.println(MSTweight);
            		}
        			break;
        		}
        		temp1++;
        	}
            num1++;
        }

        //QuickSort(new_weight_info, 0, vortexnum - 2);

		return MSTweight;

	}

	public static void isConnected (int[][] new_weight_info, int[] visited, int start, int end, Stack sta){
		int edgei_n = -1;
		int vortexnum = new_weight_info.length;
		int vortexi;
//		System.out.println(start);
//		System.out.println(end);
		visited[start] = 1;
		sta.push(start);
		if(!sta.isEmpty()){	
			if((int)sta.peek() == end){
				return;
			}else{
				for (int i = 0; i < vortexnum; i++){
					//if((i != start) && ((new_weight_info[i][0] == start) || (new_weight_info[i][1] == start))){
					if(((new_weight_info[i][0] == start) || (new_weight_info[i][1] == start))){
						if(new_weight_info[i][0] == start){
							vortexi = new_weight_info[i][1];
						}else{
							vortexi = new_weight_info[i][0];
						}
					//if((i != start) && (mst.adj[start][i] != 0)){
						if(vortexi == end){
							sta.push(vortexi);
							return;
						}
						if(visited[vortexi] == 0){
							visited[vortexi] = 1;
							isConnected(new_weight_info, visited, vortexi, end, sta);
							if((int)sta.peek() == end){
								return;
							}
							sta.pop();
						}
					}
				}
				//sta.pop();

			}
		}
		return;
	}

	//recomputeMST(u,v,weight,G);
	public static int recomputeMST(int u, int v, int weight, EdgeWeightedGraph G, int[][] new_weight_info, int MSTweight){
		int vortexnum = G.V();
		int edgenum = G.E();
		int edgei_n = -1;
		Stack<Integer> sta = new Stack<Integer>();
		Edge newedge = new Edge(u, v, weight);
		int[] visited = new int [vortexnum];
		int newMSTweight = MSTweight;
		int current_weight = 0;
		int temp1, temp2, temp3;

		//Check if the new Edge is too large, if it is, there is no need to renew MST
		if (weight > new_weight_info[vortexnum - 2][2]){
			return newMSTweight;
		}

		//When the new edge is added, the original MST is connected once. In the connected loop, remove the largest edge
		//mst.addEdge(u, v, weight);
		isConnected(new_weight_info, visited, u, v, sta);

		//choose the largest edge

		int edgemax = -1;
		int vmax = -1, wmax = -1;
		int current_v = (int) sta.pop(), previous_v;
		int posi = 0, posi_max = 0;
		/*
		System.out.println("uv");
		System.out.println(u);
		System.out.println(v);
		System.out.println("Stack");
		System.out.println(current_v);
	*/
		while(!sta.isEmpty()){
			previous_v = current_v;
			current_v = (int) sta.pop();
			//System.out.println(current_v);
			for (int i = 0; i < vortexnum - 1; i++){
				if(((previous_v == new_weight_info[i][0]) && (current_v == new_weight_info[i][1])) || ((previous_v == new_weight_info[i][1]) && (current_v == new_weight_info[i][0]))){
					current_weight = new_weight_info[i][2];
					posi = i;
					break;
				}
			}
			if(edgemax < current_weight){
				edgemax = current_weight;
				vmax = previous_v;
				wmax = current_v;
				posi_max = posi;
			}

		}


		//we only need to replace if edgemax > weight
		

		if(edgemax > weight){
			/*
			new_weight_info[posi_max][0] = u;
			new_weight_info[posi_max][1] = v;
			new_weight_info[posi_max][2] = weight;
			*/
			//mst.addEdge(u, v, weight);
			newMSTweight = MSTweight - edgemax + weight;
		
			//resort new_weight_info
			//remove the deleted edge
			
			/*
			System.out.println("Replaced");
			System.out.printf("%d %d %d %n", vmax, wmax, edgemax);
			System.out.println("Replacing");
			System.out.printf("%d %d %d %n", u, v, weight);
			*/


						int i = vortexnum -2;
			for (i = vortexnum-2; i >=0; i--){
				if(((new_weight_info[i][0] == vmax) && (new_weight_info[i][1] == wmax) || (new_weight_info[i][0] == wmax) && (new_weight_info[i][1] == vmax)) && (new_weight_info[i][2] == edgemax)){
					new_weight_info[i][0] = 0;
					new_weight_info[i][1] = 0;
					new_weight_info[i][2] = 0;
					break;
				}
			}

			for (int j = i; j < vortexnum - 2; j++){
				new_weight_info[j][0] = new_weight_info[j+1][0];
				new_weight_info[j][1] = new_weight_info[j+1][1];
				new_weight_info[j][2] = new_weight_info[j+1][2];
			}

			//add the new edge
			for (i = vortexnum-3; i >=0; i--){
				if(new_weight_info[i][2] <= weight){
					break;
				}
			}

			for (int j = vortexnum - 3; j > i; j--){
				new_weight_info[j+1][0] = new_weight_info[j][0];
				new_weight_info[j+1][1] = new_weight_info[j][1];
				new_weight_info[j+1][2] = new_weight_info[j][2];
			}
			
			new_weight_info[i+1][0] = u;
			new_weight_info[i+1][1] = v;
			new_weight_info[i+1][2] = weight;


		}
		return newMSTweight;

	}



	public static EdgeWeightedGraph parseEdges(String graph_file){
		//read edges from files
		int u, v, weight;
		EdgeWeightedGraph G;
		try{
		BufferedReader br = new BufferedReader(new FileReader(graph_file));
		/*try{
			BufferedReader br = new BufferedReader(new FileReader(graph_file));
		}catch(IOException exce){
			exce.printStackTrace();
		}*/
		String line = br.readLine();
		/*try{
 			String line = br.readLine();
		}catch(IOException exce){
  			exce.printStackTrace();
		}*/
		
		String[] split = line.split(" ");
		int vortexnum = Integer.parseInt(split[0]);
		int edgenum = Integer.parseInt(split[1]);
		//System.out.println(edgenum);
		G = new EdgeWeightedGraph(vortexnum);

		for(int i = 0; i < edgenum; i++) {
			line = br.readLine();
			split = line.split(" ");
			u = Integer.parseInt(split[0]);
			v = Integer.parseInt(split[1]);
			weight = Integer.parseInt(split[2]);
			Edge e = new Edge(u, v, weight);
			G.addEdge(e);
			
		}


		//br.close();
 			br.close();
 			return G;
		}catch(IOException exce){
  			exce.printStackTrace();
  			throw new RuntimeException(exce);
		}
		
	}



	public static void main(String[] args){

		String graph_file_begin = "/Users/yixingli/Dropbox/CSE_6140/homework/HW1/MST/data/rmat";
		String graph_file_end = ".gr";
		String change_file_begin = "/Users/yixingli/Dropbox/CSE_6140/homework/HW1/MST/data/rmat";
		String change_file_end = ".extra";
		String output_file_begin = "/Users/yixingli/Dropbox/CSE_6140/homework/HW1/MST/results/rmat";
		String output_file_end = ".out";
		String graph_file, change_file, output_file;
		String[] file_num_string = new String[13];

		file_num_string[0] = "0406";
		file_num_string[1] = "0507";
		file_num_string[2] = "0608";
		file_num_string[3] = "0709";
		file_num_string[4] = "0810";
		file_num_string[5] = "0911";
		file_num_string[6] = "1012";
		file_num_string[7] = "1113";
		file_num_string[8] = "1214";
		file_num_string[9] = "1315";
		file_num_string[10] = "1416";
		file_num_string[11] = "1517";
		file_num_string[12] = "1618";

		for(int file_index = 0; file_index < 13; file_index++){
			System.out.println(file_index);
		graph_file = graph_file_begin + file_num_string[file_index] + graph_file_end;
		change_file = change_file_begin + file_num_string[file_index] + change_file_end;
		output_file = output_file_begin + file_num_string[file_index] + output_file_end;
		
/*
		if (args.length < 3) {
			System.err.println("Unexpected number of command line arguments");
			System.exit(1);
		}
*/

		/*
		String graph_file = args[0];
		String change_file = args[1];
		String output_file = args[2];
		*/
		
		try{
		PrintWriter output;
		output = new PrintWriter(output_file, "UTF-8");

		EdgeWeightedGraph G = parseEdges(graph_file);
		int[][] mst = new int[G.V()-1][3];
		int edgenum = G.E();
		int vortexnum = G.V();
		int[][] new_weight_info = new int [vortexnum - 1][3];
		long startMST = System.nanoTime();
		int MSTweight = computeMST(G, new_weight_info);
		long finishMST = System.nanoTime();

		//Subtract the start time from the finish time to get the actual algorithm 
		//running time
		double MSTtotal = (finishMST - startMST)/1000000.0;

		//Write to output file the initial MST weight and time
		output.println(Integer.toString(MSTweight) + " " + Double.toString(MSTtotal));
		System.out.println(Integer.toString(MSTweight) + " " + Double.toString(MSTtotal));


		//Iterate through changes file
		BufferedReader br = new BufferedReader(new FileReader(change_file));
		String line = br.readLine();
		String[] split = line.split(" ");
		int num_changes = Integer.parseInt(split[0]);
		int u, v, weight;

		while ((line = br.readLine()) != null) {
			split = line.split(" ");
			u = Integer.parseInt(split[0]);
			v = Integer.parseInt(split[1]);
			weight = Integer.parseInt(split[2]);

			//Run your recomputeMST function to recalculate the new weight of the MST 
			//given the addition of this new edge
			//Note: you are responsible for maintaining the MST in order to update 
			//the cost without recalculating the entire MST
			long start_newMST = System.nanoTime();
			int newMST_weight = recomputeMST(u, v, weight, G, new_weight_info, MSTweight);
			MSTweight = newMST_weight;
			long finish_newMST = System.nanoTime();

			double newMST_total = (finish_newMST - start_newMST)/1000000.0;

			//Write new MST weight and time to output file
			output.println(Integer.toString(newMST_weight) + " " + Double.toString(newMST_total));
			//System.out.println(Integer.toString(newMST_weight) + " " + Double.toString(newMST_total));

		}

		output.close();
		br.close();

		}catch(IOException exce){
  			exce.printStackTrace();
  			throw new RuntimeException(exce);
		}
	}
}

}
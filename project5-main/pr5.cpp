#include <climits>
#include <algorithm>
#include <queue>
#include <sys/time.h>
#include "ECLgraph.h"
#include <iostream>

static void Dijkstra(const int src, const ECLgraph& g, int* const dist, int* const parent)
{
  // initialize dist and parent array
  for (int i = 0; i < g.nodes; i++)
  {
    dist[i] = INT_MAX;
    parent[i] = 0;
  }
  dist[src] = 0;

  // set up priority queue with just source node in it
  std::priority_queue< std::pair<int, int> > pq;
  pq.push(std::make_pair(0, src));
  while (pq.size() > 0) {
    // process closest vertex
    const int v = pq.top().second;
    pq.pop();
    const int dv = dist[v];
    // visit outgoing neighbors
    for (int i = g.nindex[v]; i < g.nindex[(v + 1)]; i++) {
      const int n = g.nlist[i];
      const int d = dv + g.eweight[i];
      // check if new lower distance found
      if (d < dist[n]) {
        dist[n] = d;
        parent[n] = v;
        pq.push(std::make_pair(-d, n));

      }
    }
  }

}

/*
			ALGORITHM: BELLMANFORD
				PSUEDOCODE
		
	for i = 0, i < graph nodes , i++
		dist[i] = infinity
		parent[i] = null
	end for loop
	
	dist[src] = 0
	initialize updated
	
	while each vertex whose distance is not infinity
		set updated = false
		
		for vertex = 0, vertex < graph nodes, vertex++
			vdist = dist[vertex]
			
			for i = graph node index[vertex], i < graph node index[vertex + 1], i++
				node = graph node list[i]
				trip = vertex distance + weight of edge
				
				if trip < dist[node]
					distance to node is equal to trip total
					parent[node] = vertex
					update the list, new lower path found
					
			end for loop
		end for loop
	end while loop



				COMMENTS
	In BellmanFord Algorithm, we must assume all node distance totals are infinity.
	Set all distances to infinity except the distance of the source node, which is initialized to zero.
	Set all parent nodes to null. If the distances are infinity we cant assume any parent nodes totals.
	While each vertex, whose distance is not infinity, compute the distance to all adjacent nodes 
	and update their distance if a new lower distance was found. 
	Repeat this process until no distance was updated.
	

*/

static void BellmanFord(const int src, const ECLgraph& g, int* const dist,int* const parent)
{
	
	for( int i = 0; i < g.nodes; i++)//initialize parent array
	{
		dist[i] = INT_MAX;//infinity
		parent[i] = 0;//null
	}
	dist[src] = 0;//initialized to zero
        
	        
	bool updated;//variable to keep track of updated nodes
	while(updated)//while each vertex whose distance is not infinity   
 	{
		updated = false;//initilize to false
  
		for(int vertex = 0; vertex < g.nodes; vertex++)//process all verticies
		{
			const int vdist = dist[vertex];//obtain vertex

			
			for(int i = g.nindex[vertex]; i < g.nindex[(vertex + 1)]; i++)//visit outgoing neighbors
			{
				const int node = g.nlist[i];//obtain each node
				const int trip = vdist + g.eweight[i];//trip equals vertex distance + weight
			
			
				if(trip < dist[node])//check if new lower distnace if found
				{
					dist[node] = trip;//distance to node is equal to trip total
					parent[node] = vertex;
					updated = true;//update the list, new lower path found
        			
      			}
    		}
		}
	}		
}

 
/*
			ALGORITHM: Prim
				PSUEDOCODE
				
	for i = 0, i < graph nodes , i++
		dist[i] = infinity
		parent[i] = null
		prev[i] = null
	end for loop
	
	set up priority queue with just source node in it
	push src node to queue
	
	while(pq.size > 0)
		vertex = pq.top.second
		pq.pop
		color[vertex] = 1
		
		for i = graph node index[vertex], i < graph node index[vertex + 1], i++
			node = graph node list[i]
				trip = vertex distance + weight of edge
				
				if color[n] = 0 & dist[n] > trip
					distance to node is equal to trip total
					prev[node] = vertex
					pq.push pair (dist[n], n)
		end for loop
	end while loop	
	
	
	
				COMMENTS
	Set all distances to infinity, except the distance of the source node, which is initialized to zero.
	Create an empty priority_queue pq. Every item of pq is a pair (weight, vertex). 
	Weight is used as first item of pair as first item is by default used to compare two pairs.
	Insert source vertex into pq and make its key as 0. 
	Then run the following operations in a while loop, while priority queue has elements  
	Process cloest vertex, extract minimum key vertex from pq, and set vertex color to 1
	Nested loop that loops through all the neighbors of vertex and do following 
	If weight of edge is smaller than dist of n and n has not has it color changed            
	Update dist of n, i.e dist[n] = weight(vertex, n) 
	Insert n into the pq and set parent[n] = vertex
*/

static void MST_prim(const int src, const ECLgraph& g, int* const dist, int* const prev, int* const color)
{
 
	for( int i = 0; i < g.nodes; i++)//loop through nodes
	{
		dist[i] = INT_MAX;//Set all distances to infinity
		color[i] = 0;//color = null
		prev[i] = 0;//prev = null
	}
  
  	std::priority_queue< std::pair<int, int> > pq;//set up priority queue with just source node in it
  	pq.push(std::make_pair(0, src));//push src node to queue
  
	while(pq.size() > 0)//loop while priority queue has elements
    {
		const int vertex = pq.top().second;//process cloest vertex
		pq.pop();//pop closest vertex
		
		color[vertex] = 1;//change color of n
	  
		//visit outgoing neighbors of n for loop
		for( int i = g.nindex[vertex]; i < g.nindex[vertex + 1]; i++)
		{
			const int n = g.nlist[i];//nodes in list
			const int trip = g.eweight[n];// weight stored at each node
      
			if(color[n] == 0 && dist[n] > trip)//new lower path
			{
				dist[n] = trip;//distance to node is equal to trip total
				prev[n] = vertex;
				pq.push(std::make_pair(dist[n], n));//update the queue, new lower path found
			}
		}
    }

}


int main(int argc, char* argv[])
{
  printf("MST & SSSP (%s)\n", __FILE__);
  if (argc != 2) {fprintf(stderr, "USAGE: %s input_file_name\n", argv[0]); exit(-1);}

  // process command line
  ECLgraph g = readECLgraph(argv[1]);
  if (g.eweight == NULL) {fprintf(stderr, "ERROR: graph must have weights\n"); exit(-1);}
  printf("input: %s\n", argv[1]);
  printf("nodes: %d\n", g.nodes);
  printf("edges: %d\n", g.edges);
  const int source = 0;
  if ((source < 0) || (source >= g.nodes)) {fprintf(stderr, "ERROR: source_node_number must be between 0 and %d\n", g.nodes); exit(-1);}
  printf("source: %d\n", source);

  //MST Prim
  int* const distance = new int[g.nodes];
  int* const prev = new int[g.nodes];
  int* const color = new int[g.nodes];

  MST_prim(source,g, distance, prev, color);
  //comment this part when running with the large graphs
  for (int v = 1; v < g.nodes; v++) {
    printf("(%d ,  %d) and dist %d \n", v, prev[v],distance[v]);
  }
  timeval start, end;
  // SSSP
  int* const dist = new int[g.nodes];
  int* const parent = new int[g.nodes];
  gettimeofday(&start, NULL);
  BellmanFord(source, g, dist, parent);
  gettimeofday(&end, NULL);
  printf("BellmanFord runtime: %.6f s\n", end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) / 1000000.0);
  int maxnode = 0;
  for (int v = 1; v < g.nodes; v++) {
    if (dist[maxnode] < dist[v]) maxnode = v;
  }
  printf("vertex %d has maximum distance %d\n", maxnode, dist[maxnode]);

  // compare yout solutions
  int* const verify = new int[g.nodes];
  gettimeofday(&start, NULL);
  Dijkstra(source, g, verify, parent);
  gettimeofday(&end, NULL);
  printf("Dijkstra runtime: %.6f s\n", end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) / 1000000.0);
  for (int v = 0; v < g.nodes; v++) {
    if (dist[v] != verify[v]) {fprintf(stderr, "ERROR: verification failed for node %d: %d instead of %d\n", v, dist[v], verify[v]); exit(-1);}
  }
  printf("verification passed\n\n");
 
  // free memory
  delete [] verify;
  delete [] distance;
  delete [] dist;
  delete [] prev;
  delete [] color;
  freeECLgraph(g);
  return 0;
}

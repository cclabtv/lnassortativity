#include "Graph.cpp"
#include <utility>  //pair, make_pair
#include <stdlib.h> //rand()
#include <iostream> //cout
using std::pair, std::make_pair;
using std::cout, std::endl;


/*Generate a classical barabasi albert in linear time in the number of edges. 
Note that the linear time algoritm may generate multiedges or self loops*/
Graph* generate(int n, int d, int maxN, int maxM){
    
    
    srand(time(NULL)); //initialize the pseudorandom generation
    
    Graph* graph = new Graph(n, maxN, maxM, false);
    

    //create a clique of d+1 nodes as a start
    for(int v = 0; v < d+1; v++){
        for(int i = v+1; i < d+1; i++){
            graph->addEdge(v, i);
        }
    }
    
    int *neighbors = new int[d]; 
    int neighbor;

    
    for(int v = d+1; v<n; v++){
        for(int i = 0; i < d; i++){

            int r = rand() % (graph->m); //select a random edge
            int b = rand() % 2; //select a random endpoint of the edge

            if(!b) neighbor = graph->edgeList[r].first;
            else neighbor = graph->edgeList[r].second;

            /*
            //repeat to avoid self-loops
            while(v == neighbor){
                int r = rand() % graph->m; //select a random edge
                int b = rand() % 2; //select a random endpoint of the edge

                if(!b) neighbor = graph->edgeList[r].first;
                else neighbor = graph->edgeList[r].second;
                cout << "selected neighbor: " << neighbor << endl;
            }
            */
            neighbors[i] = neighbor;
        }
        
        //add the selected edges of v in the graph
        for(int i = 0; i < d; i++){
            graph->addEdge(v, neighbors[i]);
        }
        
    }
    delete[] neighbors;
    

    return graph;
}

/*
int main(){


    int n = 100000;
    int d = 4;
    Graph* graph = generate(n, d, n+10, n*d+50);
    
    cout << "Starting computing assortativity sparse" << endl;
    float assortativity = graph->degreeAssortativitySparse();
    cout << "Assortativity: " << assortativity << endl;
*/    
    /*
    cout << "Starting computing assortativity" << endl;
    assortativity = graph->degreeAssortativity();
    cout << "Assortativity: " << assortativity << endl;
    

    return 0;
}*/
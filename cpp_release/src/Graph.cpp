#include <utility>  //pair, make_pair, swap
#include <vector> //vector
#include "utilities.cpp" //variance
#include <iostream> //cout
#include <string> //for std::to_string
#include <ctime> //for std::time
#include <fstream> //for writing to file
#include <queue>
#include <cmath> //for exp(), log()
#include <random>
using namespace std;

/*
Graph class implemented using edgearray of fixed size (maxM)
The edgearray is used for efficient preferential attachment sampling it also uses std::vector
which is a cache friendly datastructure.
It also has the advantage of fast edgeremoval. We simply logically remove the edges by putting
at the end of the array and decreasing its logical size. 
*/
class Graph{
    
    private:
        int maxN;
        int maxM;
        bool* visited;


    public:

        int n;//graph number of nodes
        int m;//graph number of edges


        vector<pair<int, int>> edgeList; //edgelist used for node sampling in the BA
        vector<long> capacities; //capacities of edges
        int* degree; //array that contains in pos i the degree of vertex i 
        long* nodeCapacity; //array containing for each node the total capacity of its incident channels
        long totalCapacity; //the sum of all channel capacities multiplied by 2, or the sum of node capacities
        bool weighted;

        PrefixSumTree *pst; //tree that efficiently store the prefixsum of node capacities
        


        
        Graph(){}
        
        //Graph constructors
        Graph(int maxN, int maxM) {
            n = 0; m=0; totalCapacity=0; 
            this->maxN= maxN;
            this->maxM = maxM;
            
            edgeList = vector<pair<int, int>> (maxM, make_pair(0, 0));   
            capacities = vector<long> (maxM, 0);
            degree = new int[maxN];
            
            pst = new PrefixSumTree(maxN);
            

            nodeCapacity = pst->values;
            visited = new bool[maxN];
            
            for(int i = 0; i < maxN; i++){
                degree[i] = 0;
                visited[i] = false;
            }      
        }

        Graph(int N, int maxN, int maxM) : Graph(maxN, maxM){
            n = N;
            cout << "Graph successfully created!" << endl;
        }


        //Destructor
        ~Graph(){
            
            delete[] degree;
            delete[] visited;
            delete[] pst->tree;
        }
        

        //Add an edge to the graph
        void addEdge(int u, int v){
            
            
            degree[u]++;
            degree[v]++;

            edgeList[m].first = u;
            edgeList[m].second = v;
            
            m++;
        }

        /*add an edge with a specified capacity (weight) */
        void addEdgeWithCapacity(int u, int v, long capacity){

            degree[u]++;
            degree[v]++;
            //printf("%d, %d, %d, %d\n", n, m, u, v);
            edgeList[m].first = u;
            edgeList[m].second = v;


            capacities[m] = capacity;
            //cout << "capacity of the edge:" << capacity << endl;
            //cout << "Updating capacity of node " << u << " :" << nodeCapacity[u] << endl;
            pst->update(u, nodeCapacity[u]+capacity);
            //cout << "New capacity of node " << u << " = " << nodeCapacity[u] << endl;
            //cout << "Updating capacity of node " << v << " :" << nodeCapacity[v] << endl;
            pst->update(v, nodeCapacity[v]+capacity);
            //cout << "New capacity of node " << v << " = " << nodeCapacity[v] << endl;
            totalCapacity += 2*capacity;

            m++;
        }

        //remove single edge from the graph
        void removeEdge(int u, int v){
            
            //find element in the edge list
            for(int i = 0; i < m; i++){
                pair<int, int> e = edgeList[i];
                if((e.first == u && e.second == v) || (e.first == v && e.second == u)){
                    removeEdgebyRef(i);
                    return;
                }
            }
        }

        //Linear removal of edges assuming edges are ordered with increasing indexes
        void removeEdges(int *edges, int k){
            for(int j = k-1; j>=0; j--){
                removeEdgebyRef(edges[j]);
            }
        }

        //remove edge in position e in the edge list
        void removeEdgebyRef(int e){
            //update node capacities
            int u = edgeList[e].first;
            int v = edgeList[e].second;
            
            degree[u]--;
            degree[v]--;

            if(weighted){
                totalCapacity -= 2*capacities[e];
                pst->update(u, std::max( (long)0, nodeCapacity[u]-capacities[e]));
                pst->update(v, std::max((long)0, nodeCapacity[v]-capacities[e]));
            }

            //logical removal by putting it at the end and decreasing m
            swap(edgeList[e], edgeList[m-1]);
            swap(capacities[e], capacities[m-1]);
            m--;
        }
            

        int removeNode(int u){
            int t = degree[u];
            if(t == 0){
                return -1;
            }
            int *edges = new int[t];
            int j = 0;

            for(int i = 0; i < m; i++){
                if(degree[u] == 0){
                    break;
                }
                pair<int, int> e = edgeList[i];
                if(e.first == u || e.second == u){
                    edges[j] = i;
                    j++;
                }
            }
            removeEdges(edges, j);
            delete[] edges;
            return 0;
        }
            
        
        

        //Compute the degree assortativity of the graph using full sized matrix
        float degreeAssortativity(){

            //cout << "Computing degree mixing matrix" << endl;
            pair<int, float**> y = degreeMixingMatrix();
            int dmm_size = y.first;
            float** dmm = y.second;
            
            //cout << "Somma della matrice: " << sumMatrix(dmm, dmm_size) << endl;
            float r = 0.0;

            //printSquareMatrix(dmm, dmm_size);
            
            float *q = new float[dmm_size]; //distribution of the remaining degree
            
            
            //compute the marginal distribution from the joint probability
            for(int i = 0; i < dmm_size; i++){
                q[i] = 0;
                for(int j = 0; j<dmm_size; j++){
                    q[i] = q[i] + dmm[i][j];
                }
                //cout << q[i] << endl; 
            }

            //cout << "somma di q: " << sum(q, dmm_size) << endl;

            //printArray(q, dmm_size);
            for(int i = 0; i < dmm_size; i++){
                for(int j = 0; j<dmm_size; j++){
                    //cout << "r = " << r << "i = " << i << "j = " << j <<   "dm[i][j] = " << dmm[i][j] << "\n\n";     
                    r += (i * j * (dmm[i][j] - (q[i] * q[j])));
                }
            }
            float var = variance(q, dmm_size);
            //cout << "variance: "<< var << endl; 
            
            if(var == 0){
                r = 1;
            }
            else{
                r = r / var;
            }

            //deallocate
            delete[] q;

            for(int i = 0; i < dmm_size; ++i) {
                delete[] dmm[i];
            }
            delete[] dmm;
            if(r < -1){
                cout << r <<endl;
                r = -1;
            }

            

            return r;

        }

        /*return the degree mixing matrix full sized with its size: i.e. in position (i, j) 
        we have the fraction of edges with an endpoint of degree i and the other of degree j*/
        pair<int, float**> degreeMixingMatrix(){
            
            
            //compute the maximum degree in the graph
            int max_degree = 1;
            for(int i = 0; i < n; i++){
                if(degree[i] > max_degree){
                    max_degree = degree[i];
                }
            }

            
            //initialize the mixing matrix
            int size = max_degree;

            
            
            float **matrix = new float*[size];
            for(int i = 0; i < size; ++i) {
                matrix[i] = new float[size];
            }
            

            for(int i = 0; i < size; i++) {
                for(int j = 0; j < size; j++) {
                    matrix[i][j] = 0.0;
                }
            }
        
            //counting the number of edges with an endpoint of remaining degree i and the other of remaining degree j
            for(int i = 0; i < m; i++) {
                //cout << "edge: " << edgeList[i].first << " " << edgeList[i].second << endl;
                int deg_u = degree[edgeList[i].first];
                int deg_v = degree[edgeList[i].second];
                //printf("Max degree: %d\n", size);
                //printf("Degreees: %d, %d\n", deg_u, deg_v);
                matrix[deg_u-1][deg_v-1]++;
                matrix[deg_v-1][deg_u-1]++;
            }
            //printSquareMatrix(matrix, size);

            //normalizing over the number of edges
            for(int i = 0; i < size; i++) {
                for(int j = 0; j<size; j++){
                    matrix[i][j] = matrix[i][j] / (2*m);
                }
                
            }

            return make_pair(size, matrix);
        }

        //Compute the degree assortativity of the graph using sparse matrix
        float degreeAssortativitySparse(){

            //cout << "Computing degree mixing matrix" << endl;
            pair<int, SparseMatrix> y = degreeMixingMatrixSparse();
            int dmm_size = y.first;
            SparseMatrix dmm = y.second;
            
            //cout << "Somma della matrice: " << sumMatrix(dmm, dmm_size) << endl;
            float r = 0.0;

            
            float *q = new float[dmm_size]; //distribution of the remaining degree
            
            //compute the marginal distribution from the joint probability
            for(int i = 0; i < dmm_size; i++){
                q[i] = 0;
                for(int j = 0; j<dmm_size; j++){
                    q[i] = q[i] + dmm.get(i, j);
                }
                //cout << q[i] << endl; 
            }

            //cout << "somma di q: " << sum(q, dmm_size) << endl;

            //printArray(q, dmm_size);
            for(int i = 0; i < dmm_size; i++){
                for(int j = 0; j<dmm_size; j++){
                    //cout << "r = " << r << "i = " << i << "j = " << j <<   "dm[i][j] = " << dmm[i][j] << "\n\n";     
                    r += (i * j * (dmm.get(i, j) - (q[i] * q[j])));
                }
            }
            float var = variance(q, dmm_size);
            //cout << "variance: "<< var << endl; 
            
            if(var == 0){
                r = 1;
            }
            else{
                r = r / var;
            }

            //deallocate
            delete[] q;

            if(r < -1){
                cout << r <<endl;
                r = -1;
            }

            

            return r;

        }

        /*return the degree mixing spaarse matrix (hashmap dictionary) with its size: i.e. in position (i, j) 
        we have the fraction of edges with an endpoint of degree i and the other of degree j*/
        pair<int, SparseMatrix> degreeMixingMatrixSparse(){
            
            
            //creating the sparse matrix
            SparseMatrix *matrix = new SparseMatrix();

            //compute the maximum degree in the graph
            int max_degree = 0;
            for(int i = 0; i < n; i++){
                if(degree[i] > max_degree){
                    max_degree = degree[i];
                }
            }

            
            //initialize the mixing matrix
            int size = max_degree;

        
            //counting the number of edges with an endpoint of remaining degree i and the other of remaining degree j
            for(int i = 0; i < m; i++) {
                
                int deg_u = degree[edgeList[i].first];
                int deg_v = degree[edgeList[i].second];
                int count = matrix->get(deg_u-1, deg_v-1);
                matrix->insert(deg_u-1, deg_v-1, count+1);
                matrix->insert(deg_v-1, deg_u-1, count+1);
            }
            //printSquareMatrix(matrix, size);


            //normalizing over the number of edges


            for(auto& value : matrix->data) {
                value.second = value.second / (2*m);
            }
            

            return make_pair(size, *matrix);
        }


    //@TODO implement computation of the diameter
    int diameter(){
        vector<vector<int>> adjList;
        adjList.resize(n);
        for(int i = 0; i < n; i++){
            adjList[i].resize(degree[i]+1);
        }
        //initialize the adjacency list
        for(int i = 0; i < m; i++){
            adjList[edgeList[i].first].push_back(edgeList[i].second);
            adjList[edgeList[i].second].push_back(edgeList[i].first);
        }
        //do the n bfss
        int max = 0;
        int dist;
        for(int i = 0; i < n; i ++){
            dist = bfs(adjList, i);
            if(dist>max){
                max = dist;
            }
        }

        return max;
    } 

    /*return a pair that represent the interval in which the diameter lies i.e. 
    the longest bfs and 2times the shortest bfs*/
    pair<int, int> diameterApx(){
        vector<vector<int>> adjList;
        adjList.resize(n);
        for(int i = 0; i < n; i++){
            adjList[i].resize(degree[i]+1);
        }
        //initialize the adjacency list
        for(int i = 0; i < m; i++){
            adjList[edgeList[i].first].push_back(edgeList[i].second);
            adjList[edgeList[i].second].push_back(edgeList[i].first);
        }

        int min = n;
        int max = 0;
        int dist, u;
        //do log(n) bfs
        for(int i = 0; i < log2(n)+10; i ++){
            u = rand() % n;
            dist = bfs(adjList, u);
            if(dist < min) min = dist;
            if(dist > max) max = dist;
        }
        return make_pair(max, 2*min);


    }
    int bfs(vector<vector<int>> adjList, int node){
        
        int max = 0;
        queue<pair<int, int>> q;
        q.push(make_pair(node, 0));
        visited[node] = true;
        int v, dist;
        pair<int, int> x;

        while(!q.empty()){
            x = q.front();
            v = x.first;
            dist = x.second;
            if(dist > max){
                max = dist;
            }
            q.pop();
            for(auto & node : adjList[v]){
                if(!visited[node]){
                    q.push(make_pair(node, dist+1));
                    visited[node] = true;
                }
            }
        }

        //reset the visited vector
        for(int i = 0; i < n; i++){
            visited[i] = false;
        }

        return dist;
    }

    //save the graph as a file where each line is an edge separated by a comma
    void saveToFile(std::string filename){
        std::ofstream output_file(filename + ".txt");

        for(int i = 0; i<m; i++){
            output_file << edgeList[i].first << "," << edgeList[i].second << "\n";
        }

        output_file.close();
    }
};

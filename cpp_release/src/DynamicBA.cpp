#include "Graph.cpp"
#include <utility>  //pair, make_pair
#include <stdlib.h> //rand()
#include <iostream> //cout
//#include <random> //std::random_device, std::mt19937, std::uniform_real_distribution<> 
using std::string, std::to_string, std::vector, std::pair;





class DynamicBa{

public:
    
    
    int d;//parameter of the preferential attachment

    Graph *graph;
   
    SparseMatrix *deathRateDic; //dictates how a given edge disappear with the probability stored in deathRateDic.get(u, v)

    std::mt19937_64::result_type seed; // Field to hold the seed
    std::mt19937_64 gen; // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis{0, 1}; // Uniform distribution between 0.0 and 1.0
    
    vector<vector<pair<int, int>>> newEdges;
    vector<vector<pair<int, int>>> deletedEdges;
    int round;

    float beta; //hyperparameter of the probability of edge/node disappearence

    enum class Sel {UNIFORM, DEG_DIST, EXP, EXP_MAXDEG, CAPACITY}; //criteria on how to choose the probability on the edges
    
    //constructors
    DynamicBa(){}

    DynamicBa(int d, int maxN, int maxM, float beta):
        d(d), round(0), beta(beta){
        
        
        //generate the random seed
        seed = std::random_device{}();
        gen.seed(seed);

        
        newEdges.resize(maxN);
        deletedEdges.resize(maxN);

        graph = new Graph(maxN, maxM);
    

        //create a clique of d+1 nodes as a start
        for(int v = 0; v < d+1; v++){
            for(int i = v+1; i < d+1; i++){
                graph->addEdge(v, i);
            }
        }


        
        for(int i = 0; i < graph->m; i++){
            newEdges[round].push_back(graph->edgeList[i]);
            deletedEdges[round].push_back(graph->edgeList[i]);
        }

        
        deathRateDic = new SparseMatrix();
    }





    //add new edges to the node v following the pref. attachment rule
    void prefAttach(int v){
        //cout << "prefAttach" << endl;
        int neighbor;
        int *neighbors = new int[d];

        //random generators
        std::uniform_int_distribution randNode(0, graph->n);
        std::uniform_int_distribution randEdge(0, graph->m);
        std::uniform_int_distribution bi(0, 1);

        for(int i = 0; i < d; i++){
            if(graph->m == 0){ //if the graph has depleted all of its edges, choose d random neighbours
                neighbor = randNode(gen);
                while(v == neighbor){
                    neighbor = randNode(gen);
                }
            }
            else{
                int r = randEdge(gen); //select a random edge
                int b = bi(gen); //select a random endpoint of the edge

                if(!b) neighbor = graph->edgeList[r].first;
                else neighbor = graph->edgeList[r].second;

                
                //repeat to avoid self-loops
                while(v == neighbor){
                    int r = randEdge(gen); //select a random edge
                    int b = bi(gen); //select a random endpoint of the edge

                    if(!b) neighbor = graph->edgeList[r].first;
                    else neighbor = graph->edgeList[r].second;
                    //cout << "selected neighbor: " << neighbor << endl;
                }
            }
            neighbors[i] = neighbor;
        }
        
        //add the selected edges of v in the graph
        for(int i = 0; i < d; i++){

            graph->addEdge(v, neighbors[i]);

            //add the edges to the list of edges added anewEdgest round "round"
            newEdges[round].push_back(make_pair(v, neighbors[i]));
        }
        delete [] neighbors;
    }


    //update death rate dictionary for edges
    void updateDeathRateDic(Sel sel){
   
        float prob;
        if(sel == Sel::UNIFORM){
            prob = beta/graph->n;
            for(int i = 0; i<graph->m; i++){
                deathRateDic->insert(graph->edgeList[i].first,graph->edgeList[i].second, prob);
            }
       }
       
        else if(sel == Sel::DEG_DIST){
            for(int i = 0; i<graph->m; i++){
                float min, max;
                if(graph->degree[graph->edgeList[i].first] > graph->degree[graph->edgeList[i].second]){
                    min = (float) graph->degree[graph->edgeList[i].second];
                    max = (float) graph->degree[graph->edgeList[i].first];
                }else{
                    min = (float) graph->degree[graph->edgeList[i].first];
                    max = (float) graph->degree[graph->edgeList[i].second];
                }
                prob =  1.0/(((max-min)+1)*min);
                deathRateDic->insert(graph->edgeList[i].first,graph->edgeList[i].second, prob);
            }
        }    
        //this use a probability inversely proportional to the differnece of the degrees
        else if(sel == Sel::EXP){
            int u, v;
            for(int i = 0; i<graph->m; i++){
                u = graph->edgeList[i].first;
                v = graph->edgeList[i].second;
                prob = exp((-1)*beta*abs(graph->degree[u]-graph->degree[v]));
                deathRateDic->insert(u, v, prob);
            }
        }
        //this use a probability inversely proportional to the maximum degree 
        else if(sel == Sel::EXP_MAXDEG){
            float prob;
            int u, v;
            for(int i = 0; i<graph->m; i++){
                u = graph->edgeList[i].first;
                v = graph->edgeList[i].second;
                int m = std::max(graph->degree[u], graph->degree[v]);
                prob = exp((-1)*beta*m);
                deathRateDic->insert(u, v, prob);
            }
        }
        else if(sel == Sel::CAPACITY){
            float prob, t;
            int u, v;
            for(int i = 0; i<graph->m; i++){
                u = graph->edgeList[i].first;
                v = graph->edgeList[i].second;
                prob = exp((-1)*beta * graph->totalCapacity/graph->capacities[i]);
                deathRateDic->insert(u, v, prob);
            }
        }
    }
    
    //remove edges according to the deathratedic
    void removeEdges(){
        
        int count = 0;
        int *removedEdges = new int[graph->m]; //save the position of the edges to remove

        for(int i = 0; i<graph->m; i++){
            if(deathRateDic->get(graph->edgeList[i].first, graph->edgeList[i].second) > dis(gen)){
                removedEdges[count] = i;
                //add the edges to the list of edges deleted at round "round"
                deletedEdges[round].push_back(graph->edgeList[i]);
                count++;
            }
        }

        graph->removeEdges(removedEdges, count);


    }

    //Add new d edges to nodes 'node' in graph 'graph' according to the BA model
    void updateNodes(){
        
        
        prefAttach(graph->n);


        graph->n++;

    }
        
        
      

    

    void saveDegDistToFile(string str){
        std::ofstream output_file("../data/degrees_" + str + ".txt");

        for(int i = 0; i<graph->n; i++){
            if(graph->degree[i] > 0){
                output_file << graph->degree[i] << "\n";
            }
        }

        output_file.close();

    }


    void writeHistory(){
         std::ofstream output_file("../data/graph_history_" +  to_string(seed) + ".csv");
         output_file << "newEdges,removedEdges\n";
         string str = "[";
         for(auto & p : newEdges[0]){
            str = str + "(" + to_string(p.first) + "," + to_string(p.second) + "),";
         }
         if(str.back() == ',') str.pop_back();
         str = str + "],[]";
         output_file << str << "\n";

         for(int i = 1; i < round; i++){
            str = "[";
            for(auto & p : newEdges[i]){
                str = str + "(" + to_string(p.first) + "," + to_string(p.second) + "),";
            }
            if(str.back() == ',') str.pop_back();
            str = str + "],[";
            for(auto & p : deletedEdges[i]){
                str = str + "(" + to_string(p.first) + "," + to_string(p.second) + "),";
            }
            if(str.back() == ',') str.pop_back();
            str = str + "]";

            output_file << str << "\n";
         }
          
    }

    //Destructor
    ~DynamicBa(){
        delete graph;
    }

};


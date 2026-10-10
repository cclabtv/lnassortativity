#include "BarabasiAlbert.cpp"
//#include <random> //std::random_device, std::mt19937, std::uniform_real_distribution<> 
using std::string, std::to_string, std::vector, std::pair;





class DynamicBa{

public:
    
    
    int d;//parameter of the preferential attachment

    Graph *graph;
    
    int* rateActVec; //a rate activation vector for nodes: a node activate i.e. add new edges using preferential attachment
   
    SparseMatrix *deathRateDic; //dictates how a given edge disappear with the probability stored in deathRateDic.get(u, v)
    
    float nodeDeathRate; //probability that tells the probability of a node to disappear at each round
    float newNodeRate; //rate of new nodes entering the network at each round that is 1/p in expectaiton

    std::mt19937_64::result_type seed; // Field to hold the seed
    std::mt19937_64 gen; // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis{0, 1}; // Uniform distribution between 0.0 and 1.0
    std::lognormal_distribution<double> logdist; //lognormal distribution initalized in init
    //vectors to keep track of each new or removed edges at each round
    vector<vector<pair<int, int>>> newEdges;
    vector<vector<pair<int, int>>> deletedEdges;
    int round;

    float beta; //hyperparameter of the probability of edge/node disappearence

    int actualNodeNumber;
    
    //criteria for generating random capacities, none indicates unweighted graphs
    enum class Prob {POW_LAW, LOGNORMAL, GPD, NONE}; 
    //GPD is the generalized pareto distribution 

    Prob capacityDist;
    float *distParams; //parameters of the probability distribution

    enum class Sel {UNIFORM, DEG_DIST, EXP, EXP_MAXDEG, CAPACITY}; //criteria on how to choose the probability on the edges
    
    //constructors
    DynamicBa(){}

    DynamicBa(int n, int d, int maxN, int maxM, Prob p, float *pparams, float beta):
        d(d), round(0), capacityDist(p), distParams(pparams), actualNodeNumber(d+1), beta(beta){
        
        nodeDeathRate = -1.0, newNodeRate = -1.0; //set them to -1.0 if not initialized
        
        //generate the random seed
        seed = std::random_device{}();
        gen.seed(seed);

        //set the parameters of the lognormal distribution and initialize the graph
        if(p == Prob::LOGNORMAL){
            logdist.param(std::lognormal_distribution<double>::param_type(pparams[0], pparams[1]));
            graph = new Graph(d+1, maxN, maxM, true);
    

            //create a clique of d+1 nodes as a start
            for(int v = 0; v < d+1; v++){
                for(int i = v+1; i < d+1; i++){
                    long c = 1+log(logdist(gen));
                    graph->addEdgeWithCapacity(v, i, c);
                }
            }
        }
        if(p == Prob::POW_LAW){

            graph = new Graph(d+1, maxN, maxM, true);
            
            //create a clique of d+1 nodes as a start
            for(int v = 0; v < d+1; v++){
                for(int i = v+1; i < d+1; i++){
                    //generate random uniform variable between 0 and 1
                    double u = dis(gen);

                    float alpha = pparams[0];
                    float xmin = pparams[1];

                    //generate powerlaw from uniform variable
                    long c = (long) xmin * std::pow(1 - u, 1 / (1 - alpha));
                    graph->addEdgeWithCapacity(v, i, c);
                }
            }
        }

        if(p == Prob::GPD){
            //gpd = new generalized_pareto_distribution(pparams[0], pparams[1], pparams[2]);
            graph = new Graph(d+1, maxN, maxM, true);
            long c;
            double u = dis(gen);
            float mu = pparams[0];
            float sigma = pparams[1];
            float xi = pparams[2];


             //create a clique of d+1 nodes as a start
            for(int v = 0; v < d+1; v++){
                for(int i = v+1; i < d+1; i++){
                    //generate a gpd distributed variable
                    if (xi == 0) {
                        c = (long) mu - sigma * std::log(u);
                    } else {
                        c = (long) mu + ((sigma / xi) * (std::pow(1 - u, -xi) - 1));
                }
                    graph->addEdgeWithCapacity(v, i, c);

                }
            }
        }

        if(p == Prob::NONE){
            graph = generate(n, d, maxN, maxM);
            actualNodeNumber = n;
        }
        
        newEdges.resize(maxN);
        deletedEdges.resize(maxN);
        
        for(int i = 0; i < graph->m; i++){
            newEdges[round].push_back(graph->edgeList[i]);
            deletedEdges[round].push_back(graph->edgeList[i]);
        }

        


        rateActVec = new int[maxN];
        deathRateDic = new SparseMatrix();
    }

    DynamicBa(int n, int d, int maxN, int maxM, Prob p, float *pparams, float beta, float newNodeRate, float nodeDeathRate)
        : DynamicBa(n, d, maxN, maxM, p, pparams, beta){
            this->nodeDeathRate = nodeDeathRate;
            this->newNodeRate = newNodeRate;
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

    void weightedPrefAttach(int v){
        /*Nodes are sampled by generating a random number r between 0 and the sum of all the node capacities.
        Then the selected node will be the smallest i such that r>i in the prefix sum array.
        This can be done through binary search O(log(n)).
        */
        int neighbor;
        std::uniform_int_distribution<int> randNode(0, graph->n);
        std::uniform_int_distribution<long> randVal(0, graph->totalCapacity);
        long c = 0;
        
        /*instead of selecting d neighbours we select a random geometric r.v. with expectation d*/
        float p = 1/(float)d;
        
        std::exponential_distribution ex((-1)*log(1-p));
        int x = ex(gen);
        int *neighbors = new int[x];

        for(int i = 0; i < x; i++){
            if(graph->m == 0){ //if the graph has depleted all of its edges, choose d random neighbours
                neighbor = randNode(gen);
                while(v == neighbor){
                    neighbor = randNode(gen);
                }
            }
            else{
                
                long r = randVal(gen); //generate a random value
                
                neighbor = graph->pst->search(r); //search the bin in the prefix sum tree corresponding to that value

                
                //repeat to avoid self-loops
                while(v == neighbor){
                    long r = randVal(gen); //generate a random value

                    neighbor = graph->pst->search(r); //search the bin in the prefix sum tree corresponding to that value
                }
            }
            neighbors[i] = neighbor;
        }

        
        //add the selected edges of v in the graph
        for(int i = 0; i < x; i++){
            
            if(capacityDist == Prob::LOGNORMAL){
                c = 1+log(logdist(gen));
                //cout << c << endl;
            }

            if(capacityDist == Prob::POW_LAW){
                //generate random uniform variable between 0 and 1
                double u = dis(gen);

                float alpha = distParams[0];
                float xmin = distParams[1];

                //generate powerlaw from uniform variable
                c = (long) xmin * std::pow(1 - u, 1 / (1 - alpha));

                
            }

            if(capacityDist == Prob::GPD){
               //generate random uniform variable between 0 and 1
                double u = dis(gen);

                float mu = distParams[0];
                float sigma = distParams[1];
                float xi = distParams[2];
                //generate gpd from uniform variable
                if (xi == 0) {
                        c = (long) mu - sigma * std::log(u);
                    } else {
                        c = (long) mu + ((sigma / xi) * (std::pow(1 - u, -xi) - 1));
                }
            }

            graph->addEdgeWithCapacity(v, neighbors[i], c);
            //add the edges to the list of edges added at round "round"
            newEdges.at(round).push_back(make_pair(v, neighbors[i]));
        }
        delete [] neighbors;
    }



    //update the activation rate vector for nodes
    void updateRateActVec(){
        //cout << "updateRateActVec" << endl;
        for(int i = 0; i < graph->n; i++){
            if(graph->m!=0){
                rateActVec[i] = 50/graph->m;
            }
        }
    }

    //update death rate dictionary for edges
    void updateDeathRateDic(Sel sel){
    //cout << "updateDeathRateDic" << endl;
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
        //cout << "removeEdges" << endl;
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
        
        //case of the standard unweighted preferential attachment
        if(capacityDist == Prob::NONE){
            /*
            for(int i = 0; i < graph->n; i++){
                if(rateActVec[i] > dis(gen)){
                    prefAttach(i);
                }
            }
            */
            prefAttach(graph->n);
            
        }

        //weighted case
        else{
            weightedPrefAttach(graph->n);
            
        }

        graph->n++;
        actualNodeNumber++;

        if(newNodeRate > 0.0){
            //using an exponential random variable for generating a geometric rv
            //this can be done by discretizing an exponential with parameter -log(1-p)
            

            std::exponential_distribution ex((-1)*log(1-newNodeRate));
            int x = ex(gen);
            //cout << newNodeRate << "x = " << x << endl;

            

            for(int i = 1; i<x; i++){
                weightedPrefAttach(graph->n);
                graph->n++;
                actualNodeNumber++;
            }

        }
        
        
        //sample the nodes for removal: each node as a probability nodeDeathRate to disappear
        if(nodeDeathRate > 0.0){
            
            //sample a binomial random variable  with parameter nodeDeathRate
            // when n is large enough instead of using a binomial we can
            // use the poisson approximation  since p is very small
            std::uniform_int_distribution<int> randNode(0, graph->n);

            int k;

            if(graph->n > 50){
                std::poisson_distribution pois(graph->n*nodeDeathRate);
                k = pois(gen);
            }
            else{
                std::binomial_distribution bin(graph->n, nodeDeathRate);
                k = bin(gen);
            }
            //cout << "k = " << k << "\tn = " << actualNodeNumber << "\tm = " << graph->m << endl;

            for(int i = 0; i<k; i++){
                //int u = randNode(gen);
                //cout << "removed node: "<< u << "\tround: " << round << endl;
                
                //continue in case the node is already removed from the network
                while(graph->removeNode(randNode(gen)) == -1){}
                actualNodeNumber--;
            }
            
        }
        

        //for now if we input the probability of node disappearence as 0 we 
        //use that each node disappear with probability inversely proportional to a function 
        // of its degree
        /*
        if(nodeDeathRate == 0.0){
             
            for(int i = 0; i < graph->n; i++){
                if(graph->degree[i] > 0 &&  exp((-1)*beta*graph->degree[i]) > dis(gen)){
                    graph->removeNode(i);
                    actualNodeNumber--;
                }
            }
        }*/
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
         std::ofstream output_file("../data/historycalGraph" +  to_string(std::time(0)) + ".csv");
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
        delete [] rateActVec;
    }

};


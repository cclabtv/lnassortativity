#include "DynamicBA.cpp"

//this use a probability inversely proportional to the maximum degree 

int simulate(int n, int d, int rounds, float beta, int skip, float p, float xi){

   
    DynamicBa::Sel sel = DynamicBa::Sel::CAPACITY;
    //parameters for powerlaw
    //float params[2] = {2.05, 1};
    float q = 0; //setting the probability for a node to disappear equal to 0 for now

    //parameters for gpd
    float mu = 0.0; //location
    float sigma = 10000; //scale
    //xi is the shape of the gpd 
   /* when xi is equal -1 the distribution is a uniform in [0, sigma] 
   while increasing the parameter will make the distribution increasingly similar to a powerlaw
   */
    float params[3] = {mu, sigma, xi};
    int maxN;
    
    if(p != 0){maxN = n + rounds + (1/p)*(rounds)+50;}
    else{maxN = n + rounds +50;}

    int maxM = maxN*2*d+100;
    DynamicBa *dynBa = new DynamicBa(n, d, maxN, maxM, DynamicBa::Prob::GPD, params, beta, p, q);
    std::ofstream file;
    string str = "";
    
    //convert beta and xi to two decimal string
    string sbeta = to_string(beta);
    sbeta = sbeta.substr(0, sbeta.find(".")+3);
    string sxi = to_string(xi);
    sxi = sxi.substr(0, sxi.find(".")+3);

    str = str + "d_" + to_string(d) + "_beta_" + sbeta + "_xi_" + sxi +
    "_rounds_" + to_string(rounds) + "_" + to_string(std::time(0));
    cout << str << endl;

    file.open("../data/assortativityWcap_" + str + ".csv");

    //first line is the seed of the simulation
    file << "#" << dynBa->seed << "\n";
    file << "n,m,assortativity,diameter\n";
    float r;
    pair<int, int> diam;
    for(int i = 0; i < rounds; i++){
       if(beta > 0){
            dynBa->updateDeathRateDic(sel); 
        }
       dynBa->updateNodes();
       if(beta > 0){
            dynBa->removeEdges();   
        }
       
       dynBa->round++;

        if(i % skip == 0){
            r = dynBa->graph->degreeAssortativity();
            diam = make_pair(0, 0);
            diam = dynBa->graph->diameterApx(); //make_pair<int, int>(0, 0);
            cout << "Assortativity: " << r << "\tDiam: [" << diam.first << "," << diam.second << "]" << endl;
            file <<  dynBa->actualNodeNumber << "," << dynBa->graph->m << "," << r << "," <<
                diam.first << " " << diam.second << "\n";
        }
    }
    dynBa->saveDegDistToFile(str);
    dynBa->writeHistory();
    
    return 0;
}


void print_help() {
    std::cout << "Usage: simulate [-h] [n d rounds beta skip p xi]\n"
              << "Options:\n"
              << "\t-h        Show this help message\n"
              << "\td         The number of neighbors in the preferential attachment\n"
              << "\trounds    Number of iteration\n"
              << "\tbeta      Parameter in the probability of edge disappearenc if 0 no edge disappearence\n"
              << "\tskip      Frequency at which assortativity is computed 1 every round, 10 every 10 rounds etc.\n"
              << "\tp         Probability of a node in the network to make new edges\n"
              << "\txi        Shape of the Generalized Pareto Distribution for choosing the edge weights\n";
              //<< "\tq         Probability of a node (and all its incident edges) to disappear from the network\n";
}

int main(int argc, char* argv[]) {
    if (argc == 1 || std::string(argv[1]) == "-h") {
        print_help();
        return 0;
    }
    if (argc != 7) {
        std::cerr << "Invalid number of arguments. Use -h for help.\n";
        return 1;
    }
    int d = std::stoi(argv[1]);
    int n = 2*d;
    int rounds = std::stoi(argv[2]);
    float beta = std::stof(argv[3]);
    int skip = std::stof(argv[4]);
    float p = std::stof(argv[5]);
    float q = std::stof(argv[6]);
    return simulate(n, d, rounds, beta, skip, p, q);
}




#include "DynamicBA.cpp"
#include <string>

//this use a probability inversely proportional to the maximum degree 

int simulate(int n, int d, int rounds, float beta, int skip){

   
    DynamicBa::Sel sel = DynamicBa::Sel::EXP_MAXDEG;
    DynamicBa *dynBa = new DynamicBa(d, n+rounds+50, (n+rounds)*100*d, beta);
    std::ofstream file;
    string str = "";

    //convert beta to two decimal string
    string sbeta = to_string(beta);
    sbeta = sbeta.substr(0, sbeta.find(".")+3);

    str = str + "d_" + to_string(d) + "_beta_" + sbeta + 
    "_rounds_" + to_string(rounds) + "_" + to_string(dynBa->seed);
    cout << str << endl;

    file.open("../data/assortativity_" + str + ".csv");

    //first line is the seed of the simulation
    file << "#" << dynBa->seed << "\n";
    file << "n,m,assortativity,diameter\n";
    float r;
    pair<int, int> diam;
    for(int i = 0; i < rounds; i++){
       
       dynBa->updateDeathRateDic(sel); 
       dynBa->updateNodes();
       dynBa->removeEdges();
        //cout << i << endl;
       dynBa->round++;

        if(i % skip == 0){
            r = dynBa->graph->degreeAssortativitySparse();
            
            //cout << "Assortativity: " << r << endl; //print assortativity on screen
            
            //save to file number of nodes, number of edges and assortativity as a csv
            file <<  dynBa->graph->n << "," << dynBa->graph->m << "," << r << "\n"; 
        }
    }
    dynBa->saveDegDistToFile(str);
    dynBa->writeHistory();
    
    return 0;
}


void print_help() {
    std::cout << "Usage: simulate [-h] [n d rounds beta]\n"
              << "Options:\n"
              << "\t-h        Show this help message\n"
              << "\td         The number of neighbors in the preferential attachment\n"
              << "\trounds    Number of iteration\n"
              << "\tbeta      Parameter in the probability of edge disappearence\n"
              << "\tskip      Frequency at which assortativity is computed 1 every round, 10 every 10 rounds etc.\n";
}

int main(int argc, char* argv[]) {
    if (argc == 1 || std::string(argv[1]) == "-h") {
        print_help();
        return 0;
    }
    if (argc != 5) {
        std::cerr << "Invalid number of arguments. Use -h for help.\n";
        return 1;
    }
    int d = std::stoi(argv[1]);
    int n = 2*d;
    int rounds = std::stoi(argv[2]);
    float beta = std::stof(argv[3]);
    int skip = std::stof(argv[4]);
    return simulate(n, d, rounds, beta, skip);
}




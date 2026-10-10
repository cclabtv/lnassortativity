#include "DynamicBA.cpp"

int simulate(int d, int rounds){
    int n = 2*d;
    DynamicBa *dynBa = new DynamicBa(n, d, 200+rounds, (1000+rounds)*100*d, DynamicBa::Prob::NONE, NULL, 0);
    std::ofstream file;
    

    file.open("../data/assortativityBA_d_" + to_string(d) + "_rounds_" + to_string(rounds) 
        + "_" + to_string(std::time(0)) + ".csv");
    
    file << "n,m,assortativity\n";
    float r;
    //int diam;

    for(int i = 0; i < rounds; i++){
       
        dynBa->prefAttach(dynBa->graph->n);
        dynBa->graph->n++;
        //cout << i << endl;
       
            r = dynBa->graph->degreeAssortativity();
            //diam = dynBa->graph->diameter();
            cout << "Assortativity: " << r << endl;
            file <<  dynBa->graph->n << "," << dynBa->graph->m << "," << r << "\n";
       
    }
    
    return 0;
}

void print_help() {
    std::cout << "Usage: simulate [-h] [d rounds]\n"
              << "Options:\n"
              << "\t-h        Show this help message\n"
              << "\td         The number of neighbors in the preferential attachment\n"
              << "\trounds    Number of iteration\n";
}

int main(int argc, char* argv[]) {
    if (argc == 1 || std::string(argv[1]) == "-h") {
        print_help();
        return 0;
    }
    if (argc != 3) {
        std::cerr << "Invalid number of arguments. Use -h for help.\n";
        return 1;
    }
    int d = std::stoi(argv[1]);
    int rounds = std::stoi(argv[2]);
    return simulate(d, rounds);
}

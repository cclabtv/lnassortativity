#include <iostream> //cout, endl
#include <unordered_map> //unordered_map
#include <boost/functional/hash.hpp> //boost::hash<std::pair<int, int>>
#include <utility> // for std::pair
#include <algorithm> //for std::min, std::max
#include <cmath> //for exp(), log(), log2(), ceil
using std::unordered_map;
using std::cout, std::endl;

// Define the SparseMatrix class
class SparseMatrix {

public:
    // Use an unordered_map with a pair of integers as the key (representing the row and column indices)
    // and a float as the value (representing the value at that position in the matrix).
    // The hash function from the Boost library is used to hash the pair of integers.
    std::unordered_map<std::pair<int, int>, float, boost::hash<std::pair<int, int>>> data;


    // Insert a value into the matrix at position (i, j).
    // Only non-zero values are stored, as zero values are the default in a sparse matrix.
    void insert(int i, int j, float val) {
        if (val != 0) {
            data[{i, j}] = val;
        }
        else{
            data.erase({i, j});
        }
    }

    // Get the value at position (i, j) in the matrix.
    // If the position does not exist in the data map, it means that the value at that position is zero,
    // so zero is returned. Otherwise, the stored value is returned.
    double get(int i, int j) {
        if(data.find({i, j}) != data.end()) {
            return data[{i, j}];
        } else {
            return 0.0;
        }
    }
};

/*Efficient data structure for storing the prefix sum of an array with update time, query time and search time O(log(n)).
It is structured as a binary tree of 2n nodes where each node store the prefix sum of the leaves in its left subtree, and the leaves 
represents the original value in the array. This is used to implement an efficient "weighted" sample, where you extract elements in 
an array according to their internal values*/
class PrefixSumTree{
    public:

        //length of the original vector, the tree will be of length 2n
        int size;

        //number of internal nodes in the tree
        int internalNodes;
        /*
        tree represented as a vector in which i is parent of 2i (left) and 2i+1 (right)
        the root is in position 1 not 0

        the internal node of the tree are 2^(ceiling(log_2(n)) - 1)
        */
        long* tree;

        /*
        values represent the original values for which the prefix sum is computed
        that are the leaves of the prefixsum-tree
        */
        long* values; 
        


        PrefixSumTree(int n) : size(n){

            internalNodes = pow(2, ceil(log2(n)));
            tree = new long[size+internalNodes+1]();
            values = tree + internalNodes;

        }

        /*
        not implemented because not used 
        int prefSum(int i){
            
        }
        */



        //update the tree when the value in position i is updated to the value val
        void update(int i, long val){
            

            int j = i+internalNodes;
            long diff = val - values[i];
            int father;

            tree[j] = val;

            while(j > 1){
                
                father = j / 2;
                if(j%2 == 0){
                    tree[father] += diff;
                    
                }


                j = father;
            }
        }

        /*
        return the position i in the original vector s.t. x is between 
        prefixsum[i-1] and prefixsum[i]
        */
        int search(long x){
            
            int i = 1;
            long t = x;


            while(i < internalNodes){
                
                if(t <= tree[i]){
                    i = 2*i;
                }
                
                else{
                    t=t-tree[i];
                    i = 2*i +1;
                }   
                
            }
            return i-internalNodes;
        }

    void print(){
        for(int i= 1; i < internalNodes+size+1; i++){
            printf("tree[%d] = %ld\n", i, tree[i]);
        }
    }

};

//simple function to compute the variance of a probability vector
float variance(float* prob_vector, int size) {
    float mean = 0.0;
    for(int i = 0; i < size; i++) {
        mean += prob_vector[i]*i;
    }

    float var = 0.0;
    for(int i = 0; i < size; i++) {
        float diff = i - mean;
        var += prob_vector[i]*diff*diff;
    }

    return var;
}

float sum(float *arr, int size){
    float sum = 0.0;
    for(int i = 0; i < size; i++) {
        sum += arr[i];
    }
    return sum;
}

//compute the sum of a square matrix
float sumMatrix(float **matrix, int size){
    float s = 0.0;
    for(int i = 0; i < size; i++) {
        s += sum(matrix[i], size);
    }
    return s;
}


void printSquareMatrix(float **matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void printArray(long *arr, int size){
    for (int i = 0; i < size; i++) {
        cout << arr[i] << " ";
    }
    cout << endl;
}




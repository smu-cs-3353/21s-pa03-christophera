#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include "graphing.h"

using namespace boost;

int main (int argc, char * argv[]) {
    std::string input = argv[1];
    std::string output = argv[2];
    graphing myGraph(input);
    myGraph.calc_clusters(output);
    return 0;
}
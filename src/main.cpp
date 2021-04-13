#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include "graphing.h"

using namespace boost;

int main (int argc, char * argv[]) {
    std::string name = argv[1];
    graphing myGraph(name);
    myGraph.calc_clusters();
    return 0;
}
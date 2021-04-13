#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include "graphing.h"

using namespace boost;

int main(int,char*[])
{
    std::string name = "input.txt";
    graphing myGraph(name);
    myGraph.calc_clusters();
    return 0;
}
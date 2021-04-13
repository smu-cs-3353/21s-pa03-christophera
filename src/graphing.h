//
// Created by Chris on 4/8/2021.
//

#ifndef COMMUNITYDETECTION_GRAPHING_H
#define COMMUNITYDETECTION_GRAPHING_H

#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>

using namespace boost;
struct Edge_with_score{
    double score = 0;
};
struct vertex_with_path{
    int s_paths = 0;
};
class graphing {
private:
    typedef adjacency_list<setS, vecS, undirectedS,vertex_with_path,Edge_with_score> Graph;
    Graph g;
    void calc_edge_betweeness(Graph& g);
    std::pair<int,int> * edge_array;
    int count_vertices;
    void clear_edge_scores(Graph& g);
    void remove_highest_edge_scores(Graph& g);

public:
    graphing();
    graphing(std::pair<int,int> Edges[], int vertices_num);
    explicit graphing(std::string fileName);
    void printVertices();
    void printEdges();
    void calc_clusters();
};
using namespace boost;
#endif //COMMUNITYDETECTION_GRAPHING_H

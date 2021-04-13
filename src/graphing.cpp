//
// Created by Chris on 4/8/2021.
//


#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/bc_clustering.hpp>
#include "graphing.h"

graphing::graphing(std::string fileName) {
    std::string currentLine;
    std::ifstream file ("data/"+fileName);
    if (file.is_open())
    {

        getline (file,currentLine);
        int lines = std::stoi(currentLine);
        edge_array = new std::pair<int,int>[lines];

        std::vector<int> vec;
        int a = 0;
        while ( getline (file,currentLine) )
        {
            int vert1= currentLine[0]-65;
            int vert2= currentLine[4]-65;//ascii values -65 (A==0)
            vec.push_back(vert1);
            vec.push_back(vert2);
            edge_array[a] = std::pair<int,int>(vert1,vert2);
            a++;
        }
        file.close();

        //Unique nodes are found to determine the size of our graph(num of vertices)
        sort(vec.begin(), vec.end());
        std::vector<int>::iterator it;
        it = std::unique(vec.begin(), vec.end());
        vec.resize(distance(vec.begin(),it));
        count_vertices = vec.size();

        Graph temp(count_vertices);
        g = temp;

        // add the edges to the graph object
        for (int i = 0; i < lines; ++i)
            add_edge(edge_array[i].first, edge_array[i].second, g);

    }
    else
    {
        std::cout << "Unable to open file"<<std::endl;
        exit(1);
    }
}

void graphing::printVertices() {
    typedef graph_traits<Graph>::vertex_descriptor Vertex;

    // get the property map for vertex indices
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);

    std::cout << "vertices(g) = ";
    typedef graph_traits<Graph>::vertex_iterator vertex_iter;
    std::pair<vertex_iter, vertex_iter> vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        Vertex v = *vp.first;
        std::cout << index[v] <<  " ";
    }
    std::cout << std::endl;
}

void graphing::printEdges() {
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);
    std::cout << "edges(g) = " << std::endl;
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        std::cout << "(" << index[source(*ei, g)]
                  << "," << index[target(*ei, g)] << ") ";

        Graph::edge_descriptor e = *ei;
        std::cout << "score: " << g[e].score<<std::endl;

    }
    std::cout << std::endl;


}

template < class Graph > struct exercise_vertex
{
    //state vars
    Graph& g;
    const char* name;

    //Constructor
    exercise_vertex(Graph& g_, const char name_[]) : g(g_), name(name_) {}

    //vertex descriptor
    typedef typename graph_traits< Graph >::vertex_descriptor Vertex;
    void operator()(const Vertex& v) const
    {
        typename property_map< Graph, vertex_index_t >::type vertex_id
                = get(vertex_index, g);
        std::cout << "vertex: " << name[get(vertex_id, v)] << std::endl;

        // Write out the outgoing edges
        std::cout << "\tout-edges: ";
        typename graph_traits< Graph >::out_edge_iterator out_i, out_end;
        typename graph_traits< Graph >::edge_descriptor e;
        for (boost::tie(out_i, out_end) = out_edges(v, g); out_i != out_end; ++out_i)
        {
            e = *out_i;
            Vertex src = source(e, g), targ = target(e, g);
            std::cout << "(" << name[get(vertex_id, src)] << ","
                      << name[get(vertex_id, targ)] << ") ";
        }
        std::cout << std::endl;

        // Write out the incoming edges
        std::cout << "\tin-edges: ";
        typename graph_traits< Graph >::in_edge_iterator in_i, in_end;
        for (boost::tie(in_i, in_end) = in_edges(v, g); in_i != in_end; ++in_i)
        {
            e = *in_i;
            Vertex src = source(e, g), targ = target(e, g);
            std::cout << "(" << name[get(vertex_id, src)] << ","
                      << name[get(vertex_id, targ)] << ") ";
        }
        std::cout << std::endl;

        // Write out all adjacent vertices
        std::cout << "\tadjacent vertices: ";
        typename graph_traits< Graph >::adjacency_iterator ai, ai_end;
        for (boost::tie(ai, ai_end) = adjacent_vertices(v, g); ai != ai_end;
             ++ai)
            std::cout << name[get(vertex_id, *ai)] << " ";
        std::cout << std::endl;
    }

};


void graphing::calc_clusters()
{
    int vertices_count = num_vertices(g);
    auto edge_count = num_edges(g);
    double vertex_degrees[vertices_count];
    Graph::vertex_iterator v, vend;
    int i =0;
    for (boost::tie(v, vend) = vertices(g); v != vend; ++v)
    {
        vertex_degrees[i] = degree(*v,g);
        i++;
    }

    Graph g_copy;
    copy_graph(g,g_copy);

    // create A-P matrix to later calculate modularites of components
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);
    //A Matrix, elements contain 1 if the edge between (u,v) exists 0 otherwise
    double A_Matrix[vertices_count][vertices_count] = {{0}};
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        A_Matrix[index[source(*ei, g)]][index[target(*ei, g)]] = 1;
        A_Matrix[index[target(*ei, g)]][index[source(*ei, g)]] = 1;

    }
    //P Matrix, expected # of edges, degree(u) * degree(v)  / 2 * num_edges
    double P_Matrix[vertices_count][vertices_count] = {{0}};

    for (int k = 0; k < vertices_count; k++) {
        for (int j = 0; j < vertices_count; j++) {
            P_Matrix[k][j] = (vertex_degrees[k]*vertex_degrees[j]) / (2*edge_count);
        }
    }
    // A-P Matrix, subtract matrices
    double A_P_Matrix[vertices_count][vertices_count] = {{0}};
    for (int k = 0; k < vertices_count; k++) {
        for (int j = 0; j < vertices_count; j++) {
            A_P_Matrix[k][j] = A_Matrix[k][j] - P_Matrix[k][j];
        }
    }

    using Component = int;
    using Mapping = std::map<Graph::vertex_descriptor, Component>;
    Mapping mappings;
    int n = boost::connected_components(g_copy, boost::make_assoc_property_map(mappings));

    std::map<int,double> modularity_at_each_removal;
    modularity_at_each_removal.insert({0,0});
    int laps = 1;
    while(n!=num_vertices(g_copy))
    {
        clear_edge_scores(g_copy);
        calc_edge_betweeness(g_copy);
        remove_highest_edge_scores(g_copy);
        double modularity =0;
        //finds subcomponents of the graph after deleting edges
        n = boost::connected_components(g_copy, boost::make_assoc_property_map(mappings));
        std::cout<<std::endl<<"current communities:"<<std::endl;
        for (Component c = 0; c<n; ++c) {
            std::vector<int> current_vertices;
            std::cout<<"community "<< c<<": ";
            //place the vertices of the current component into an array
            for (auto& mapping : mappings)
                if (mapping.second == c)
                {
                    current_vertices.push_back(mapping.first);
                    std::cout << " " << mapping.first;
                } std::cout<<std::endl;
            //TODO calculate the modularity between the current vertices,we need to calculate modularites of edges and vertices and see what
            //TODO we get when we add them (and /4)
            //partition AP Matrix and add up those partitions
            for(int a = 0; a <current_vertices.size(); a++)
            {
                for(int b = 0; b<current_vertices.size();b++)
                {
                   modularity+= A_P_Matrix[current_vertices[a]][current_vertices[b]];
                }
            }

        }
        modularity = modularity / (edge_count*2.0);
        //place the modularity into a map
        modularity_at_each_removal.insert({laps, modularity});
        laps++;
        std::cout<<"modularity: " << modularity;std::cout<<std::endl<<std::endl;
    }

    //map our Max modularities to a multimap to sort by highest of all
    std::multimap<double,int> MM;
    typedef std::multimap<double,int>::iterator MMAPIterator;
    for (auto& it : modularity_at_each_removal) {
        MM.insert({ it.second, it.first });
    }
    //the last key value will have the number of removals that produced the most modularity
    std::pair<MMAPIterator, MMAPIterator> result = MM.equal_range(MM.rbegin()->first);
    auto it = result.first;

    int optimal_removals = 0;
    optimal_removals=it->second;

    //remove the highest edge scores X amount of times
    for(int j =0; j<optimal_removals;j++)
    {
        clear_edge_scores(g);
        calc_edge_betweeness(g);
        remove_highest_edge_scores(g);
    }
    std::cout<<"optimal removals: " << optimal_removals<<std::endl;
    //finds subcomponents of the graph after deleting edges
    using Component = int;
    using Mapping = std::map<Graph::vertex_descriptor, Component>;
    int n2 = boost::connected_components(g, boost::make_assoc_property_map(mappings));
    std::ofstream out_file ("data/communities.txt");
    for (Component c = 0; c<n2; ++c) {
        std::cout << "component " << c << ":";
        out_file << "Community " << c << ":";
        for (auto& mapping : mappings)
            if (mapping.second == c)
            {
                std::cout << " " << mapping.first;
                out_file << " " << mapping.first;
            }
        std::cout << "\n";
        out_file << "\n";
    }
    out_file.close();
}



template <class ParentDecorator>
struct print_parent {
    print_parent(const ParentDecorator& p_) : p(p_) { }
    template <class Vertex>
    void operator()(const Vertex& v) const {
        std::cout << "parent[" << v << "] = " <<  p[v]  << std::endl;
    }
    ParentDecorator p;
};


template <class NewGraph, class Tag>
struct graph_copier
        : public boost::base_visitor<graph_copier<NewGraph, Tag> >
{
    typedef Tag event_filter;

    graph_copier(NewGraph& graph) : new_g(graph) { }

    template <class Edge, class Graph>
    void operator()(Edge e, Graph& g) {
        boost::add_edge(boost::source(e, g), boost::target(e, g), new_g);
    }
private:
    NewGraph& new_g;
};

template <class NewGraph, class Tag>
inline graph_copier<NewGraph, Tag>
copy_graph(NewGraph& g, Tag) {
    return graph_copier<NewGraph, Tag>(g);
}

//traverse through the tree recording the # of shortest paths at each vertex
void graphing::calc_edge_betweeness(Graph& g) {

    typedef Graph::vertex_descriptor Vertex;
    typename property_map< Graph, vertex_index_t >::type vertex_id
            = get(vertex_index, g);
    Graph G_copy(num_vertices(g));
    // Array to store predecessor (parent) of each vertex. This will be
    // used as a Decorator (actually, its iterator will be).
    std::vector<Vertex> p(boost::num_vertices(g));
    // VC++ version of std::vector has no ::pointer, so
    // I use ::value_type* instead.
    typedef std::vector<Vertex>::value_type *Piter;

    // Array to store distances from the source to each vertex .  We use
    // a built-in array here just for variety. This will also be used as
    // a Decorator.
    boost::graph_traits<Graph>::vertices_size_type d[num_vertices(g)];
    std::fill_n(d, num_vertices(g), 0);

    typename graph_traits< Graph >::adjacency_iterator ai, ai_end;

    // The source vertex
    typedef graph_traits<Graph>::vertex_iterator vertex_iter;
    std::pair<vertex_iter, vertex_iter> vp;
    std::vector<int> visited_nodes;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        Vertex s = *vp.first;
        p[s] = s;
        boost::breadth_first_search
                (g, s,
                 boost::visitor(boost::make_bfs_visitor
                                        (std::make_pair(boost::record_distances(d, boost::on_tree_edge()),
                                                        std::make_pair
                                                                (boost::record_predecessors(&p[0],
                                                                                            boost::on_tree_edge()),
                                                                 copy_graph(G_copy, boost::on_examine_edge()))))));



        //Place the distances into a node and get the max distance found(cause i dont know how to effectively use vertices_size_type!)
        std::vector<int> distances_temp;
        distances_temp.reserve(num_vertices(g));
        for(int i =0; i <num_vertices(g); i++)
        {
            distances_temp.push_back(d[i]);
        }
        auto max_levels = max_element(std::begin(distances_temp), std::end(distances_temp));

        //create and array of vectors with the size of the amount of levels we'll need [max distance]
        std::vector<int> level_and_nodes[*max_levels + 1];
        //place each vertex on a level depending on their distance from the source vertex
        for(int i =0; i <num_vertices(g); i++)
        {
            level_and_nodes[d[i]].push_back(i);
        }
        /*
        // prints level_and_nodes(debugging)
        for (int i = 0; i < *max_levels+1; i++) {

            std::cout << "Elements at index "
                 << i << ": ";

            // Displaying element at each column,
            // begin() is the starting iterator,
            // end() is the ending iterator
            for (auto it_print = level_and_nodes[i].begin();
                 it_print != level_and_nodes[i].end(); it_print++) {

                // (*it) is used to get the
                // value at iterator is
                // pointing
                std::cout << *it_print << ' ';
            }
            std::cout << std::endl;
        }
        */
        // PLACE NODES IN GROUPS BASED ON THEIR DISTANCE AND EACH GROUP ABOVE GIVES THEIR CURRENT # of S_Paths TO ADJECENT NODES BELOW 1 LEVEL DISTANCE
        // ROOT STARTS WITH 1
        int nodes_and_spath_temp[num_vertices(g)] ={0};
        nodes_and_spath_temp[level_and_nodes[0][0]]++;  //grant source node initial value of 1 to send to level below

        std::pair<vertex_iter, vertex_iter> vertex_it = vertices(g);
        for (int i = 0; i < *max_levels+1; i++)
        {
            for (auto floor1_node = level_and_nodes[i].begin(); floor1_node != level_and_nodes[i].end(); floor1_node++)
            {
                //transfer flow down 1 level, dont do it for last level
                if(i != *max_levels)
                {
                    std::vector<int> adj_nodes;
                    //get adjacent nodes to the current node
                    Vertex v = vertex_it.first[*floor1_node];
                    typename graph_traits< Graph >::adjacency_iterator ai2, ai_end2;
                    for (boost::tie(ai2, ai_end2) = adjacent_vertices(v, g); ai2 != ai_end2;++ai2)
                        adj_nodes.push_back(get(vertex_id, *ai2));

                    for (auto floor2_node = level_and_nodes[i + 1].begin(); floor2_node != level_and_nodes[i+1].end(); floor2_node++)
                    {
                        if(std::find(adj_nodes.begin(), adj_nodes.end(), *floor2_node) != adj_nodes.end())//if vertex is adj to source vertex send the flow
                        {
                            nodes_and_spath_temp[*floor2_node] += nodes_and_spath_temp[*floor1_node];
                        }
                    }
                }
            }
        }

        double flow_at_each_node[num_vertices(g)];
        std::fill_n(flow_at_each_node, num_vertices(g), 1.0);
        std::map<std::pair<int,int>,double> edges_and_their_flow;
        //calculates edge scores (works surprisingly well)
        for (int i = *max_levels; i >= 0 ; i--) {
            if(i!=0)
            for (auto last_floor_node = level_and_nodes[i].begin(); last_floor_node != level_and_nodes[i].end(); last_floor_node++)
            {
                int distanceRequired = distances_temp[*last_floor_node]-1;
                std::vector<int> flow_targets;
                //get adjacent nodes to the current node
                Vertex v = vertex_it.first[*last_floor_node];
                typename graph_traits< Graph >::adjacency_iterator ai2, ai_end2;
                for (boost::tie(ai2, ai_end2) = adjacent_vertices(v, g); ai2 != ai_end2;++ai2)
                    flow_targets.push_back(get(vertex_id, *ai2));
                for(auto it2 = flow_targets.begin(); it2 != flow_targets.end();)
                {
                    //if adj node isn't 1 distance above current node remove it (we aren't sending flow there)
                    if(distances_temp[*it2]!=distanceRequired)
                    {
                        it2 = flow_targets.erase(it2);
                    }
                    else
                        ++it2;
                }
                for(auto it2 = flow_targets.begin(); it2 != flow_targets.end();it2++)
                {
                    double sPaths_curr_node = nodes_and_spath_temp[*last_floor_node];
                    double sPaths_targ_node = nodes_and_spath_temp[*it2];
                    double ratio = sPaths_targ_node/sPaths_curr_node;
                    flow_at_each_node[*it2] += flow_at_each_node[*last_floor_node]*ratio;

                    std::pair<int,int> edge_points{*last_floor_node,*it2};
                    double edge_value = flow_at_each_node[*last_floor_node]*ratio;
                    edges_and_their_flow.insert(std::pair<std::pair<int,int>,double>{edge_points,edge_value});

                }
            }
        }
        /*
        //print out edge scores(debugging)
        const char name[] = "ABCDEFG";
        std::cout << "mymap contains:\n";
        for (auto it=edges_and_their_flow.begin(); it!=edges_and_their_flow.end(); ++it)
            std::cout <<"("<< name[it->first.first] <<","<< name[it->first.second]<<") => " << it->second << '\n';
        */

        //add the scores from our local container(edges_and_their_flow) to the bundled property of the appropriate edge of our graph
        for (auto it=edges_and_their_flow.begin(); it!=edges_and_their_flow.end(); ++it)
        {
            g[boost::edge(it->first.first,it->first.second,g).first].score +=it->second;
        }
        std::fill_n(d, num_vertices(g), 0);
    }

}
// set all edge scores to 0, expect to use this when removing edges and having to recalculate edge scores
void graphing::clear_edge_scores(graphing::Graph& g) {
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    {
        g[*ei].score=0;
    }
}

void graphing::remove_highest_edge_scores(graphing::Graph& g) {
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);
    graph_traits<Graph>::edge_iterator ei, ei_end;
    std::map<std::pair<int,int>,double> edges_and_their_flow;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    {
        std::pair<int,int> vertex_pair (index[source(*ei, g)],index[target(*ei, g)]);
        edges_and_their_flow.insert(std::pair<std::pair<int,int>,double>{vertex_pair,g[*ei].score});
    }
    // Use a Multimap to sort the map by value
    std::multimap<double,std::pair<int,int>> MM;
    typedef std::multimap<double,std::pair<int,int>>::iterator MMAPIterator;
    for (auto& it : edges_and_their_flow) {
        MM.insert({ it.second, it.first });
    }
    //the last key value will have the edges with the largest edge values
    std::pair<MMAPIterator, MMAPIterator> result = MM.equal_range(MM.rbegin()->first);
    std::cout << "removing Edges with highest betweeness..." <<std::endl;
    for (MMAPIterator it = result.first; it != result.second; it++)
    {
        std::cout << "removing ("<<it->second.first<<","<<it->second.second<<")" << std::endl;
        remove_edge(it->second.first,it->second.second,g);
    }
    std::cout << std::endl;

}

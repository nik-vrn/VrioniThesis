#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <Structs/Graphs/dynamicGraph.h>
#include <Structs/Graphs/forwardStarImpl.h>
#include <Structs/Graphs/adjacencyListImpl.h>
#include <Structs/Graphs/packedMemoryArrayImpl.h>
#include <boost/utility/enable_if.hpp>
#include <iomanip>
#include <Utilities/geographic.h>
#include <Utilities/graphGenerators.h>
#include <Utilities/timer.h>
#include <boost/program_options.hpp>

struct node;
struct edge;

typedef DynamicGraph< PackedMemoryArrayImpl, node, edge>   pmaGraph;
typedef DynamicGraph< ForwardStarImpl, node, edge>         fsGraph;
typedef DynamicGraph< AdjacencyListImpl, node, edge>       Graph;
typedef Graph::NodeIterator                                NodeIterator;
typedef Graph::EdgeIterator                                EdgeIterator;
typedef Graph::NodeDescriptor                              NodeDescriptor;

typedef GraphReader< Graph> AdjReader;
typedef GraphReader< pmaGraph> PmaReader;
typedef GraphReader< fsGraph> FsReader;


struct ShortcutData{

    NodeIterator head;
    unsigned int weight;
};

struct node: DefaultGraphItem
{
    node():x(0),y(0),dist(std::numeric_limits<unsigned int>::max()),pqitem(0),pred(0),timestamp(0)
    {
    }
    int id;

    unsigned int x, y;
    unsigned int dist;
    unsigned int pqitem;

    void* pred;
    unsigned int timestamp;

    std::vector<ShortcutData> forwardEdges;
    std::vector<ShortcutData> backwardEdges;


};

struct edge: DefaultGraphItem
{
    edge():weight(0)
    {
    }

    int id;
    unsigned int weight;
};


#include <Algorithms/basicGraphAlgorithms.h>
#include <Algorithms/ShortestPath/Astar/uniDijkstra.h>

#include "../incl/vector_io.h"
#include "../incl/timer.h"
#include "../incl/contraction_hierarchy.cpp"
#include "../incl/min_max.h"
#include "../incl/inverse_vector.h"


namespace po = boost::program_options;

std::string basePath, map;


const char* graphTypeLabel[4] = { "", "Adjacency", "PMA", "ForwardStar"};

unsigned int graphVariant;

unsigned int algTimestamp;

unsigned int format;
using RoutingKit::ContractionHierarchy;


//Template function to calculate weights of graph edges
template < typename GraphType>
void calcWeights( GraphType& G)
{
    typedef typename GraphType::NodeIterator NodeIterator;
    typedef typename GraphType::EdgeIterator EdgeIterator;
    typedef typename GraphType::InEdgeIterator InEdgeIterator;
    typedef typename GraphType::SizeType SizeType;

    NodeIterator u,v,lastnode;
    EdgeIterator e,lastedge;
    InEdgeIterator k;

    unsigned int max = 0;
    //unsigned int max_ux, max_uy, max_vx, max_vy, u_id, v_id;

    std::stringstream sstr;
    sstr << "Calculating weights of " << G.getNumEdges() << " edges";

    ProgressBar edge_progress( G.getNumEdges(),sstr.str());
    Timer timer;
    timer.start();

    // Iterate over nodes and edges to calculate weights
    for( u = G.beginNodes(), lastnode = G.endNodes(); u != lastnode; ++u)
    {
	    for( e = G.beginEdges(u), lastedge = G.endEdges(u); e != lastedge; ++e)
	    {

            v = G.target(e);
            k = G.getInEdgeIterator(e);

            e->weight = euclideanDistance( u->x, u->y, v->x, v->y);
            k->weight = e->weight;

            if( e->weight > max)
            {
                //max_ux = u->x;
                //max_uy = u->y;
                //max_vx = v->x;
                //max_vy = v->y;
                //u_id = G.getRelativePosition(u);
                //v_id = G.getRelativePosition(v);
                max = e->weight;
            }
            ++edge_progress;
            //std::cout << u->lat << " " << v->lat << std::endl;
        }
    }

    std::cout << "\tMax weight: " << max << "\tTime:\t" << timer.getElapsedTime() << std::endl;
    //std::cout << "( " << max_ux << ", " << max_uy << ") " << u_id << "->" << v_id << " ( " << max_vx << ", " << max_vy << ")\n";
}



//Function for preprocessing the graph
void preprocessing(Graph& G, std::string& ch_file){

    try{
        auto ch = ContractionHierarchy::build(G, [](std::string msg){std::cout << msg << std::endl;});
        check_contraction_hierarchy_for_errors(ch);
        ch.save_file(ch_file);
    }catch(std::exception&err){
        std::cerr << "Stopped on exception : " << err.what() << std::endl;
    }
}


int main( int argc, char* argv[])
{
    std::string ch_file = "ch.dat";
    graphVariant = 0;
    format = 0;

    AdjReader* reader = 0;
    PmaReader* pmaReader = 0;
    FsReader* fsReader = 0;

    basePath = std::string(getenv("HOME")) + "/Projects/Graphs/";
    map = "luxembourg";

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("graphtype,g", po::value< unsigned int>(), "graph type. All[0], Adjacency List[1], Packed Memory Graph[2], forward Star[3]. Default:0")
		("format,f", po::value< unsigned int>(), "map format. DIMACS10[0], DIMACS9[1]. Default:0")
        ("map,m", po::value< std::string>(), "input map. The name of the map to read. Default:'luxembourg'. Maps should reside in '$HOME/Projects/Graphs/DIMACS{9,10}/' and should consist of 2 files, both with the same map name prefix, and suffixes 'osm.graph' and 'osm.xyz' in the case of DIMACS10 and, '.gr' and '.co' in the case of DIMACS9.")
        ("queries,q", "The input queries are in the file queries.dat. The file must be placed in '$HOME/Projects/Graphs/.../$MAP/'")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if ( vm.empty())
    {
        std::cout << desc << "\n";
        return 0;
    }


    if ( vm.count("graphtype"))
    {
        graphVariant = vm["graphtype"].as<unsigned int>();
    }

    if ( vm.count("map"))
    {
        map = vm["map"].as<std::string>();
    }

    if ( vm.count("format") && ( format = vm["format"].as<unsigned int>()) == 1)
    {
        basePath = basePath + "DIMACS9/" + map + "/";

        reader = new DIMACS9Reader<Graph>( basePath + map + ".gr",
                                           basePath + map + ".co");

        pmaReader = new DIMACS9Reader<pmaGraph>( basePath + map + ".gr",
                                                 basePath + map + ".co");

        fsReader = new DIMACS9Reader<fsGraph>( basePath + map + ".gr",
                                               basePath + map + ".co");
    }

    else
    {
        basePath = basePath + "DIMACS10/" + map + "/";

        reader = new DIMACS10Reader<Graph>( basePath + map + ".osm.graph",
                                            basePath + map + ".osm.xyz");

        pmaReader = new DIMACS10Reader<pmaGraph>( basePath + map + ".osm.graph",
                                                  basePath + map + ".osm.xyz");

        fsReader = new DIMACS10Reader<fsGraph>( basePath + map + ".osm.graph",
                                                  basePath + map + ".osm.xyz");
    }


    Graph G;
    pmaGraph pmaG;
    fsGraph  fsG;

    Timer timer;
    timer.start();
    G.read(reader);

     // Assign IDs to the nodes
    int nodeId = 0;
    for( NodeIterator u = G.beginNodes(), lastnode = G.endNodes(); u != lastnode; ++u)
        u->id = nodeId++;


    if(format == 0)
        calcWeights(G);
    timer.stop();
    std::cout << "Graph has " << (double)G.memUsage() / 1048576.0 << " Mbytes (without landmarks). Time spent to read:\t"
              << timer.getElapsedTime() << "sec" << std::endl;
    if(format == 0)
        calcWeights(G);

    long long  time_max = 0;
    long long  time_sum = 0;

    Timer prep_timer;
    prep_timer.start();

    //Calling preprocessing function
    preprocessing(G, ch_file);

    prep_timer.stop();
    long long time = prep_timer.getElapsedTimeInMilliSec();
    time_max = std::max(time_max, time);

    std::cout << "max running time: " << time_max << "musec" << std::endl;

    G.clear();

   return 0;
}

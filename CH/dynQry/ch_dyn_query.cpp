#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <map>
#include <Structs/Graphs/dynamicGraph.h>
#include <Structs/Graphs/forwardStarImpl.h>
#include <Structs/Graphs/adjacencyListImpl.h>
#include <Structs/Graphs/packedMemoryArrayImpl.h>
#include <boost/utility/enable_if.hpp>
#include <iomanip>
#include <chrono>
#include <Utilities/geographic.h>
#include <Utilities/graphGenerators.h>
#include <Utilities/timer.h>
#include <boost/program_options.hpp>
#include <vector>
#include <algorithm>
#include <random>



struct node;
struct edge;

typedef DynamicGraph< PackedMemoryArrayImpl, node, edge>   pmaGraph;
typedef DynamicGraph< ForwardStarImpl, node, edge>         Graph;//fsGraph;
typedef DynamicGraph< AdjacencyListImpl, node, edge>       AdjGraph;//Graph;
typedef Graph::NodeIterator                                NodeIterator;
typedef Graph::EdgeIterator                                EdgeIterator;
typedef Graph::NodeDescriptor                              NodeDescriptor;

typedef GraphReader< AdjGraph> AdjReader;
typedef GraphReader< pmaGraph> PmaReader;
typedef GraphReader< Graph> FsReader;


struct ShortcutData{

    ShortcutData(): weight(0)
    {}

    ShortcutData(NodeIterator head, unsigned int weight):
    head(head), weight(weight)
    {}

    ShortcutData( const ShortcutData& sd):
    head(sd.head), weight(sd.weight)
    {}

    void operator=( const ShortcutData& sd)
    {
      head=sd.head;
      weight=sd.weight;
    }

    NodeIterator head;
    unsigned int weight;

};

struct node: DefaultGraphItem
{
    node():x(0),y(0),dist(std::numeric_limits<unsigned int>::max()),pqitem(0),pred(0),timestamp(0),rank(0)
    {
    }
    int id;

    unsigned int x, y;
    unsigned int dist,distBack;
    unsigned int pqitem,pqitemBack;

    void* pred;
    void* succ;
    unsigned int timestamp;
    unsigned int timestampBack;

    int rank; //gamma rank
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
    unsigned int shortcutWeight;
};


#include <Algorithms/basicGraphAlgorithms.h>
#include <Algorithms/ShortestPath/Astar/uniDijkstra.h>
#include <Algorithms/basicGraphAlgorithms.h>
#include <Algorithms/ShortestPath/bidDijkstra.h>
//#include <Algorithms/ShortestPath/bidirectionalDijkstra.h>

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

using namespace RoutingKit;

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
    for( NodeIterator u = G.beginNodes(), lastnode = G.endNodes(); u != lastnode; ++u)
    {
        for( EdgeIterator e = G.beginEdges(u), lastedge = G.endEdges(u); e != lastedge; ++e)
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

template <typename GraphType>
class ChUpdater {
public:
    typedef typename GraphType::NodeIterator NodeIterator;
    typedef typename GraphType::EdgeIterator EdgeIterator;
    typedef typename GraphType::InEdgeIterator InEdgeIterator;
    typedef typename GraphType::SizeType SizeType;

    ChUpdater(GraphType& G, ContractionHierarchy& ch) : G(G), ch(ch) {}


    struct EdgeItem
    {
        EdgeItem() : isShortcut(false), edge(nullptr) {}
        EdgeItem(EdgeIterator edge, bool isShortcut) : edge(edge), isShortcut(isShortcut) {}
        EdgeItem(const EdgeItem& eIt) : edge(eIt.edge), isShortcut(eIt.isShortcut) {}

        EdgeIterator edge;
        ShortcutData shtEdge;
        bool isShortcut;
    };


    void update( EdgeIterator affectedEdge, int weightDecrease) {

        NodeIterator sourceNode = G.source( affectedEdge);
        NodeIterator targetNode = G.target( affectedEdge);

        if( ( affectedEdge->weight + weightDecrease) < 0)
        {
            std::cout << "Warning negative weight, aborting...\n";
            return;
        }

        int oldWeight = affectedEdge->weight;
        affectedEdge->weight += weightDecrease;

        if( affectedEdge->shortcutWeight > affectedEdge->weight)
          affectedEdge->shortcutWeight  = affectedEdge->weight;
        else
          return;

        std::queue<EdgeItem> Q;

        EdgeItem eIt;
        eIt.isShortcut = false;
        eIt.edge = affectedEdge;

        Q.push( eIt);
        while( !Q.empty())
        {
            eIt = Q.front();
            Q.pop();

            if( eIt.isShortcut == false)
            {
                //typedef PriorityQueue<int, EdgeIterator, HeapStorage>      PriorityQueueType;
                //PriorityQueueType pq;

                //pq.insert( sourceNode->rank, sourceNode, &(sourceNode->pqitem));

                EdgeIterator uv =  eIt.edge;
                //pq.popMin();
                //EdgeIterator e;
                NodeIterator u = G.source(affectedEdge);
                NodeIterator v = G.target(affectedEdge);

                for ( ShortcutData& shU : u->forwardEdges)
                {
                    NodeIterator x = shU.head;

                    for ( ShortcutData& shX : x->forwardEdges)
                    {
                        NodeIterator y = shX.head;
                        if( y == v)
                        {
                            examine(Q, uv, shU, shX);
                        }
                    }

                    for (EdgeIterator e = G.beginEdges(x), lastedge = G.endEdges(x); e != lastedge; ++e)
                    {
                        NodeIterator y = G.target(e);
                        if( y == v)
                        {
                            examine(Q, uv, shU, e);
                        }
                    }
                }

                for (EdgeIterator e = G.beginEdges(u), lastedge = G.endEdges(u); e != lastedge; ++e)
                {
                    NodeIterator x = G.target( e);

                    for ( ShortcutData& shX : x->forwardEdges)
                    {
                        NodeIterator y = shX.head;
                        if( y == v)
                        {
                            examine(Q, uv, e, shX);
                        }
                    }

                    for (EdgeIterator h = G.beginEdges(x), lastedge = G.endEdges(x); h != lastedge; ++h)
                    {
                        NodeIterator y = G.target( e);
                        if( y == v)
                        {
                            examine(Q, uv, e, h);
                        }
                    }
                }

            for ( ShortcutData& shV : v->backwardEdges)
            {
                NodeIterator x = shV.head;

                for ( ShortcutData& shX : x->backwardEdges)
                {
                    NodeIterator y = shX.head;
                    if( y == u)
                    {
                        examine(Q, uv, shX, shV);
                    }
                }

                for( EdgeIterator e = G.beginEdges(u), lastedge = G.endEdges(u); e != lastedge; ++e)
                {
                    NodeIterator y = G.target(e);

                    if( y == u)
                    {
                        InEdgeIterator k = G.getInEdgeIterator(e);
                        examine(Q, uv, e, shV);
                    }
                }
            }

            for( EdgeIterator e = G.beginEdges(u), lastedge = G.endEdges(u); e != lastedge; ++e)
            {
                NodeIterator x = G.target(e);

                for ( ShortcutData& shV : v->backwardEdges)
                {
                    NodeIterator y = shV.head;
                    if( y == u)
                    {
                       InEdgeIterator k = G.getInEdgeIterator(e);
                       examine(Q, uv, e, shV);
                    }
                }
            }
          }

        }
    }

    void examine( std::queue<EdgeItem>& Q, EdgeIterator& uv, EdgeIterator& ux, EdgeIterator& xv)
    {
        if( (uv->shortcutWeight + ux->shortcutWeight) >= xv->shortcutWeight)
          return;
        else
        {
            xv->shortcutWeight = uv->shortcutWeight + ux->shortcutWeight;
            EdgeItem eIt;
            eIt.isShortcut = false;
            eIt.edge = xv;
            Q.push( eIt);
          }
    }

    void examine( std::queue<EdgeItem>& Q, EdgeIterator& uv, ShortcutData& ux, ShortcutData& xv)
    {
        if ((uv->shortcutWeight + ux.weight) >= xv.weight)
          return;
        else
        {
            xv.weight = uv->shortcutWeight + ux.weight;
            EdgeItem eIt;
            eIt.isShortcut = true;
            eIt.shtEdge = xv;
            Q.push(eIt);
          }
    }

    void examine(std::queue<EdgeItem>& Q, EdgeIterator& uv, EdgeIterator& ux, ShortcutData& xv)
    {

        if ((uv->shortcutWeight + ux->shortcutWeight) >= xv.weight)
          return;
        else
        {
            xv.weight = uv->shortcutWeight + ux->shortcutWeight;
            EdgeItem eIt;
            eIt.isShortcut = true;
            eIt.shtEdge = xv;
            Q.push(eIt);

          }
       }

    void examine(std::queue<EdgeItem>& Q, EdgeIterator& uv, ShortcutData& ux, EdgeIterator& xv)
    {
        if ((uv->shortcutWeight + ux.weight) >= xv->shortcutWeight)
          return;
        else
        {
            xv->shortcutWeight = uv->shortcutWeight + ux.weight;
            EdgeItem eIt;
            eIt.isShortcut = false;
            eIt.edge = xv;
            Q.push(eIt);
        }
    }

private:
    GraphType& G;
    ContractionHierarchy& ch;
};

std::vector<int> generate_random_binary_array(int arrayLength, int numberOfOnes) {
  // std::cout << "/* ROUTINE: generate_random_binary_array */" << std::endl;

    std::vector<int> randomBinaryArray;

    randomBinaryArray.reserve(100);
    for (int i = 0; i < 100; ++i)
    {
      randomBinaryArray.push_back(0);
    }

    if (numberOfOnes <= 0) {
        //Construct an ALL-ZEROS vector if numberOfOnes is less than or equal to 0
        randomBinaryArray = std::vector<int>(arrayLength, 0);
    } else if (numberOfOnes >= arrayLength) {
        // Construct an ALL-ONES vector if numberOfOnes is greater than or equal to arrayLength
        randomBinaryArray = std::vector<int>(arrayLength, 1);
    } else {
        // Reserve space,Insert 0s at the end of the vector,Insert 1s at the end of the vector
        randomBinaryArray.reserve(arrayLength);
        randomBinaryArray.insert(randomBinaryArray.end(), arrayLength - 100, 0);
        randomBinaryArray.insert(randomBinaryArray.end(), numberOfOnes, 1);

        // Use random_device to obtain a seed for the random number engine
        std::random_device rd;

        //Create a random number engine and seed it with rd()
        std::default_random_engine eng(rd());

        // Shuffle the vector starting from the 100th element to the end
        std::shuffle(randomBinaryArray.begin() + 100, randomBinaryArray.end(), eng);
    }

    return randomBinaryArray;
}

std::ofstream openOutputFile(const std::string& filename) {
    std::ofstream outputFile(filename, std::ios::out | std::ios::app);
    if (!outputFile.is_open()) {
        std::cout << "Error opening the file: " << filename << std::endl;
        exit(1);
    }
    return outputFile;
}

std::vector<int> readFromFile(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<int> data;

    if (file.is_open()) {
        int value;
        while (file >> value) {
            data.push_back(value);
        }

        file.close();
    } else {
        std::cout << "Error opening file: " << filename << std::endl;
    }

    return data;
}

int main( int argc, char* argv[])
{
    graphVariant = 0;
    format = 0;

    std::string ch_file="ch.dat";
    std::string source_file = "src.txt";
    std::string target_file = "trg.txt";
    std::string distance_file = "res.txt";

    AdjReader* reader = 0;
    PmaReader* pmaReader = 0;
    FsReader* fsReader = 0;

    basePath = std::string(getenv("HOME")) + "/Projects/Graphs/";
    map = "luxembourg";

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("chfile,c", po::value<std::string>(), "contraction hierarchy file. Default:'ch.dat'")
        ("graphtype,g", po::value< unsigned int>(), "graph type. All[0], Adjacency List[1], Packed Memory Graph[2], forward Star[3]. Default:0")
        ("format,f", po::value< unsigned int>(), "map format. DIMACS10[0], DIMACS9[1]. Default:0")
        ("src,s", po::value< std::string>(), "input file containing source nodes. Default: 'src.txt'")
        ("trg,t", po::value< std::string>(), "input file containing target nodes. Default: 'trg.txt'")
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

    if (vm.count("chfile"))
    {
        ch_file = vm["chfile"].as<std::string>();
    }

    if ( vm.count("graphtype"))
    {
        graphVariant = vm["graphtype"].as<unsigned int>();
    }

    if ( vm.count("map"))
    {
        map = vm["map"].as<std::string>();
    }

    if ( vm.count("src"))
    {
        source_file = vm["src"].as<std::string>();
    }

    if ( vm.count("trg"))
    {
        target_file = vm["trg"].as<std::string>();
    }

    if ( vm.count("format") && ( format = vm["format"].as<unsigned int>()) == 1)
    {
        basePath = basePath + "DIMACS9/" + map + "/";

        reader = new DIMACS9Reader<AdjGraph>( basePath + map + ".gr",
                                           basePath + map + ".co");

        pmaReader = new DIMACS9Reader<pmaGraph>( basePath + map + ".gr",
                                                 basePath + map + ".co");

        fsReader = new DIMACS9Reader<Graph>( basePath + map + ".gr",
                                               basePath + map + ".co");
    }

    else
    {
        basePath = basePath + "DIMACS10/" + map + "/";

        reader = new DIMACS10Reader<AdjGraph>( basePath + map + ".osm.graph",
                                            basePath + map + ".osm.xyz");

        pmaReader = new DIMACS10Reader<pmaGraph>( basePath + map + ".osm.graph",
                                                  basePath + map + ".osm.xyz");

        fsReader = new DIMACS10Reader<Graph>( basePath + map + ".osm.graph",
                                                  basePath + map + ".osm.xyz");
    }


    Graph adjG;
    pmaGraph pmaG;


    Timer timer;
    timer.start();

    //adjacencyList
    //Graph  G;
    //G.read(reader);

    //forward star
    Graph  G;
    fsReader->read( G);

    std::map<int, NodeIterator> nodeMap;

    unsigned int timestamp = 0;

    // Assign IDs to the nodes and filling the map
    int nodeId = 0;
    for( NodeIterator u = G.beginNodes(), lastnode = G.endNodes(); u != lastnode; ++u)
    {
        u->id = nodeId++;
        nodeMap[u->id] = u;
    }

    std::vector<unsigned> tail;
    std::vector<unsigned> head;
    std::vector<unsigned> weights;
    NodeIterator u,v,endNode;
    EdgeIterator e,endEdge;

      for( u = G.beginNodes(), endNode = G.endNodes(); u != endNode; ++u)
      {

  	    for( e = G.beginEdges(u), endEdge = G.endEdges(u); e != endEdge; ++e)
  	    {
  	    	v = G.target(e);
  	    	tail.push_back( u->id);
  	    	head.push_back( v->id);
  	    	weights.push_back( e->weight);
        //  std::cout << u->id <<  " -> " << v->id << std::endl;
  	    }
  	}

    // Calling load_file from ContractionHierarchy
     ContractionHierarchy ch = ContractionHierarchy::load_file(ch_file);

     Dijkstra<Graph> dijkstra(G,&timestamp);
     BidirectionalDijkstra<Graph> bidirectionalDijkstra(G,&timestamp);
     ChUpdater up(G, ch);
     ContractionHierarchyQuery ch_query(ch);

   // Calculate and populate forwardEdges and backwardEdges
   //std::cout<<ch.forward.first_out.size();
    if( false)
    {
    for( unsigned x = 0; x < G.getNumNodes(); ++x)
    for (unsigned xy = ch.forward.first_out[x]; xy < ch.forward.first_out[x+1]; ++xy)
    {

    	unsigned shweight = ch.forward.weight[xy];
    	int orArc1Id = ch.forward.shortcut_first_arc[xy];    	
    	//int tailNodeId =  ch.forward.tail[i];
    	int tailNodeId = x;    	
    	int headNodeId = ch.forward.head[xy];
    	
    	unsigned y = ch.forward.head[xy];
				        std::cout << x <<" "  << y<< " "<< ch.forward.weight[xy];
				    	std::cout <<" "  << ch.rank[x]  <<std::endl;
				    	std::cout<< ch.rank[y]<<std::endl;
				    	{int a;std::cin>>a;}	
				    	
    	//auto a = find_arc_given_sorted_head( first_out, head, ch.order[x], ch.order[y]);
    	
    	NodeIterator u = nodeMap[tailNodeId];
	NodeIterator v = nodeMap[headNodeId];    
	
	EdgeIterator e = G.getEdgeIterator( u,v);
	unsigned orweight = e->weight;	
    	
    	u->forwardEdges.emplace_back(v,shweight);
           unsigned distance = dijkstra.runQuery(u, v);  
           std::cout << "Error: " << "distance: " << distance << " " << "shweight: " << shweight << " " << "orweight: " << orweight << std::endl;  
            
          if (distance != shweight) {
            std::cout << "Error: " << "distance: " << distance << " " << "shweight: " << shweight << std::endl;
            {int a;std::cin>>a;}
          }            	
    }   
    }
    
    //exit(1);
   
   
    if( false)
    {
    for (unsigned i = 0; i < ch.forward.first_out.size(); ++i)
    {
    	//int shArcId = ch.forward.first_out[i];
    	if( !ch.forward.is_shortcut_an_original_arc.is_set(i))
    	{
    	unsigned weight = ch.forward.weight[i];
    	int orArc1Id = ch.forward.shortcut_first_arc[i];
    	int orArc2Id = ch.forward.shortcut_second_arc[i];
    	int tailNodeId = tail[orArc1Id];
    	int headNodeId = head[orArc2Id];
    	NodeIterator u = nodeMap[tailNodeId];
	NodeIterator v = nodeMap[headNodeId];    	
    	
    	u->forwardEdges.emplace_back(v,weight);
           unsigned distance = dijkstra.runQuery(u, v);    
            
          if (distance != weight) {
            std::cout << "Error: " << "distance: " << distance << " " << "shweight: " << weight << std::endl;
            {int a;std::cin>>a;}
          }            	
        }
    }
    }

   if(false)
   {
    for (const auto& pair : nodeMap)
    {
        NodeIterator u = pair.second;

        for (unsigned i = ch.forward.first_out[u->id]; i < ch.forward.first_out[u->id + 1]; ++i)
        {
            unsigned v ;//= input_arc_id[ch.forward.shortcut_second_arc[i]];//ch.forward.head[i];
            
	    	std::cout <<" " << v;
	    	std::cout <<" " << head[v];	    	
		std::cout <<" " << ch.forward.head[i];
		std::cout << " " ;
             v= head[v];
            /*if(ch.forward.is_shortcut_an_original_arc.is_set()){
				return input_weight[ch.forward.shortcut_first_arc[a]];
			}else{
				return link(get_backward_weight(ch.forward.shortcut_first_arc[a]), get_forward_weight(ch.forward.shortcut_second_arc[a]));
			}*/
			
			//std::cout << ch.forward.head.size() << " " << ch.forward.shortcut_first_arc.size();
            
            unsigned weight = ch.forward.weight[i];
            u->forwardEdges.emplace_back(nodeMap[v],weight);
            unsigned distance = dijkstra.runQuery(u, nodeMap[v]);


	if( ch.forward.is_shortcut_an_original_arc.is_set(i))
	{
		std::cout  <<" original ";
	}else
	{
	std::cout  <<" shortcut ";
	}
	
		std::cout << ch.forward.shortcut_first_arc[i];    //arc id
		std::cout <<" "<< ch.forward.shortcut_second_arc[i];	 //tai;
		
	    	int orArc1Id = ch.forward.shortcut_first_arc[i];	
	    	std::cout <<" " << tail[orArc1Id];
		std::cout <<" " << head[orArc1Id];	    	
	    	std::cout <<" " << tail[i];
		std::cout <<" " << head[i];		
		
	    	int orArc2Id = ch.forward.shortcut_second_arc[i];		
		std::cout <<" w> " << weights[ch.forward.shortcut_first_arc[i]];
		std::cout <<" " << weights[ch.forward.shortcut_second_arc[i]];	
		std::cout <<" " << weights[tail[orArc1Id]];	
		std::cout <<" " << weights[head[orArc1Id]];	// ch.forward.shortcut_first_arc[i]] ==  head[orArc1Id]
		std::cout <<" " << weights[tail[orArc2Id]];	
		std::cout <<" " << weights[head[orArc2Id]];									
		std::cout << " " <<weight;	
		std::cout << " " <<distance;	
			std::cout << " " ;
		
		{int a;std::cin>>a;}	
	
          if (distance != weight) {
            ;//std::cout << "Error: " << "distance: " << distance << " " << "shweight: " << weight << std::endl;
            //{int a;std::cin>>a;}
          }
            //std::cout << "F";
        }

        for (unsigned i = ch.backward.first_out[u->id]; i < ch.backward.first_out[u->id + 1]; ++i)
        {
            unsigned v = ch.backward.head[i];
            unsigned weight = ch.backward.weight[i];
            u->backwardEdges.emplace_back(nodeMap[v],weight);
        }
    }
    }
    
     std::map<int,int> rank2Id;
            for (unsigned j = 0; j < ch.rank.size(); ++j)
            {
	        ;//rank2Id
            	//std::cout<< " " <<ch.rank[j]; 
            }    
            
    
    for (const auto& pair : nodeMap)
    {
        NodeIterator u = pair.second;
        unsigned uR = ch.rank[u->id];
    	rank2Id[uR] = u->id;
    } 
    
    for (const auto& pair : nodeMap)
    {
        NodeIterator u = pair.second;
        unsigned uR = ch.rank[u->id];

        for (unsigned i = ch.forward.first_out[uR]; i < ch.forward.first_out[uR + 1]; ++i)
        {
  
            unsigned vR = ch.forward.head[i];
            int vId = rank2Id[vR];
            
            unsigned weight = ch.forward.weight[i];
            u->forwardEdges.emplace_back(nodeMap[vId],weight);
   
        }

        for (unsigned i = ch.backward.first_out[uR]; i < ch.backward.first_out[uR + 1]; ++i)
        {
            unsigned vR = ch.backward.head[i];
            int vId = rank2Id[vR];
            
            unsigned weight = ch.backward.weight[i];
            u->backwardEdges.emplace_back(nodeMap[vId],weight);

        }
    }    

    if(format == 0)
        calcWeights(G);
    timer.stop();
    std::cout << "Graph has " << (double)G.memUsage() / 1048576.0 << " Mbytes (without landmarks). Time spent to read:\t"
              << timer.getElapsedTime() << "sec" << std::endl;
    if(format == 0)
        calcWeights(G);

    for( NodeIterator u = G.beginNodes(), lastnode = G.endNodes(); u != lastnode; ++u)
    {
      for( EdgeIterator e = G.beginEdges(u), lastedge = G.endEdges(u); e != lastedge; ++e)
      {
        e->shortcutWeight = e->weight;
      }
    }

    std::vector<unsigned int>sourceIds;
    std::vector<unsigned int>targetIds;

    int numQueries = 1000;

    int MaxId = G.getNumNodes();

    for( int i=0; i<numQueries; i++)
    {
      int sid =  rand() % MaxId;
      int tid =  rand() % MaxId;

      if( sid == tid)
      {
        i--;
        continue;
      }
      sourceIds.push_back(sid);
      targetIds.push_back(tid);
    }


    std::vector<unsigned>distances(numQueries,0);

    std::cout << "Loaded " << numQueries << " test queries" << std::endl;

    std::cout << "Running test queries ... \n" << std::flush;

    long long time_max = 0;
    long long time_sum = 0;

    std::vector<EdgeIterator> normalEdgesUp;
    std::vector<ShortcutData> shortcutEdgesUp;
    std::vector<EdgeIterator> edges;

    for (NodeIterator u = G.beginNodes(), lastnode = G.endNodes(); u != lastnode; ++u)
        for (EdgeIterator e = G.beginEdges(u), lastedge = G.endEdges(u); e != lastedge; ++e)
            edges.push_back(e);


    std::queue<NodeIterator> Q;
    NodeIterator s = G.beginNodes();
    Q.push( s);
    std::map<int ,bool> isVisited;
    isVisited[s->id] = true;

    while( !Q.empty())
    {
        NodeIterator u = Q.front();
        NodeIterator v = Q.front();
        Q.pop();

        for (const ShortcutData& shortcut : u->forwardEdges) {
            NodeIterator v = shortcut.head;
            v->rank = std::max( u->rank + 1, v->rank);
            if( !isVisited[v->id])
            {
                isVisited[v->id] = true;
                Q.push(v);
            }
        }

        for (const ShortcutData& shortcut : v->backwardEdges) {
            NodeIterator u = shortcut.head;
            u->rank = std::max( v->rank + 1, u->rank);
            if( !isVisited[u->id])
            {
                isVisited[u->id] = true;
                Q.push(u);
            }
        }

        for (EdgeIterator e = G.beginEdges(u), lastedge = G.endEdges(u); e != lastedge; ++e)
        {
            NodeIterator v = G.target( e);
            if( !isVisited[v->id])
            {
                isVisited[v->id] = true;
                Q.push(v);
            }
        }
    }

    int arrayLength = 1000;
    int numberOfOnes = 100;

    std::vector<int> result = generate_random_binary_array(arrayLength , numberOfOnes);
/*
    for (int i = 0; i < result.size(); ++i) {
        std::cout << result[i];
    }
     std::cout << std::endl;
*/
    std::ofstream outputFile = openOutputFile("output.txt");
    std::ofstream edgesFile = openOutputFile("edges.txt");
    std::ofstream totalDistance = openOutputFile("totalDistance.txt");
    
    // Calculate the weight decrease factor
    double weightDecreaseFactor = 0.1;
    std::cout << "weightDecreaseFactor: " << weightDecreaseFactor << std::endl;

    int numEdges = G.getNumEdges();
    int numUpdates = 1000;

    std::ofstream distanceFile = openOutputFile(distance_file);
    std::ofstream sourceFile = openOutputFile(source_file);
    std::ofstream targetFile = openOutputFile(target_file);

    std::vector<int> sources = readFromFile(source_file);
    std::vector<int> targets = readFromFile(target_file);
    std::vector<int> edgesData = readFromFile("edges.txt");
/*
//Update
    for (unsigned i = 0; i < numUpdates; ++i) {
   
      int randomEdgeIndex = rand() % numEdges;
      EdgeIterator modifiedEdge = edges[randomEdgeIndex];

      int weightDecrease = - static_cast<int>(modifiedEdge->weight * weightDecreaseFactor) + modifiedEdge->weight;

      if (weightDecrease <= 0) {
         std::cout << "Weight decrease is zero or negative. Skipping update " << (i + 1) << std::endl;
         continue;
       }
       outputFile << "Update " << (i + 1) << " - Decrease: " << weightDecrease << std::endl;
       edgesFile << randomEdgeIndex << std::endl;

       Timer timer;
       timer.start();

       up.update(modifiedEdge, weightDecrease);

       timer.stop();
       long long time = timer.getElapsedTimeInMilliSec() * 1000;
       time_max = std::max(time_max, time);
       time_sum += time;
   }
   std::cout << "max running time: " << time_max << "musec" << std::endl;
   std::cout << "avg running time: " << time_sum  << "musec" << std::endl;
   
*/
//Query&Update
for (unsigned i = 0; i < result.size(); ++i) {

    int randomEdgeIndex = edgesData[i];

    if (result[i] == 0) {

      std::cout << "Executing Query " << (i + 1) << "..." << std::endl;

      NodeIterator source = nodeMap[sources[i]];
      NodeIterator target = nodeMap[targets[i]];

      if (sources[i] == 0 || targets[i] == 0) {
       std::cout << "Invalid source or target node ID in query " << (i + 1) << std::endl;
     }
      else {
	  unsigned kitDstance = ch_query.reset().add_source(source->id).add_target(target->id).run().get_distance(); 

          Timer timer;
          timer.start();

          unsigned chDistance = bidirectionalDijkstra.runQuery(source, target);
          //unsigned distance = dijkstra.runQuery(source, target);
	 
          if (kitDstance != chDistance) {
           std::cout << "Error: " <<" " << "chDistance: " << chDistance << " " << "kitDstance: " << kitDstance << std::endl;
          }

          timer.stop();
          long long time = timer.getElapsedTimeInMilliSec() * 1000;
          time_max = std::max(time_max, time);
          time_sum += time;
		
          std::cout << "avg running time: " << time_sum << "musec" << std::endl;
          totalDistance << chDistance << std::endl;

      }
  } else {
      // Update
      std::cout << "Executing Update " << (i + 1) << "..." << std::endl;

      if (randomEdgeIndex < edges.size()) {
      EdgeIterator modifiedEdge = edges[randomEdgeIndex];

      int weightDecrease = -static_cast<int>(modifiedEdge->weight * weightDecreaseFactor) + modifiedEdge->weight;

      if (weightDecrease <= 0) {
          std::cout << "Weight decrease is zero or negative. Skipping update " << (i + 1) << std::endl;
          continue;
      }

      Timer timer;
      timer.start();

      up.update(modifiedEdge, weightDecrease);

      timer.stop();
      long long time = timer.getElapsedTimeInMilliSec() * 1000;
      time_max = std::max(time_max, time);
      time_sum += time;

      std::cout << "avg running time: " << time_sum  << "musec" << std::endl;

    } else {
          std::cout << "EdgeIndex out of bounds. Exiting the loop. Index: " << std::endl;
          break;
          }
        }
    }
    sourceFile.close();
    targetFile.close();
    distanceFile.close();

/*
//Query
   for(unsigned i=0; i < numQueries; ++i){

     NodeIterator source = nodeMap[sourceIds[i]];
     NodeIterator target = nodeMap[targetIds[i]];

     if (sourceIds[i] == 0 || targetIds[i] == 0)
      {
        std::cout << "Invalid source or target node ID in query " << (i + 1) << std::endl;
      }
      else {
         Timer timer;
         timer.start();

         unsigned chDistance = bidirectionalDijkstra.runQuery(source, target);
          //unsigned distance = dijkstra.runQuery(source, target);
	 unsigned kitDstance = ch_query.reset().add_source(source->id).add_target(target->id).run().get_distance(); 

          if (kitDstance != chDistance) {
            std::cout << "Error: " << " " << "chDistance: " << chDistance << " " << "kitDistance: " << kitDstance << std::endl;
          }

          timer.stop();
          long long time = timer.getElapsedTimeInMilliSec() * 1000;
          time_max = std::max(time_max, time);
          time_sum += time;
          
          //DistanceFile << "Query " << (i + 1) << ": Source Node ID = " << sourceIds[i] << ", Target Node ID = " << targetIds[i] << std::endl;
          sourceFile << sourceIds[i] << std::endl;
          targetFile << targetIds[i] << std::endl;
          distanceFile << chDistance << std::endl;

          }
        }
        
   std::cout << "done" << std::endl;
   std::cout << "max running time query: " << time_max << "musec" << std::endl;
   std::cout << "avg running time query: " << time_sum << "musec" << std::endl;
*/
    G.clear();

    return 0;
}

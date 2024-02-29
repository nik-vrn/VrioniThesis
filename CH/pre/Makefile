INCLUDEDIR='/home/vrioni/Projects/pgl/include/'
BOOSTINCLUDEDIR='/usr/include/boost/'
BOOSTLIBDIR='/usr/local/lib'
ROUTINGKITLIBDIR='/home/vrioni/RoutingKit/lib'

all: compile

compile:
	g++ ch.cpp  -O3 -I$(INCLUDEDIR) -I$(BOOSTINCLUDEDIR) -L$(BOOSTLIBDIR) -L$(ROUTINGKITLIBDIR) -DNDEBUG -lroutingkit -lboost_program_options  -static

debug:
	g++ ch.cpp -O0 -g -I$(INCLUDEDIR) -I$(BOOSTINCLUDEDIR) -L$(BOOSTLIBDIR) -L$(ROUTINGKITLIBDIR) -Wall -lroutingkit -lboost_program_options 
	
debug-noassert:
	g++ ch.cpp -O0 -g -I$(INCLUDEDIR) -I$(BOOSTINCLUDEDIR) -L$(BOOSTLIBDIR) -L$(ROUTINGKITLIBDIR) -Wall -lroutingkit -lboost_program_options -DNDEBUG


run:
	./a.out

valgrind:
	valgrind --leak-check=full ./a.out -s 10

profile: debug
	valgrind --tool=callgrind --callgrind-out-file=profile  ./a.out -s 10
	kcachegrind profile &

clean: 
	rm *.out

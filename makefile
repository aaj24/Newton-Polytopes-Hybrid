G=mpic++ -std=c++11 -O3 -w

main: main.cpp PolyDecomp.o convex_hull.o priority_update.o
	$(G) -fopenmp $^ -o $@
PolyDecomp.o: PolyDecomp.cpp PolyDecomp.hpp
	$(G) -fopenmp PolyDecomp.cpp PolyDecomp.hpp -c
convex_hull.o: convex_hull.cpp convex_hull.hpp
	$(G) convex_hull.cpp convex_hull.hpp -c
priority_update.o: priority_update.cpp priority_update.hpp
	$(G) priority_update.cpp priority_update.hpp -c
ArrayX.o: ArrayX.cpp ArrayX.hpp
	$(G) $^ -c
clean:
	rm *.o *.gch main

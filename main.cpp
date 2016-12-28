//
//  main.cpp
//  CMPS 299
//
//  Created by Ali Jbaily on 2/13/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#include "PolyDecomp.hpp"
#include <fstream>
#include <time.h>

vector<Point> circularGenerator(double c, int r);

using namespace std;

int main(int argc, const char * argv[]) {
    int provided;
    MPI_Init(NULL, NULL);
    vector<Point> ConvexHull;
    fstream file;
    file.open(argv[1]);
    if (!file.is_open())
        return 9;
    
    coord_t x, y;
    while (file >> x >> y)
        ConvexHull.push_back({x, y}); // read input
    
    file.close();
    int m = (int)ConvexHull.size() - 1, threads = atoi(argv[2]);
    
    init(ConvexHull);

    HybridPolyDecomp(threads);

//	if(b) MPI_Abort(MPI_COMM_WORLD, 0);

    clean();
    
    MPI_Finalize();
    return 0;
}

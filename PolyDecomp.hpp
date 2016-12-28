//
//  PolyDecomp.hpp
//  CMPS 299
//
//  Created by Ali Jbaily on 3/8/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#ifndef PolyDecomp_hpp
#define PolyDecomp_hpp

#include <iostream>
#include <set>
#include <map>
#include "convex_hull.hpp"
#include <math.h>
#include <algorithm>
#include <vector>
#include <climits>
#include <string.h>
#include "mpi.h"
#include "omp.h"
#include "priority_update.hpp"
#include "ArrayX.cpp"

struct edge {
    int n;
    Point e;
};

struct PointAndInt{
    Point p;
    int i;
};

typedef struct {
    unsigned int u;
    set<pair<int, int>> S;
} IntAndSet;

typedef struct {
    unsigned int u;
    ArrayX<pair<int, int>> a;
} IntAndX;

typedef struct {
    IntAndX** A;
    unsigned int u;
} Alg2Type;

typedef vector<Point> Polygon;

int gcd (int, int);
bool PolyDecomp();
bool PolyDecompArray();
void PolyDecompDFS(int, Point);
bool PolyDecompDFS();
bool MultiThreadedPolyDecompDFS(int threads);
bool HybridPolyDecomp(int);
Alg2Type HybridPolyDecompNum(int threads);
Alg2Type HybridPolyDecompNumX(int threads);
Point** CalcIP(vector<Point>);
bool belongToIP(Point*);
vector<PointAndInt> getLayer(int);
void init(Polygon p);
void clean();
void printPoint(Point *);

#endif /* PolyDecomp_hpp */
//
//  PolyDecomp.cpp
//  CMPS 299
//
//  Created by Ali Jbaily on 3/8/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#include "PolyDecomp.hpp"

#define CHUNKSIZE 1000

Point** IP, *v0;
vector<Point> latticePoints;
Polygon ConvexHull;

int ymin = INT32_MAX, ymax = 0, xmin = INT32_MAX, xmax = 0, m, world_rank, world_size, threads, atomicLevel;
edge* edge_sequence;
atomic<bool> B;
map<Point, atomic<int>> flags;

IntAndX** newA;

void printPoint(const Point *p){
    cout << p->x << ", " << p->y << endl;
}

vector<PointAndInt> getLayer(int threads) {
    
    vector<PointAndInt> oldSet, newSet;
    Point p;
    map<Point, int> mp;
    oldSet.push_back({*v0, 0});
    newSet.push_back({*v0, 0});
    threads--;
    
    
    for (int i = 0; i < m - 1; i++) {
        for (auto u : oldSet) {
            int step;
            if (edge_sequence[i].n == 1 || threads >= (edge_sequence[i].n - 1)) step = 1;
            else step = floor((edge_sequence[i].n - 1)/threads);
            for (int k = step; k <= edge_sequence[i].n; k+=step) {
                p = u.p + edge_sequence[i].e * k;
                if (belongToIP(&p) && mp[{p.x,p.y}] == 0) {
                    flags[p] = i;
                    newSet.push_back({p, i + 1});
                    threads--;
                    if (threads == 0) return newSet;
                }
            }
        }
        
        oldSet = newSet;
    }
    
    return newSet;
}

void applyFlags(int limit) {
    if (limit > m) limit = m;
    set<Point> oldSet, newSet;
    Point p;
    oldSet.insert(*v0);
    newSet.insert(*v0);
    flags[*v0] = -1;
    
    for (int i = 0; i < limit; i++) {
        for (auto u : oldSet) {
            for (int k = 1; k <= edge_sequence[i].n; k++) {
                p = u + edge_sequence[i].e * k;
                if (belongToIP(&p)) {
                    newSet.insert(p);
                    flags[p] = -1;
                }
            }
        }
        oldSet = newSet;
    }
    atomicLevel = limit - 1;
}

void PolyDecompDFS(int i, Point *v, int stop) {
    if (B) return;
    int end;
    i < m - 1 ? end = edge_sequence[i].n : end = edge_sequence[i].n - 1;
    for (int k = 0; k <= end && !B; k++) {
        Point p = *v + edge_sequence[i].e*k;
        if (p == *v0 && k > 0) {B = true; return;}
        if (!belongToIP(&p)) return;
        else if (k > 0 && i + 1 < m)
            PolyDecompDFS(i+1, &p, m);
        else if(k == 0 && i + 1 < stop)
            PolyDecompDFS(i+1, &p, stop);
    }
}

void PolyDecompDFSAtomic(int i, Point *v, int stop) {
    if (B) return;
    int end;
    i < m - 1 ? end = edge_sequence[i].n : end = edge_sequence[i].n - 1;
    for (int k = 0; k <= end && !B; k++) {
        Point p = *v + edge_sequence[i].e*k;
        if (p == *v0 && k > 0) {B = true;return;}
        if (!belongToIP(&p)) return;
        else {
            if(flags[p] == -1) {
                priority_update(&flags[p], i);
                if (i + 1 <= atomicLevel)
                    PolyDecompDFSAtomic(i+1, &p, m);
                else if (i + 1 < m)
                    PolyDecompDFS(i+1, &p, m);
            }
            else if (flags[p] > i){
                
                if(k > 0) {
                    if (flags[p] < stop) {
                        stop = flags[p];
                        atomic_compare_exchange_strong(&flags[p], &stop, i);
                    }
                    if (i + 1 <= atomicLevel)
                        PolyDecompDFSAtomic(i+1, &p, stop);
                    else if (i + 1 < m)
                        PolyDecompDFSAtomic(i+1, &p, stop);
                }
                
                else if(k == 0 && i + 1 < stop) {
                    if (i + 1 <= atomicLevel)
                        PolyDecompDFSAtomic(i+1, &p, stop);
                    else if (i + 1 < m)
                        PolyDecompDFS(i+1, &p, stop);
                }
            }
            else if (k == 0 && i + 1 < stop) { // k == 0 so it can continue progressing from that point and not fail when i > flag[p]
                if (i + 1 <= atomicLevel)
                    PolyDecompDFSAtomic(i+1, &p, m);
                else if (i + 1 < m)
                    PolyDecompDFS(i+1, &p, m);
            }
        }
    }
}

bool HybridPolyDecomp(){
    applyFlags(4);
    double start, end, totTime;
    int tid;
    omp_set_num_threads(threads);
    start = omp_get_wtime();
    Point p;
#pragma omp parallel shared(totTime) private(tid, p)
    {
        tid = omp_get_thread_num();
        double tstart, tend;
        tstart = omp_get_wtime();
        
        int start, end, length, startIndex, endIndex;
        start = world_rank*(edge_sequence[0].n/world_size);
        world_rank == world_size - 1 ? end = edge_sequence[0].n : end = ((world_rank+1)*(edge_sequence[0].n/world_size)) - 1;
        
        length = end - start + 1;
    
        startIndex = tid*(length/threads);
        startIndex += start;
        tid == threads -1 ? endIndex = length - 1: endIndex = ((tid+1)*(length/threads)) - 1;
        endIndex += start;
        for (int k = startIndex; k <= endIndex && !B; k++) {
            p = *v0 + edge_sequence[0].e * k;
            priority_update(&flags[p], 1);
            PolyDecompDFSAtomic(1, &p, m);
        }
        tend = omp_get_wtime();
        totTime += (tend - tstart);

    }
    end = omp_get_wtime();
    return B;
}

bool HybridPolyDecomp(int t){
    threads = t;
    double start, end;
	start = MPI_Wtime();
    B = HybridPolyDecomp();
    end = MPI_Wtime();
	cout << "rank " << world_rank << " took " << end - start << " seconds with result " << boolalpha << B << endl;
	return B;
}

bool between(Point *p, Point *p1, Point *p2) {
    if(p1->y < p->y && p->y < p2->y) // between
        return true;
    else if((p1->x <= p->x && p1->y == p->y) || (p->x <= p2->x && p->y == p2->y)) // extremities
        return true;
    else if(p1->x <= p->x && p->x <= p2->x && p1->y == p2->y) // on the same line
        return true;
    else
        return false;
}

int getSector(Point *p){
    for (int r = 0; r < world_size; r++) {
        int start = r*((xmax-xmin)/world_size) + xmin, end;
        r == world_size - 1 ? end = xmax : end = ((r+1)*((xmax-xmin)/world_size)) - 1 + xmin;
        if (p->x >= start && p->x <= end)
            return r;
    }
    cout << "-1\n";
    return -1;
}

void getLatticePoints() {
    int start = world_rank*((xmax-xmin)/world_size), end;
    world_rank == world_size - 1 ? end = xmax-xmin: end = ((world_rank+1)*((xmax-xmin)/world_size)) - 1;
    if (end < start) end = start;
    for (coord_t y = ymin; y <= ymax; y++)
        for (coord_t x = xmin + start; x <= xmin + end; x++) {
            Point v = {x, y};
            if (belongToIP(&v))
                latticePoints.push_back(v);
        }
}

Alg2Type PolyDecompNum(int threads) {
    
    omp_set_num_threads(threads);
    
    int start = world_rank*((xmax-xmin)/world_size) + xmin, end;
    world_rank == world_size - 1 ? end = xmax : end = ((world_rank+1)*((xmax-xmin)/world_size)) - 1 + xmin;
    if (end < start) end = start;
    
    int length = end - start + 1;
    
    newA = new IntAndX*[ymax+1];
    IntAndX** oldA = new IntAndX*[ymax+1];
    
    mutex bmx;
    vector<int> *buffer = new vector<int>[world_size];
    
    int rank = getSector(v0);
    
    
#pragma omp parallel shared(bmx, buffer, newA, oldA)
    {
        int tid = omp_get_thread_num();
        
#pragma omp for
        for (int y = 0; y <= ymax; y++) {
            newA[y] = new IntAndX[end-start+1];
            oldA[y] = new IntAndX[end-start+1];
        }
        
        int tstart = tid*(length/threads) + start, tend;
        tid == threads - 1 ? tend = end : tend = ((tid+1)*(length/threads)) - 1 + start;
        if (tend < tstart) tend = tstart;
        
#pragma omp single
        {
            if(world_rank == rank) {
                oldA[(int)v0->y][((int)v0->x)-start].u = 1;
                newA[(int)v0->y][((int)v0->x)-start].u = 1;
            }
        }
        
        for (int i = 0; i < m; i++) {
            for (int y = ymin; y <= ymax; y++)
                for (int x = tstart; x <= tend; x++) {
                    Point v = {x,y};
                    if (oldA[y][x-start].u > 0) {
                        for (int k = 1; k <= edge_sequence[i].n; k++) {
                            Point vp = v + edge_sequence[i].e*k;
                            if (belongToIP(&vp)) {
                                int xp = vp.x, yp = vp.y, sector = getSector(&vp);
                                if (sector == world_rank) {
#pragma omp atomic
                                    newA[yp][xp-start].u += oldA[y][x-start].u;
                                    newA[yp][xp-start].a.atomicAdd({k,i});
                                }
                                else {
                                    int elems[5] = {k,i,oldA[y][x-start].u, xp,yp};
                                    bmx.lock();
                                    buffer[sector].insert(buffer[sector].end(), elems, elems + 5);
                                    bmx.unlock();
                                }
                                
                            }
                        }
                    }
                }
            
#pragma omp barrier
            
#pragma omp single
            {
                /*****************      Send        *****************/
                for (int j = 0; j < world_size; j++) {
                    if(j != world_rank) {
                        int count = buffer[j].size();
                        MPI_Send(&count, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                        
                        for (int chunk = 0; chunk < count/CHUNKSIZE; chunk++)
                            MPI_Send(&buffer[j][chunk*CHUNKSIZE], CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD);
                        
                        if(count%CHUNKSIZE > 0)
                            MPI_Send(&buffer[j][(count/CHUNKSIZE)*CHUNKSIZE], count%CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD);
                    }
                }
                
                /****************************************************/
                
                
                /*****************       Receive        *****************/
                for (int j = 0; j < world_size; j++) {
                    if(j != world_rank) {
                        
                        int count = -1, *array;
                        MPI_Recv(&count, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        array = new int[count];
                        
                        for (int chunk = 0; chunk < count/CHUNKSIZE; chunk++)
                            MPI_Recv(&array[chunk*CHUNKSIZE], CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        
                        if(count%CHUNKSIZE > 0)
                            MPI_Recv(&array[(count/CHUNKSIZE)*CHUNKSIZE], count%CHUNKSIZE, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        
                        /******** After Receiving ********/
                        
                        for (int k = 0; k < count; k+=5) {
                            int x = array[k+3], y = array[k+4];
                            newA[y][x-start].u += array[k+2];
                            newA[y][x-start].a.add({array[k], array[k+1]});
                        }
                        delete array;
                    }
                    /*********************************/
                }
                /*******************************************************/
                
                for (int sec = 0; sec < world_size; sec++)
                    if(sec != world_rank)
                        buffer[sec].erase(buffer[sec].begin(), buffer[sec].end());
            }
            
            for (int y = ymin; y <= ymax; y++)
                for (int x = tstart; x <= tend; x++)
                    oldA[y][x-start].u = newA[y][x-start].u;
            
//#pragma omp single
            MPI_Barrier(MPI_COMM_WORLD);
            
        }
    }
    
    //    for (int i = 0; i <= ymax; i++)
    //        delete [] oldA[i];
    delete oldA;
    
    //    cout << newA[(int)v0->y][(int)v0->x].a.size << " elements at v0" << endl;
    
    return {newA, newA[(int)v0->y][(int)v0->x - start].u};
}

Alg2Type HybridPolyDecompNum(int threads) {
    double start, end;
    start = MPI_Wtime();
    Alg2Type ret = PolyDecompNum(threads);
    end = MPI_Wtime();
    
    int rank = getSector(v0);
    
    if(world_rank == rank)
        cout << "p" << world_rank << " finished in " << end - start << " with " << ret.u << " summands" << endl;
    else
        cout << "p" << world_rank << " finished in " << end - start << endl;
    
    return ret;
}

Point** CalcIP() {
    int convexSize = (int)ConvexHull.size();
    ymax = ymin = ConvexHull[0].y;
    for (int i = 1; i < convexSize; i++) {
        if (ConvexHull[i].y > ymax)
            ymax = ConvexHull[i].y;
        if (ConvexHull[i].y < ymin)
            ymin = ConvexHull[i].y;
        if (ConvexHull[i].x > xmax)
            xmax = ConvexHull[i].x;
        if (ConvexHull[i].x < xmin)
            xmin = ConvexHull[i].x;
    }
    
    IP = new Point*[ymax + 1];
    for (int i = 0; i <= ymax; i++)
        IP[i] = new Point[2];
    
    double slope[convexSize];
    for (int i = 0; i < convexSize - 1; i++)
        slope[i]=(ConvexHull[i].y-ConvexHull[i+1].y)/(ConvexHull[i].x - ConvexHull[i+1].x);
    
    for (intmax_t y = ymin; y <= ymax; y++) {
        int count = 0;
        double x[2] = {-1, -1};
        for (int p = 0; p < convexSize && count < 2; p++) {
            const int y1 = ConvexHull[p].y;
            const int y2 = ConvexHull[p+1].y;
            if (y >= min(y1, y2) && y <= max(y1, y2)) {
                if (abs(slope[p]) < 1e-2) {
                    x[0] = ConvexHull[p].x;
                    x[1] = ConvexHull[(p+1)%convexSize].x;
                    count = 2;
                }
                else {
                    double value = ((y-ConvexHull[p].y)/slope[p])+ConvexHull[p].x;
                    if (x[0] != value)
                        x[count++] = value;
                }
            }
        }
        const int x1 = (const int)x[0], x2 = (const int)x[1];
        IP[y][0].x = ceil(min(x1, x2));
        IP[y][1].x = floor(max(x1, x2));
        IP[y][0].y = IP[y][1].y = y;
    }
    
    if (IP[ymin][0].x == -1)
        IP[ymin][0].x = IP[ymin][1].x;
    if (IP[ymax][0].x == -1)
        IP[ymax][0].x = IP[ymax][1].x;
    
    return IP;
}

bool belongToIP(Point* p) {
    return p->y <= ymax && p->y >= ymin && p->x >= IP[(int)p->y][0].x && p->x <= IP[(int)p->y][1].x;
}

int gcd (int a, int b) {
    int c;
    while (a != 0) {
        c = a; a = b%a;  b = c;
    }
    return b;
}

void init(Polygon p) {
    ConvexHull = p;
    m = (int)ConvexHull.size() - 1;
    edge_sequence = new edge[m];
    
    for (int i = 0; i < m; i++) {
        coord_t dx = ConvexHull[i+1].x - ConvexHull[i].x, dy = ConvexHull[i+1].y - ConvexHull[i].y;
        int n = gcd(dx, dy);
        n > 0 ? edge_sequence[i] = {n, dx/n, dy/n} : edge_sequence[i] = {-n, dx/-n, dy/-n};
    }
    
    v0 = &ConvexHull[0];
    CalcIP();
    
    MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &world_size);
    
}

void clean() {
    delete [] edge_sequence;
    for (int i = 0; i <= ymax; i++)
        delete [] IP[i];
    delete [] IP;
    delete newA;
}

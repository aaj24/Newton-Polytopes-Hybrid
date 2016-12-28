//
//  convex_hull.hpp
//  CMPS 299
//
//  Created by Ali Jbaily on 2/13/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#ifndef convex_hull_hpp
#define convex_hull_hpp

#include <stdio.h>
#include <algorithm>
#include <vector>

using namespace std;

typedef double coord_t;         // coordinate type
typedef double coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2

struct Point {
    coord_t x, y;
    
    bool operator <(const Point &p) const {
        return x < p.x || (x == p.x && y < p.y);
    }
    
    Point operator +(const Point &p) const {
        return Point {x + p.x, y + p.y};
    }
    
    Point operator -(const Point &p) const {
        return Point {x - p.x, y - p.y};
    }
    
    Point operator *(const int m) const {
        return Point {x * m, y * m};
    }
    bool operator ==(const Point &p) const {
        return x == p.x && y == p.y;
    }
    bool operator !=(const Point &p) const {
        return !(*this == p);
    }
};

coord2_t cross(const Point &, const Point &, const Point &);

vector<Point> convex_hull(vector<Point>);

#endif /* convex_hull_hpp */

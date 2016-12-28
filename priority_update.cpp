//
//  priority_update.cpp
//  CMPS 299
//
//  Created by Ali Jbaily on 4/13/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#include "priority_update.hpp"

using namespace std;

void priority_update(atomic<int>* addr, int newval){
    int oldval = *addr;
    while (newval > oldval)
        if (atomic_compare_exchange_strong(addr, &oldval, newval))
            return;
        else
            oldval = *addr;
}
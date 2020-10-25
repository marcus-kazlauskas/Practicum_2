//
//  main.cpp
//  Program_2 Работа №5, Вариант №2, Задание №7
//
//  Created by Kozlov Mikhail on 04.05.16.
//  Copyright © 2016 Kozlov Mikhail. All rights reserved.
//

#include <iostream>
#include <iomanip>

#include "functions.h"

using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    
    int elCount = 10;
    double xInit = 0;
    double xEnd = 1;
    double tInit = 0;
    double tEnd = 1;
    double e = 0.01;
    double hx = (xEnd-xInit)/(elCount);
    double hxConst = hx;
    double ht = (tEnd-tInit)/(elCount);
//    double htConst = ht;
    int numh = 1; //full elCounth = elCount*numh
    int numt = 1; //full elCountt = elCount*numt
    double u[elCount+1];
    double *u1 = new double[elCount*numh*2+1];
    double *u1Short = new double[elCount*numh+1];
    double *u2 = new double[elCount*numh+1];
//    double u1h[elCount+1];
//    double u2h[elCount+1];
    
    for (int i = 0; i <= elCount; i++){
        u[i] = sol(xInit+hx*i, tEnd);
        u2[i] = 0;
        u1Short[i] = 1;
    }
    
//    cout << psi(tInit+ht*1)+hx*(cos(0)-2*psi1(tInit+ht*1))+pow(ht, 2)/2*(psi2(tInit+ht*1)) << endl;
    
    while (norm1(u2, u1Short, elCount, numh)/(pow(2, 2)-1) > e){
        hx /= 2;
        ht /= 2;
        numh *= 2;
        numt *= 2;

        delete [] u1;
        delete [] u1Short;
        delete [] u2;
        
        u1 = new double[elCount*numh*2+1];
        u1Short = new double[elCount*numh+1];
        u2 = new double[elCount*numh+1];
        
        cout << numh*2 << " " << hx/2 << endl;
        
        for (int l = 0; l <= elCount*numh; l++){
            u2[l] = fi(xInit+hx*l);
        }
        
        for (int l = 0; l <= elCount*numh*2; l++){
            u1[l] = fi(xInit+hx/2*l);
        }

        difSch(u2, hx, ht, numh, numt, xInit, tInit, elCount);
        difSch(u1, hx/2, ht/2, numh*2, numt*2, xInit, tInit, elCount);
        
        for (int i = 0; i <= elCount*numh; i++){
            u1Short[i] = u1[i*2];
        }
        
//        for (int i = 0; i <= elCount*num; i++){
//            u1h[i] = u1Short[i*num];
//            u2h[i] = u2[i*num];
//        }
    }
    
    cout << "Solution of the equation with error = " << e << ":" << endl;
    
    for (int i = 0; i <= elCount; i++){
        cout << xInit+hxConst*i << " " << u[i] << " " << u1Short[i*numh] << " " << fabs(u[i]-u2[i*numh]) << " " << fabs(u[i]-u1Short[i*numh]) <<  endl;
    }
    
    return 0;
}

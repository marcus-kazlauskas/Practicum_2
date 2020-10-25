//
//  functions.h
//  Program_2
//
//  Created by Kozlov Mikhail on 04.05.16.
//  Copyright © 2016 Kozlov Mikhail. All rights reserved.
//

#ifndef functions_h
#define functions_h

#include <stdlib.h>
#include <math.h>

double norm1(double *y2, double *y1, int elCount, int num){
    double max = fabs(y2[0]-y1[0]);
    
    for (int i = 1; i <= elCount*num; i++){
        if (fabs(y2[i]-y1[i]) > max){
            max = fabs(y2[i]-y1[i]);
        }
    }
    
    return max;
}

double norm2(double *y2, double *y1){
    double sum = 0;
    
    for (int i = 1; i <= 10; i++){
        sum += fabs(y2[i]-y1[i]);
    }
    
    return sum;
}

double norm3(double *y2, double *y1){
    double sum = 0;
    
    for (int i = 1; i <= 10; i++){
        sum += pow(y2[i]-y1[i], 2);
    }
    
    return pow(sum, 0.5);
}

inline double a(double x, double t){
    return 1;
}

inline double b(double x, double t){
    return cos(x);
}

inline double fi(double x){
    return log(1+pow(x, 2))+sin(x);
}

inline double psi(double t){
    return log(1+pow(t, 2));
}

inline double psi1(double t){
    return 2*t/(1+pow(t, 2));
}

inline double psi2(double t){
    return 2/(1+pow(t, 2))-pow(2*t, 2)/pow(1+pow(t, 2), 2);
}

inline double sol(double x, double t){
    return log(1+pow(x-t, 2))+sin(x);
}

void difSch(double *u, double hx, double ht, int numh, int numt, double xInit, double tInit, int elCount){
    double uBuf[elCount*numh+1];
    
    for (int l = 0; l <= elCount*numh; l++){
        u[l] = fi(xInit+hx*l);
    }
    
    for (int n = 1; n <= elCount*numt; n++){
        for (int l = 2; l <= elCount*numh; l++){
            uBuf[l] = u[l]+ht/(2*hx)*(-u[l-2]+4*u[l-1]-3*u[l])+
                      pow(ht, 2)/(2*pow(hx, 2))*(u[l-2]-2*u[l-1]+u[l])+ht*cos(xInit+hx*l)+pow(ht, 2)/2*sin(xInit+hx*l);
        }
        
        u[0] = psi(tInit+ht*n);
//        u[1] = u[0];
        u[1] = psi(tInit+ht*n)+hx*(cos(0)-2*psi1(tInit+ht*n))+pow(hx, 2)/2*(psi2(tInit+ht*n)+0-0-sin(0));
        
        for(int l = 2; l <= elCount*numh; l++){
            u[l] = uBuf[l];
        }
    }
}

#endif /* functions_h */

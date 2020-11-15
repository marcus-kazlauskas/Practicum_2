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

double norm1(double *y2, double *y1, int elCount, int num){     // p-норма при p=infinity, используется для определения точности численного решения
    double max = fabs(y2[0]-y1[0]);
    
    for (int i = 1; i <= elCount*num; i++){
        if (fabs(y2[i]-y1[i]) > max){
            max = fabs(y2[i]-y1[i]);
        }
    }
    
    return max;
}

double norm2(double *y2, double *y1){                           // манхэттенское расстояние, в задаче не используется
    double sum = 0;
    
    for (int i = 1; i <= 10; i++){
        sum += fabs(y2[i]-y1[i]);
    }
    
    return sum;
}

double norm3(double *y2, double *y1){                           // евклидова норма, в задаче не используется
    double sum = 0;
    
    for (int i = 1; i <= 10; i++){
        sum += pow(y2[i]-y1[i], 2);
    }
    
    return pow(sum, 0.5);
}

inline double a(double x, double t){                            // a(x,t)=1 >= 0
    return 1;
}

inline double b(double x, double t){
    return cos(x);
}

inline double fi(double x){                                     // u(x,0)=fi(x); 0 <= x <= 1
    return log(1+pow(x, 2))+sin(x);
}

inline double psi(double t){                                    // u(0,t)=psi(t); 0 <= t <= 1
    return log(1+pow(t, 2));
}

inline double psi1(double t){                                   // psi'(t); u'(0,t)=psi'(t)
    return 2*t/(1+pow(t, 2));
}

inline double psi2(double t){                                   // psi''(t); u''(0,t)=psi''(t)-b(0,t)+psi'(t)-b'_t(0,t)+b'_x(0,t)
    return 2/(1+pow(t, 2))-pow(2*t, 2)/pow(1+pow(t, 2), 2);
}

inline double sol(double x, double t){                          // аналитическое решение u(x,t)
    return log(1+pow(x-t, 2))+sin(x);
}

void difSch(double *u, double hx, double ht, int numh, int numt, double xInit, double tInit, int elCount){              // вычисление  вычисление сеточной ф-ции с помощью разностной схемы
    double uBuf[elCount*numh+1];
    
    for (int l = 0; l <= elCount*numh; l++){
        u[l] = fi(xInit+hx*l);
    }
    
    for (int n = 1; n <= elCount*numt; n++){
        for (int l = 2; l <= elCount*numh; l++){
            uBuf[l] = u[l]+ht/(2*hx)*(-u[l-2]+4*u[l-1]-3*u[l])+                                                         // вычисление следа ф-ции через предложенную разностную схему
                      pow(ht, 2)/(2*pow(hx, 2))*(u[l-2]-2*u[l-1]+u[l])+ht*cos(xInit+hx*l)+pow(ht, 2)/2*sin(xInit+hx*l);
        }
        
        u[0] = psi(tInit+ht*n);                                                                                         // u(0,t)=psi(t)
//        u[1] = u[0];                                                                                                  // аппроксимация ф-ции вблизи граничного условия через разложение Тейлора
        u[1] = psi(tInit+ht*n)+hx*(cos(0)-psi1(tInit+ht*n))+pow(hx, 2)/2*(psi2(tInit+ht*n)+0-0-sin(0));                 // u(0+hx,t)=u(0,t)+u'_x(0,t)*hx+u''_xx(0,t)*hx^2/2+O(hx^3)
                                                                                                                        // которое преобразовано с помощью данной разностной схемы
        for(int l = 2; l <= elCount*numh; l++){                                                                         // u(0+hx,t)=psi(t)+hx*(b(0)-psi'(t))+h^2/2*(psi''(t)+0-0+b'(0))
            u[l] = uBuf[l];
        }
    }
}

#endif /* functions_h */

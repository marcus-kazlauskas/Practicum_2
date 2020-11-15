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
    
    int elCount = 10;                                   // количество точек разбиения отрезка
    double xInit = 0;                                   // начало отрезка
    double xEnd = 1;                                    // конец отрезка
    double tInit = 0;                                   // начальное время
    double tEnd = 1;                                    // конечное время
    double e = 0.001;                                   // погрешность
    double hx = (xEnd-xInit)/(elCount);                 // шаг вычислений по x
    double hxConst = hx;
    double ht = (tEnd-tInit)/(elCount);                 // шаг вычислений по t; при ht = hx выполняется признак устойчивочти ht <= 2*hx
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
        u[i] = sol(xInit+hx*i, tEnd);                   // след аналитического решения
        u2[i] = 0;
        u1Short[i] = 1;
    }
    
//    cout << psi(tInit+ht*1)+hx*(cos(0)-2*psi1(tInit+ht*1))+pow(ht, 2)/2*(psi2(tInit+ht*1)) << endl;
    
    while (norm1(u2, u1Short, elCount, numh)/(pow(2, 2)-1) > e){        // численное решение, норма делится на 2^2-1,
        hx /= 2;                                                        // потому что сетка удваивается каждую итерацию,
        ht /= 2;                                                        // а разностная схема обеспечивает 2-й порядок аппроксимации
        numh *= 2;
        numt *= 2;

        delete [] u1;
        delete [] u1Short;
        delete [] u2;
        
        u1 = new double[elCount*numh*2+1];
        u1Short = new double[elCount*numh+1];
        u2 = new double[elCount*numh+1];
        
        cout << numh*2 << " " << hx/2 << endl;                          // величина дополнительного разбиения отрезка (помимо изначальный 10 частей)
        
        for (int l = 0; l <= elCount*numh; l++){                        // след ф-ции при t=0, u(x,0)=fi(x), 0 <= x <= 1
            u2[l] = fi(xInit+hx*l);
        }
        
        for (int l = 0; l <= elCount*numh*2; l++){                      // след ф-ции при t=0 на удвоенной сетке
            u1[l] = fi(xInit+hx/2*l);
        }

        difSch(u2, hx, ht, numh, numt, xInit, tInit, elCount);          // вычисление следа u(x,t) с помощью разностной схемы, аппроксимирующей задачу
        difSch(u1, hx/2, ht/2, numh*2, numt*2, xInit, tInit, elCount);  // то же самое на удвоенной одновременно и по x и по t сетке
        
        for (int i = 0; i <= elCount*numh; i++){
            u1Short[i] = u1[i*2];                                       // каждое второе значение ф-ции выкидывается, чтобы потом считать погрешность
        }
        
//        for (int i = 0; i <= elCount*num; i++){
//            u1h[i] = u1Short[i*num];
//            u2h[i] = u2[i*num];
//        }
    }
    
    cout << "Solution of the equation with error = " << e << ":" << endl;
    
    for (int i = 0; i <= elCount; i++){
        cout << xInit+hxConst*i << " " << u[i] << " " << u1Short[i*numh] << " " << fabs(u[i]-u2[i*numh]) << " " << fabs(u[i]-u1Short[i*numh]) <<  endl;
    }   // вывод: x; аналитическое решение u(x,1); численное решение u(x,1); погрешность предпосленей итерации; погрешность итогового решения
    
    return 0;
}

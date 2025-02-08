#pragma once
#include<iostream>
#include<cstdio> 
#include<math.h>
#include"Read.h"
using namespace std;

class Calc{
    public: 
    static void capway(double* m, int n, double dr, double* y);
    static double spline(int n,double dr,double* m, double* y,double t);
    static double dist(point p, double x, double y, double z);
    static double dist(point p, point q);
    static double spline1(int n, double dr, double* y, double t);
};

#pragma once
#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
#include<map>
#include<algorithm>
using namespace std;

struct point{
    double x;
    double y;
    double z;
};

class MyRead{
    public:
    int point_num = 0;
    point* point_vec = NULL;
    double lx, ly, lz;
    int nx, ny, nz;
    double* V = NULL;
    double cutoff;
    double dr;
    int mesh;
    int l;
    double* distr = NULL;
    map<string, string> arg2meth;
    
    void ReadInput(string pathname);
    void ReadPoints();
    void ReadV();
    void ReadDistr();
};
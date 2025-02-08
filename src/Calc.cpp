#include "Calc.h"
using namespace std;

void Calc::capway(double* m, int n, double dr, double* y){
	double a = dr;
	double b = 4*dr;
	double c[n];
	c[0] = 0;
	m[0] = 0;
	for(int i = 0;i < n-1;++i){
		c[i] = dr;
		m[i] = (y[i+1]-y[i-1])/dr;
	}
	c[n-1] = -1;
	m[n-1] = 0;

	c[0] = c[0] / b;
	m[0] = m[0] / b;

	for(int i = 1;i < n;++i){
		double la = 1.0/(b-a*c[i-1]);
		c[i] = c[i]*la;
		m[i] = (m[i]-a*m[i-1])*la;
	}
	for(int i = n-2;i--;i>0){
		m[i] = m[i]-c[i]*m[i+1];
	}

}

double Calc::spline(int n,double dr,double* m, double* y,double t){
	int index = int(t/dr);
	if(index >= n-1){
		return 0.0;
	}
	double a = y[index];
	double b = (y[index+1]-y[index])/dr-dr*m[index]/2.0-dr*(m[index+1]-m[index])/6.0;
	double c = m[index]/2.0;
	double d = (m[index+1]-m[index])/(6.0*dr);
	double la = t-index*dr;
	return (a+b*la+c*la*la+d*la*la*la);
}

double Calc::dist(point p, double x, double y, double z){
	return sqrt((p.x-x)*(p.x-x)+(p.y-y)*(p.y-y)+(p.z-z)*(p.z-z));
}

double Calc::dist(point p, point q){
	return sqrt((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y)+(p.z-q.z)*(p.z-q.z));
}

double Calc::spline1(int n, double dr, double* y, double t){
	int index = int (t / dr);
	if(index >= 80)//!!!
		return 0;
	double la = (t-index*dr)/dr;
	return ((1.0-la)*y[index] + la*y[index+1]);
}

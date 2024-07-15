//============================================================================
// Name        : helloworld.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//
//	Riddhi Singh, May, 2014
//  The Pennsylvania State University
//  rus197@psu.edu
//
//  Adapted by Tori Ward, July 2014
//  Cornell University
//  vlw27@cornell.edu
//
//  Adapted by Jonathan Herman and David Hadka, Sept-Dec 2014
//  Cornell University and The Pennsylvania State University
//
//  Adapted by Julianne Quinn, July 2015 as DPS problem
//  Cornell University
//  jdq8@cornell.edu
//============================================================================

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sstream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/tools/roots.hpp>
#include "./../borg/moeaframework.h"
#include "./../boostutil.h"
#include "./../borg/borgms.h"

#define LL 12
#define nSamples 100

double inf_pool [20000][LL];
double lat_pool [20000][LL];

namespace ublas = boost::numeric::ublas;
namespace tools = boost::math::tools;
using namespace std;

int na = 5;
int nb = 5;
int nvars = na*(2*nb+2);
int nobjs = 4;
int nconstrs = 0;

ublas::vector<double> tot_e(nSamples);
ublas::vector<double> min_e(nSamples);
ublas::vector<double> wd(nSamples);
ublas::vector<double> sp(nSamples);

ublas::vector<double> C(na*nb);
ublas::vector<double> R(na*nb);
ublas::vector<double> W(na*2);
ublas::vector<double> W1(na);
ublas::vector<double> W2(na);

ublas::vector<double> xvar(nb);

ublas::vector<double> s1(LL+1);
ublas::vector<double> s2(LL+1);

// some parameters
double b1 = 18.5;
double b2 = 0.2265;
double db1 = 47.4;
double db2 = 0.07;
double smin1 = 50e8;
double smax1 = 247e8;
double smin2 = 6e8;
double smax2 = 57e8;
double kkhead[10] = {205,14,122,12.5,18.7,99.3,16,16,73,9.2};
double kkcsa[10] = {0.0119,0.0151,0.1062,0.1112,0.1310,0.2426,0.2651,0.2849,0.3047,0.3262};
double rmin1[12] = {313,177,207,369,365,367,435,258,340,347,419,353};
double rmax1[12] = {639,596,631,888,927,867,945,1018,714,700,631,594};
double rmin2[12] = {293,259,233,350,589,599,472,461,494,646,587,402};
double rmax2[12] = {513,422,662,1168,1122,1171,1109,1079,863,1203,866,493};
double demand[12] = {550,500,350,750,1100,900,800,750,750,800,650,600};
double inf_mm[12] = {5.2322, 5.3028, 5.5483, 5.8276, 6.2475, 6.7888, 7.1154, 6.8763, 6.9024, 6.7491, 6.0982, 5.4033};
double std_mm[12] = {0.3484, 0.4975, 0.4516, 0.2204, 0.3439, 0.3519, 0.4682, 0.4025, 0.4327, 0.3932, 0.2920, 0.2380};

double RBFpolicy(ublas::vector<double> xvar,ublas::vector<double> C, ublas::vector<double> R, ublas::vector<double> W );
double inhead1(double ss);
double inhead2(double ss);
double outhead1(double rr);
double outhead2(double rr);

void headwater_model(double* vars, double* objs, double* constrs)
{
	// read RBF parameters
	int mark = 0;
	for (int i = 0; i < na; i++){
		for (int j = 0; j < nb; j++ ){
			C(i*nb+j) = vars[mark];
			mark++;
		}
		for (int j = 0; j < nb; j++ ){
			R(i*nb+j) = 10*vars[mark];
			mark++;
		}
		for (int j = 0; j < 2; j++ ){
			W(i*2+j) = vars[mark];
			mark++;
		}
	}

	// normalization of W
	double total1,total2;
	for (int i = 0; i < na; i++){
		total1 = total1 + W(i*2);
		total2 = total2 + W(i*2+1);
	}

	if (total1 != 0){
		for (int i = 0; i < na; i++){
		  W1(i) = W(i*2)/total1;
		}
	}
	else {
		for (int i = 0; i < na; i++){
		  W1(i) = 1/na;
		}
	}

	if (total2 != 0){
		for (int i = 0; i < na; i++){
			W2(i) = W(i*2+1)/total2;
			}
		}
	else {
		for (int i = 0; i < na; i++){
			W2(i) = 1/na;
		}
	}

	std::fill(tot_e.begin(), tot_e.end(), 0);
	std::fill(min_e.begin(), min_e.end(), 0);
	std::fill(wd.begin(), wd.end(), 0);
	std::fill(sp.begin(), sp.end(), 0);

	int linesToUse [nSamples];
	double tmpinf, tmplat, tmpd1, tmpr1, tmpa1, tmpd2, tmpr2, tmpa2, tmpe, tmpme;

	srand (time(NULL));
	for (int s=0; s < nSamples; s++) {
		//pick a random number based on time
		//choose 100 of 10,000 available inflow value lines
		linesToUse[s] = rand() % 20000;
	}


	for (int s = 0; s < nSamples; s++) {
		double *in_flow = new double [LL];
		double *lat_flow = new double [LL];
		int index = linesToUse[s];
		for (int i = 0; i < LL; i++) {
			in_flow[i] = inf_pool[index][i];
			lat_flow[i] = lat_pool[index][i];
		}

		s1(0) = 120e8;
		s2(0) = 25e8;

		tmpme = 1e8;

		for (int i = 0; i < LL; i++){
			tmpe = 0;

			tmpinf = exp(in_flow[i]*std_mm[i%12] + inf_mm[i%12]);
			tmplat = b1 + b2*tmpinf + lat_flow[i%12]/4*(db1 + db2*tmpinf);
			tmpinf = tmpinf*30*24*3600;
			tmplat = tmplat*30*24*3600;

			xvar(0) = (s1(i) - smin1)/(smax1 - smin1);
			xvar(1) = (s2(i) - smin2)/(smax2 - smin2);
			xvar(2) = (double)(i%12)/12;
//      xvar(3) = (in_flow[i]>0)*0.33 + 0.34;
//      xvar(4) = (lat_flow[i]>0)*0.33 + 0.34;
			xvar(3) = (in_flow[i]+2.0)/4.0;
			xvar(4) = (lat_flow[i]+2.0)/4.0;


			tmpd1 = RBFpolicy(xvar, C, R, W1)*(rmax1[i%12]-rmin1[i%12]) + rmin1[i%12];
			tmpd1 = tmpd1*30*24*3600;
			tmpa1 = max(0.0, s1(i)+tmpinf-smax1-tmpd1);
			tmpr1 = min(tmpd1+tmpa1, s1(i)+tmpinf-smin1);
			s1(i+1) = s1(i) + tmpinf - tmpr1;

			tmpd2 = RBFpolicy(xvar, C, R, W2)*(rmax2[i%12]-rmin2[i%12]) + rmin2[i%12];
			tmpd2 = tmpd2*30*24*3600;
			tmpa2 = max(0.0, s2(i)+tmplat+tmpr1-smax2-tmpd2);
			tmpr2 = min(tmpd2+tmpa2, s2(i)+tmplat+tmpr1-smin2);
			s2(i+1) = s2(i) + tmpr1 + tmplat - tmpr2;

			tmpe += 8.5*1e3*min(tmpr1,tmpd1)*(inhead1(s1(i))-outhead1(tmpr1));
			tmpe += 8.5*1e3*min(tmpr2,tmpd2)*(inhead2(s2(i))-outhead2(tmpr2));

			for (int j = 0; j < 10; j++){
				tmpe += 8.5*1e3*kkhead[j]*(tmpr1 + tmplat*kkcsa[j]);
			}
			tmpe = tmpe/1e8/1e3/3600;
			tot_e(s) += tmpe;

			tmpme = min(tmpme, tmpe);
			wd(s) += 1/(double)LL*(max(0.0, 1.2*demand[i] - tmpr2/30/24/3600)/1.2/demand[i]);
		}

		min_e(s) = tmpme;
		sp(s) = (s1(LL)+s2(LL))/(smax1+smax2);

		s1.clear();
		s2.clear();
	}
  double tmpnsample = 100.0;
  objs[0] = 0.0;
	objs[1] = 0.0;
	objs[2] = 0.0;
	objs[3] = 0.0;
	for (int s = 0; s < nSamples; s++) {
		objs[0] += -1/tmpnsample*tot_e(s);
		objs[1] += -1/tmpnsample*min_e(s);
		objs[2] += 1/tmpnsample*wd(s);
		objs[3] += -1/tmpnsample*sp(s);
	}
//  constrs[0] = max(0.0, 0.6 - (-1*objs[3]));
	tot_e.clear();
	min_e.clear();
	wd.clear();
	sp.clear();

}

double RBFpolicy(ublas::vector<double> xvar,ublas::vector<double> C, ublas::vector<double> R, ublas::vector<double> W ){
	double Y = 0;
	double tmpy = 0;
	for (int i=0; i<na; i++){
		tmpy = 0;
		for (int j=0; j<nb; j++){
			if (R(i*nb+j) != 0){
				tmpy = tmpy - (xvar(j)-C(i*nb+j))*(xvar(j)-C(i*nb+j))/R(i*nb+j)/R(i*nb+j);
			}
		}
		Y = Y + W(i)*exp(tmpy);
	}
	return Y;
}

double inhead1(double ss){
	double Y;
	ss = ss/100000000;
	Y = -0.0007252553*ss*ss + 0.5690388196*ss + 2502.7122406983;
	return Y;
}

double outhead1(double rr){
	double Y;
	rr = rr/30/3600/24;
	Y = -0.0000002260*rr*rr + 0.0034508966*rr + 2449.9359789846;
	return Y;
}

double inhead2(double ss){
	double Y;
	ss = ss/100000000;
	Y = -0.0217686085*ss*ss + 2.1651210663*ss + 1681.8568952269;
	return Y;
}

double outhead2(double rr){
	double Y;
  rr = rr/30/3600/24;
	Y = -0.0000002574*rr*rr + 0.0022264693*rr + 1619.5630307963;
	return Y;
}

int main(int argc, char* argv[])
{

	for (int i=0;i<20000;i++){
	    for (int j=0;j<LL;j++){
	      inf_pool[i][j] = 0.0;
	    }
	  }

	FILE * myfile;
	myfile = fopen("./../inflow.txt","r");

	int linenum = 0;
	int maxSize = 5000;

	if (myfile==NULL){
	perror("Error opening file");
	} else {
	char buffer [maxSize];
	while (fgets(buffer, maxSize, myfile)!=NULL){

	  linenum++;
	  if (buffer[0]!='#'){
		char *pEnd;
		char *testbuffer = new char [maxSize];
		for (int i=0; i <maxSize; i++){
		  testbuffer[i] = buffer[i];
		}

		for (int cols=0; cols < LL; cols++){
		  inf_pool[linenum-1][cols] = strtod(testbuffer, &pEnd);
		  testbuffer  = pEnd;
		}
	  }

	}
	}

	fclose(myfile);

	for (int i=0;i<20000;i++){
		for (int j=0;j<LL;j++){
		  lat_pool[i][j] = 0.0;
		}
	  }
	FILE * myfile2;
	myfile2 = fopen("./../latflow.txt","r");

	linenum = 0;

	if (myfile2==NULL){
	perror("Error opening file");
	} else {
	char buffer [maxSize];
	while (fgets(buffer, maxSize, myfile2)!=NULL){

	  linenum++;
	  if (buffer[0]!='#'){
		char *pEnd;
		char *testbuffer = new char [maxSize];
		for (int i=0; i <maxSize; i++){
		  testbuffer[i] = buffer[i];
		}

		for (int cols=0; cols < LL; cols++){
		  lat_pool[linenum-1][cols] = strtod(testbuffer, &pEnd);
		  testbuffer  = pEnd;
		}
	  }

	}
	}

	fclose(myfile2);

	unsigned int seed = atoi(argv[1]);
	srand(seed);
	int NFE = atoi(argv[2]);
	// interface with Borg-MS
	BORG_Algorithm_ms_startup(&argc, &argv);
	BORG_Algorithm_ms_max_evaluations(NFE);
	BORG_Algorithm_output_frequency(NFE/100);

	BORG_Problem problem = BORG_Problem_create(nvars, nobjs, nconstrs, headwater_model);

	for (int j=0; j < nvars; j++){
		BORG_Problem_set_bounds(problem, j, 0, 1.0);
	}

	BORG_Problem_set_epsilon(problem, 0, 1);
	BORG_Problem_set_epsilon(problem, 1, 0.1);
	BORG_Problem_set_epsilon(problem, 2, 0.001);
	BORG_Problem_set_epsilon(problem, 3, 0.001);
	char outputFilename[256];
	char runtime[256];
	FILE* outputFile = NULL;
	sprintf(outputFilename, "./sets/testDPS_S%d.set", seed);
	sprintf(runtime, "./runtime/testDPS_S%d.runtime", seed);
	BORG_Algorithm_output_runtime(runtime);
	BORG_Random_seed(seed);
	BORG_Archive result = BORG_Algorithm_ms_run(problem);

	if (result != NULL){
		outputFile = fopen(outputFilename, "w");
		if(!outputFile){
			BORG_Debug("Unable to open final output file\n");
		}
		BORG_Archive_print(result, outputFile);
		BORG_Archive_destroy(result);
		fclose(outputFile);
	}

	  BORG_Algorithm_ms_shutdown();
	  BORG_Problem_destroy(problem);

	  return EXIT_SUCCESS;
 }

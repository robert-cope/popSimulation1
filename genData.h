#include <math.h>
#include <stdio.h>


//
//    Copyright 2013 Robert Cope cope.robert.c@gmail.com
//
//    This file is part of popSimulation1.
//
//    popSimulation1 is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    popSimulation1 is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with popSimulation1.  If not, see <http://www.gnu.org/licenses/>.
//



#define k 25 

typedef struct Individual *individual;
typedef struct Individual{
    int id;
    int sex; //0 -> female -> m(aternal); 1 -> male -> p(aternal); 
    int genotype[k][2];
    int birth;
    int birthPop;
    int lastChild;
    individual m;
    individual p;
    int prev;
    double strength;
    int alive;
} Individual;

typedef struct Freq{
    double aFreqs[k][18];
} Freq;


typedef struct Freq *frequencies;

double genAlleles(int);
double nextDouble();
individual genInd(int*,individual, int);
individual genIndB(int*,individual, individual, int,int,int);
void alleleHist(individual[], int);
void printInd(individual, FILE*,int,int);
//void printInd(individual, FILE*,int, individual*, individual*, individual*, individual*,int,int,int,int);
double deathProb(int, double);
int compare_doubles(const void*,const void*);
double nextGamma(double, double);
double genG(double);
int compare_inds(const void*,const void*);
int compare_indsRev(const void*,const void*);
int compare_inds2(const void*,const void*);
int compare_ints(const void*,const void*);
frequencies getFreqs(individual[],individual[], individual[], individual[],int, int, int, int);

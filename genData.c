#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "genData.h"
#include <string.h>


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



double *alleleFreqs[k];
int nAlleles[k];
int i,j,s,t;
#define MINAGE 8

double genAlleles(int p){
    //int mMax = 60;
    double totH = 0.0;
    for (i=0;i<k;i++){
        int* tempSt = malloc(100*sizeof(int));
        int nl=0;
        int r = 100;
        int st = 10;

        while (r > 0){
        //double* tempUnifs = malloc((st-1)*sizeof(double));
        //int i;
        //for (i = 0;i<(st-1);i++){
        //    tempUnifs[i]=nextDouble();
        //}
        //qsort(tempUnifs, (st-1), sizeof(double),compare_doubles); //biggest to smallest
        //double binProb = tempUnifs[p/(100/st)];
        double gX, gY;
        if (p <= 50){
        gX = nextGamma(1.2,1);
        gY = nextGamma(1.2/p*(100-p),1);
        } else {
        gX = nextGamma(1.2/(100-p)*p,1);
        gY = nextGamma(1.2,1);
        }
        double binProb = gX / (gX+gY);
        //printf("%f\n",binProb);
        int nBin = 0;
        int c = 0;
        while(c < 100){
        if (nextDouble() < binProb){
            nBin++;
        }
        c++;
        }
        if (nBin > r) {nBin = r;}
        tempSt[nl] = nBin;
        r = r - nBin;
        nl++;
        }
        //printf("%d\n",nl);
        nAlleles[i]=nl;
        alleleFreqs[i] = malloc(nl*sizeof(double));
        //printf("%d\n",nl);
        double ss = 0.0;
        for(j=0;j<nl;j++){
            alleleFreqs[i][j]=tempSt[j]/100.0;
            printf("%d %d %f\n",i,j,alleleFreqs[i][j]);
            ss += alleleFreqs[i][j]*alleleFreqs[i][j];

        }
            printf("h: %f\n",1.0-ss);
            printf("--------------\n");
            totH += (1-ss);
    }
    printf("avg Het: %f \n",totH/k);
    return totH/k;
}

void alleleHist(individual checkInds[],int c){
    int * checkFreqs[k];
    double ssT = 0.0;
    for(i=0;i<k;i++){
        checkFreqs[i]=malloc((nAlleles[i])*sizeof(int));
        for(j=0;j<nAlleles[i];j++){checkFreqs[i][j]=0;}
        for(s=0;s<c;s++){
            int a1=checkInds[s]->genotype[i][0];
            int a2=checkInds[s]->genotype[i][1];
            if(a1==nAlleles[i]){printf("a1");}
            if(a2==nAlleles[i]){printf("a2");}
            checkFreqs[i][checkInds[s]->genotype[i][0]]++;
            checkFreqs[i][checkInds[s]->genotype[i][1]]++;
        }
        //for(j=0;j<nAlleles[i];j++){printf("%d %d - %f\n",i,j,checkFreqs[i][j]/(2.0*c));}
        double ss=0.0;
        for(j=0;j<nAlleles[i];j++){ss+=(checkFreqs[i][j]/(2.0*c)-alleleFreqs[i][j])*(checkFreqs[i][j]/(2.0*c)-alleleFreqs[i][j]);}
        //printf("----\n");
        //printf("%f\n",ss);
        ssT+=ss;
        //printf("----\n");
    }
    //printf("ssT: %f\n",ssT);
}

individual genInd(int *Nid,individual eInd,int loc){
    //Individual ind = *malloc(sizeof(Individual));
    individual ind = malloc(sizeof *ind);
    ind->birthPop=loc;
    ind->id=*Nid;
    (*Nid)++;
    ind->sex = nextDouble()<0.5?1:0;
    ind->birth = -12;
    ind->lastChild=-5;
    for(i=0;i<k;i++){
        for(j=0;j<2;j++){
        double res = nextDouble();
        //printf("-%f\n",res);
        double c = 0.0;
        int a = 0;
        ind->genotype[i][j]=-1;
        while(a<nAlleles[i]){
            c+= alleleFreqs[i][a];
            //printf("%f\n",c);
            if (res < c){
                ind->genotype[i][j]=a;
                break;
            }
            a++;
        }
        if(ind->genotype[i][j]==-1){
        ind->genotype[i][j]=a;
        }
        if(ind->genotype[i][j]==nAlleles[i]){printf("wtf");}
        //printf("%d %d - %d\n",i,j,ind->genotype[i][j]);
        }
    }
    ind->m=eInd;
    ind->p=eInd;
    ind->prev=eInd->id;
    ind->p->prev=ind->m->id;
    ind->m->prev=ind->p->id;
    ind->alive=1;
    ind->strength= nextDouble();
    //printf("-------");
    return ind;
}

individual genIndB(int *Nid,individual m, individual p, int t,int s,int loc){
    individual ind = malloc(sizeof *ind);
    ind->birthPop=loc;
    ind->id=*Nid;
    (*Nid)++;
    ind->sex=s;
    ind->birth=t;
    ind->lastChild=0;
    for(i=0;i<k;i++){
        ind->genotype[i][0]=nextDouble()<0.5?m->genotype[i][0]:m->genotype[i][1];
        ind->genotype[i][1]=nextDouble()<0.5?p->genotype[i][0]:p->genotype[i][1];
    }
    ind->m=m;
    ind->p=p;
    ind->m->prev=p->id;
    ind->p->prev=m->id;
    ind->alive=1;
    ind->strength=nextDouble();
    printf("%d (%d %d)\n",ind->id,m->id,p->id);
    //printf("%d (%d %d) %f\n",ind->id,m->id,p->id,relEstLynch(m,p,fList));
    return ind;
}

void printInd(individual ind, FILE *fp, int ti,int pop){
    char* gt;
    char* tgt;
    gt=malloc(k*8*sizeof(char));
    tgt=malloc(8*sizeof(char));
    //printf("%d\n",ind->genotype[0][0]);
    //printf("%d\n",ind->genotype[0][1]);
    sprintf(gt,"%d.%d ",ind->genotype[0][0],ind->genotype[0][1]);
    for(t=1;t<k;t++){
    
    sprintf(tgt,"%d.%d ",ind->genotype[t][0],ind->genotype[t][1]);
    strcat(gt,tgt);
    }
    int sampTime = ind->birth + (int)(ti-ind->birth)*nextDouble() + 1;
    int sizeClass = sampTime - ind->birth < 4? 2:1;
    sizeClass = sampTime - ind->birth > 8? 0 : sizeClass;
    fprintf(fp,"%d (%d %d) %d %d %d %d %d %d %s\n",ind->id, ind->m->id,ind->p->id, ind->birthPop, sampTime, sizeClass, ind->sex, ind->birth,pop,gt);
    //printf("%d %d %d \n",ind->id, ind->sex, ind->birth);
}

frequencies getFreqs(individual M1[], individual M2[], individual F1[], individual F2[],int nM1,int nM2,int nF1,int nF2){
    frequencies fStore; 
    fStore = malloc(sizeof(*fStore));
    int **sTemp = malloc(25*sizeof(int*));
    for (i=0;i<25;i++){
        sTemp[i]=malloc(18*sizeof(int));
        for (j=0;j<18;j++){
            sTemp[i][j]=0;
        }
    }
    for (i=0;i<nM1; i++){
        for (j=0;j<25;j++){
            int g1 = M1[i]->genotype[j][0];
            int g2 = M1[i]->genotype[j][1];
            sTemp[j][g1]++;
            sTemp[j][g2]++;
        }
    }
    for (i=0;i<nM2; i++){
        for (j=0;j<25;j++){
            int g1 = M2[i]->genotype[j][0];
            int g2 = M2[i]->genotype[j][1];
            sTemp[j][g1]++;
            sTemp[j][g2]++;
        }
    }
    for (i=0;i<nF1; i++){
        for (j=0;j<25;j++){
            int g1 = F1[i]->genotype[j][0];
            int g2 = F1[i]->genotype[j][1];
            sTemp[j][g1]++;
            sTemp[j][g2]++;
        }
    }
    for (i=0;i<nF2; i++){
        for (j=0;j<25;j++){
            int g1 = F2[i]->genotype[j][0];
            int g2 = F2[i]->genotype[j][1];
            sTemp[j][g1]++;
            sTemp[j][g2]++;
        }
    }
    double heSum = 0.0;
    for (i = 0; i<25; i++){
        for (j=0; j<18;j++){
            //printf("%d\n", sTemp[i][j]);
            fStore->aFreqs[i][j]=sTemp[i][j]/(2.0*(nM1+nM2+nF1+nF2));
            heSum += fStore->aFreqs[i][j]*fStore->aFreqs[i][j];
            //printf("%f\n",fStore->aFreqs[i][j]);
        }
    }
    printf("he: %f\n",1.0-0.04*heSum);
    return fStore;
}



double nextDouble(){
    return (double)rand()/(double)RAND_MAX; 
}

double deathProb(int age, double base){
    double d=base;
    if (age < MINAGE){
        d = base * MINAGE/2;
    }
    if (age > 4 * MINAGE){
        d = base * age;
    }
    if (age > 6 * MINAGE){
        d = base * age * 2;
    }
    return d;
}

     int
     compare_doubles (const void *a, const void *b)
     {
       const double *da = (const double *) a;
       const double *db = (const double *) b;
     
       return (*da > *db) - (*da < *db);
     }

     int
     compare_ints (const void *a, const void *b)
     {
       const int *da = (const int *) a;
       const int *db = (const int *) b;
     
       return (*da > *db) - (*da < *db);
     }

double genG(double d){
    double v0 = M_E / (M_E+d);
    double v1 = nextDouble(), v2 = nextDouble(), v3 = nextDouble();
    double xi;
    double eta;
    if (v1 < v0){
        xi = pow(v2, (1/d));
        eta = pow(v3*xi,(1-d));
    } else {
        xi = 1 - log(v2);
        eta = v3*pow(M_E,(-xi));
    }
    if (eta > pow(xi,(d - 1))* pow(M_E,(-xi))){
        return genG(d);
    } else {
        return xi;
    }
}
double nextGamma(double kT, double theta){
    int kf = floor(kT);
    double r = kT - kf;
    double s = 0;
    int i;
    for (i=0;i<kf;i++){
        s+=log(nextDouble());
    }
    //printf("%f %f\n",kT,s);
    return theta * (genG(r) - s);
}



int compare_inds (const void *a, const void *b)
{
  double temp = (*(const individual*)a)->strength - (*(const individual*)b)->strength;
  if (temp > 0)
    return -1;
  else if (temp < 0)
    return 1;
  else
    return 0;
}


int compare_indsRev (const void *a, const void *b)
{
  double temp = (*(const individual*)a)->strength - (*(const individual*)b)->strength;
  if (temp < 0)
    return -1;
  else if (temp > 0)
    return 1;
  else
    return 0;
}

int compare_inds2 (const void *a, const void *b)
{
  int temp = (*(const individual*)a)->id - (*(const individual*)b)->id;
  if (temp > 0)
    return -1;
  else if (temp < 0)
    return 1;
  else
    return 0;
}



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "genData.h"
#include <search.h>



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



#define DEATHRATE 0.001
#define MINAGE 8
#define INTERBIRTH 3
int i_m,j_m,s,t;
double AvgHet;
int main(int argc, char *argv[]){



    if (argc < 5) {
        printf("insufficient args\n");
        exit(0);
    }

    int het;
    if(EOF == sscanf(argv[1], "%d", &het))
    {
    //error
    }

    int propSampled;
    if(EOF == sscanf(argv[2], "%d", &propSampled))
    {
    //error
    }
    printf("%d\n",propSampled);

    int totPop;
    int m_max1 = 100, f_max1=100;
    if(EOF == sscanf(argv[3], "%d", &totPop))
    {
    //error
    }
    m_max1 = floor(totPop/2.0);
    f_max1=floor(totPop/2.0);


    char* filedest = argv[4];

    int *nextId;
    individual eInd = malloc(sizeof *eInd);
    eInd->id = -1;
    nextId=malloc(sizeof(int));
    *nextId=0;
    srand(time(NULL));
    //srand(12345);
    AvgHet = genAlleles(het);
    int testN = 100000;
    individual baseInds1[testN];
    int dead[testN];
    size_t deadcounter=0;
    for(i_m=0;i_m<testN;i_m++){
    //printf("%d\n",i_m);
    baseInds1[i_m]=genInd(nextId,eInd,-1);
    }
    alleleHist(baseInds1,testN);

    //setup population sizes
    int m_current1=0,f_current1=0;
    individual M1[2*m_max1];
    individual F1[2*f_max1];
    //fill some initial stuff
    for (i_m=0;i_m<testN;i_m++){
    //printf("%d\n",i_m);
        if(baseInds1[i_m]->sex==0){
            if(f_current1 <f_max1){
                F1[f_current1]=baseInds1[i_m];
                f_current1++;
            }
        } else {
            if(m_current1 <m_max1){
                M1[m_current1]=baseInds1[i_m];
                m_current1++;
            }
        }

    }



    //for 10 timesteps
    int ti = 10;
    for(j_m = 0;j_m<200;j_m++){
    ti++;
    //printf("%d\n",j_m);
    //kill some individuals
        //m1 & f1
        for(i_m=0;i_m<m_current1-1;i_m++){
            if(nextDouble()<deathProb(ti-M1[i_m]->birth,DEATHRATE)){
                M1[i_m]->alive=0;
                dead[deadcounter]=M1[i_m]->id;
                deadcounter++;
                for(s=i_m;s<m_current1;s++){
                    M1[s]=M1[s+1];
                }
                m_current1--;
            }
        }
        for(i_m=0;i_m<f_current1-1;i_m++){
            if(nextDouble()<deathProb(ti-F1[i_m]->birth,DEATHRATE)){
                F1[i_m]->alive=0;
                dead[deadcounter]=F1[i_m]->id;
                deadcounter++;
                for(s=i_m;s<f_current1;s++){
                    F1[s]=F1[s+1];
                }
                f_current1--;
            }
        }


        //fill the populations to their max with births
        //m1
        while(m_current1<m_max1){
            individual pot_p = M1[(int)(nextDouble()*m_current1)];
            individual pot_m = F1[(int)(nextDouble()*f_current1)];
            if (ti-pot_p->birth<MINAGE||ti-pot_m->birth<MINAGE||ti-pot_m->lastChild<INTERBIRTH){
                continue;
            }
            M1[m_current1]=genIndB(nextId,pot_m,pot_p,ti,1,1);
            m_current1++;
        }
        //f1
        while(f_current1<f_max1){
            individual pot_p = M1[(int)(nextDouble()*m_current1)];
            individual pot_m = F1[(int)(nextDouble()*f_current1)];
            if (ti-pot_p->birth<MINAGE||ti-pot_m->birth<MINAGE||ti-pot_m->lastChild<INTERBIRTH){
                continue;
            }
            F1[f_current1]=genIndB(nextId,pot_m,pot_p,ti,0,1);
            f_current1++;
        }
}
     
    FILE *fpM1;
    printf("%f\n",AvgHet);
    printf("%d\n",ti);
    printf("%s\n",filedest);
    if(fpM1 = fopen(filedest,"w")){
    fprintf(fpM1, "AvgHet: %f\n",AvgHet);
    for(i_m=0;i_m<m_current1;i_m++){
    if (nextDouble() < propSampled/100.0){
    individual tempIndX = M1[i_m];
    printInd(tempIndX,fpM1,ti,1);
    }
    }
    
    for(i_m=0;i_m<f_current1;i_m++){
    if (nextDouble() < propSampled/100.0){
    //printInd(F1[i_m],fpM1,ti);
    printInd(F1[i_m],fpM1,ti,1);
    }
    }
    //}

    }
    printf("end");
    fclose(fpM1);

}




// Author: Wes Kendall
// Copyright 2011 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// MPI_Send, MPI_Recv example. Communicates the number -1 from process 0
// to process 1.
//
#include <mpi.h>
#include <iostream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <fstream>
#include <cmath>
#include <sys/resource.h>
#include <vector>
using namespace std;

/* come ising 2d ma adesso tipo blume capel,punto di reticolo vale 0 se non c'è spin
si calcolano una volta all'inizio possibili valori di coseno di differenza angoli, così dopo basta prendere elemento che
si vuole di questo vettore. bondenergy dà solo energia tra i e j dovuta a coseno */


// PARAMETRI DA VARIARE
const int p=10; // numero angoli, pare che da 10 in su T_c sia praticamente identica a quella dell'XY
const int Liniz=16; //un po' ridondanti, valore di Lx e Ly nelle righe sopra è inutile ma devono essere definiti. si usano di fatto Liniz e Lnext
const int Lnext=24;
const int Lx=Liniz;
const int Lstep=Lnext-Liniz;
const int nmontecarlo=2000;
const int steptermalizz=1000;
const double Delta=.1; //rapporto tra coupling trasverso e quello sui piani
//const double numgiri=1;
//const double D= 1.02;// piùmeno 60 è abbastanza per rendere tutto vuoto, cioè 0, (se D negativo)
const double betacrit=0.45416474; // quella dell'XY, senza buchi
double beta=0.1; //quello variabile
const double betamin=0.3;
const double betamax=2.2;
const int numeropassibeta=10; // meglio usare un int per dire quanti valori si prendono tra betamin e betamax che un double che indica larghezza dello step
double valoribeta[numeropassibeta];
double matricebinder[numeropassibeta][2];

// PARAMETRI CHE SI CALCOLA LUI
//const int Lx=Ly;        // SE SI CAMBIA QUESTO BISOGNA CAMBIARE NEXT-PREVIOUS E PURE GENERATORE RANDOM DI SITO DA 1 A LY
int N=Lx*Lx*Lx;   //int sitivariabili=Lx*Ly*(Lx-2);
const double pi= 4.*atan(1);
const int stepdopotermalizz=nmontecarlo-steptermalizz;
int fr=(Lx-1)/16;   // serve per reticolo di corr2

// PESI CHE VENGONO SETTATI DA SETTAPESI
long int rng;
double v[p];
double vseni[p];
double pesi[p];
double pesiz[p];
double pesim[p];// pesi di boltzmann per vari possibili valori di angolo e loro inversi, per non ricalcolarli
double pesimz[p];
double pbond[2*p][2*p];

static bool cluster[Lx][Lx][Lx];
auto  m=new int[Lx][Lx][Lx];

int megaspin=1; // da 1 a p;


#include "funzioniClockStrati.h"

void randfill(int Lx)
{
  for(int k=0;k<Lx;k++)
    {
      for(int j=0;j<Lx;j++)
     {
       for(int i=0;i<Lx;i++)
       {
          m[i][j][k]=randomspin();
       }
      }
    }
}



int  main ()
{
    MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if (world_size < 2)
  {
    cout<<"non ci sono abbastanza core";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

    int tiniz=time(NULL);
    //numero
    double betastep= (betamax-betamin)/(numeropassibeta-1);
    for (int i = 0;  i < numeropassibeta; i++)
    {
      valoribeta[i] = betamin + i * betastep ;
    }
      ofstream nomequalsiasi; nomequalsiasi.open("XYbinder.txt");
      if(world_rank==world_size-1)
      {
        cout<< "dimensione del sistema e'" <<Lx<<"x"<<Lx<<"x"<<Lx<<endl
        <<"numero di sweep e' "<<nmontecarlo<<endl
        <<"numero di core e'"<<world_size<<endl;


        nomequalsiasi << "beta,Delta,mediafinale,media2finale,media4finale,binder" << endl;


      }
      int indiceL=0;

    for(int Lvar=Liniz; Lvar<Lvar /*qui dovrebbe esserci <Lnext, CAMBIARE QUANDO SI SA COME CANCELLARE E RIDEFINIRE m */+1; Lvar+=2 /*Lstep dovrebbe essere lui, non 2*/)
    {
       N=pow(Lvar,3);

      for(int indicebeta=0; indicebeta < numeropassibeta ; indicebeta++)
      {
        double mediafinale=0;   // magnetizzazione complessiva
        double media2finale=0;   // magnetizzazione complessiva
        double media4finale=0;   // magnetizzazione complessiva
        //beta=betamin;beta<=betamax;beta+=betastep
        settapesi(valoribeta[indicebeta]);


        if (world_rank < world_size-1)
        {
              rng=-time(NULL)*world_rank;     // settare entrambi in modo che dipenda da che core si sta usando
              srand (time(NULL)*world_rank);



              //cout<< "sta succedendo qualcosa in"<<world_rank<<endl;

             double medietemp[3]={0}; //media, media quadra, media^4

                //filldafile();
                randfill(Lvar);
                double unamedia;
                //ofstream nomequalsiasi; nomequalsiasi.open("dati.txt");
                for(int mosse=0;mosse<nmontecarlo; mosse++) //for(int mosse=0;mosse<nmontecarlo;mosse++)
            {
              SweepEWolff();
               if(mosse>steptermalizz)
                {
                  unamedia=mean();
                  medietemp[0]+= unamedia;
                  medietemp[1]+= pow(unamedia,2);
                  medietemp[2]+=pow(unamedia,4);
                }

            }   // finito un giro di montecarlo

                medietemp[0]/=(stepdopotermalizz-1);
                medietemp[1]/=(stepdopotermalizz-1);
                medietemp[2]/=(stepdopotermalizz-1);

                    // manda vettore magnetizzazione, vettore energie, matriciona correlazioni, usando tag diversi
              MPI_Send(
                    /* data         = */ medietemp,//&mediatemp,
                    /* count        = */ 3,
                    /* datatype     = */ MPI_DOUBLE,
                    /* destination  = */ world_size-1,
                    /* tag          = */ world_rank,
                    /* communicator = */ MPI_COMM_WORLD
                      );

        }
          else if (world_rank == world_size-1)
          {

            double ricevi[3];
           //double cosotemp=0;
          for(int i=0;i<world_size-1;i++)
          {
            MPI_Recv(                     // ultimo core riceve magnetizzazioni
            /* data         = */ ricevi,//&cosotemp,
            /* count        = */ 3,
            /* datatype     = */ MPI_DOUBLE,
            /* source       = */ i,
            /* tag          = */ i,
            /* communicator = */ MPI_COMM_WORLD,
            /* status       = */ MPI_STATUS_IGNORE);

            mediafinale+=ricevi[0];
            media2finale+=ricevi[1];
            media4finale+=ricevi[2];

          } //fine ciclo su vari core che mandano roba all'ultimo
          mediafinale/=(world_size-1);
          media2finale/=(world_size-1);
          media4finale/=(world_size-1);

          //ofstream nomequalsiasi; nomequalsiasi.open("XYbinder.txt");
          nomequalsiasi << Lx << ","<< beta << ","<< mediafinale <<","<< media2finale<< ","<<media4finale<<","<<media4finale/ (pow(media2finale,2)) << endl;
          matricebinder[indicebeta][indiceL] = media4finale/ (pow(media2finale,2));
          indiceL++;
          /*if(Lvar==Liniz) SCOMMENTARE QUANDO SI SA RIDEFINIRE m
          {
            auto  mnuovo= new int[Lnext][Lnext][Lnext];
            free(m);
            m=mnuovo;
          }*/
        }
      }//fine loop sui beta
    } //fine loop sugli L
      if(world_rank== world_size-1)
	{
	  nomequalsiasi.close();
	  cout<<endl<< "tempo trascorso e'"<<time(NULL)-tiniz<<endl;
	}
  MPI_Finalize();

      //mediafinale+=mediatemp;



    /*for(int u=1;u<Lx/2+1;u++)
    {
      cout<<magnxfinale[u]<<" , ";
       //<<sqrt(meansq[u]/numgiri-pow(magnxfinale[u]/numgiri,2))<<endl;
    }*/
    /*cout<<endl<<"le deviazioni standard sono"<<endl;
     for(int u=1;u<Lx/2+1;u++)
    {
      cout<<sqrt(meansq[u]/numgiri-pow(magnxfinale[u]/numgiri,2))<<" , ";
    }*/


  /*ofstream nomequalsiasi;
  nomequalsiasi.open("dati24con1E4step2giri.txt");
  for(int u=1;u<Lx/2+1;u++)
    {
      nomequalsiasi<<magnxfinale[u]/numgiri<<" , "
       <<sqrt(meansq[u]/numgiri-pow(magnxfinale[u]/numgiri,2))<<endl;
    }


  nomequalsiasi.close();
   */
  return 0;
}

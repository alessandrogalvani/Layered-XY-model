
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double randomreal(long *idum) // LA SUA VARIABILE DEV'ESSERE & coso, dove coso definito all'inizio coso=-time(NULL) SERVE IL MENO
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else {return temp;}
}


int randint(int max){return 1+(int) (double(max)*rand()/(RAND_MAX+1.0)); }
int randombit(){return randint(2);}
int randomspin(){return randint(p);}
int randomx(){return randint(Lx-2);}
int randomy(){return randint(Lx)-1;}




//std::uniform_int_distribution<int> cento(1,100);

int next(int i)         //se si vuole avere pbc anche su x serve Lx=Lx oppure definire nuovi next/previous
{
  if (i==Lx-1) return 0;
  else return i+1;
}
int previous(int i)
{
  if (i==0) return Lx-1;
  else return i-1;
}
int ndaijk(int i, int j, int k){return k+j*Lx+i*Lx*Lx;}
int idan(int n){return n/(Lx*Lx);}
int jdan(int n){return n%(Lx*Lx)/Lx;}
int kdan(int n){return n%Lx;}

void sitodopo(int &i,int &j, int &k)

{
  if (i<Lx-1)

    i++;

  else
    {
      i=0;
      if(j<Lx-1)
	j++;
      else
	{   j=0;
	  if(k<Lx-1) k++;
	  else k=0;
	}
    }

}// si usa per far crescere indici ijk in spinflip1 e spinflip2
/*int vicini[N][6];
void settavicini()
{
  for(int n=0;n<N;n++)
    {
      int i=idan(n);
      int j=jdan(n);
      int k=kdan(n);
      vicini[n][0]= ndaijk((i+1)%Lx, j , k );
      vicini[n][1]= ndaijk( (i+Lx-1)% Lx ,j , k);
      vicini[n][2]= ndaijk(i, (j+1)%Ly , k );
      vicini[n][3]= ndaijk(i, (j+Ly-1)% Ly, k);
      vicini[n][4]= ndaijk(i, j, (k+1)%Ly );
      vicini[n][5]= ndaijk(i, j, (k+Ly-1)%Ly );

    }
}
*/

double mean()
{
  double mediax=0; double mediay=0;

  for(int i=0;i<Lx;i++)
    {
      for(int j=0;j<Lx;j++)
      {
        for(int k=0;k<Lx;k++)
          {
             mediax+=v[abs(m[i][j][k]-1)]; mediay+=vseni[abs(m[i][j][k]-1)];
          }
      }
    }
  mediax/=N; mediay/=N;

  return mediax*v[megaspin-1]+mediay*vseni[megaspin-1];
}
double energialocale(int i, int j, int k)
{
  double E=0;
  int sito=m[i][j][k];
  if(sito!=0)
  {
    if(i>0)if(m[i-1][j][k]!=0)    E+=v[ abs(sito-m[i-1][j][k]) ] ;
    if(i<Lx-1)if(m[i+1][j][k]!=0) E+=v[ abs(sito-m[i+1][j][k]) ] ;
    if(m[i][next(j)][k]!=0)       E+=v[ abs(sito-m[i][next(j)][k]) ] ;
    if(m[i][previous(j)][k]!=0)   E+=v[ abs(sito-m[i][previous(j)][k]) ] ;
    if(m[i][j][next(k)]!=0)       E+=v[ abs(sito-m[i][j][next(k)]) ] ;
    if(m[i][j][previous(k)]!=0)   E+=v[ abs(sito-m[i][j][previous(k)]) ] ;

  }
  return E;
}
double magnx(int i)
{
  double mediax=0; double mediay=0;


  for(int j=0;j<Lx;j++)
    {
      for(int k=0;k<Lx;k++)
  {
    if(m[i][j][k]!=0) mediax+=v[m[i][j][k]-1];
    if(m[i][j][k]!=0) mediay+=vseni[m[i][j][k]-1];
  }
    }
  mediax/=(Lx*Lx); mediay/=(Lx*Lx);
  return mediax*v[megaspin-1]+mediay*vseni[megaspin-1];

}


double energiax(int i)
{
  double media=0;

  for(int j=0;j<Lx;j++)
      for(int k=0;k<Lx;k++)
      media+=energialocale(i,j,k);

  return media/(Lx*Lx);

}

void spinflip(int i, int j, int k)      // ergodico
{
  int angolonuovo;
  double expmenodeltaE=1;

  angolonuovo=randomspin();


    expmenodeltaE*= pesi[ abs(m[next(i)][j][k]-angolonuovo)];
    expmenodeltaE*= pesi[ abs(m[previous(i)][j][k]-angolonuovo)];
    expmenodeltaE*= pesi[ abs(m[i][next(j)][k]-angolonuovo)];
    expmenodeltaE*= pesi[ abs(m[i][previous(j)][k]-angolonuovo)];
    expmenodeltaE*= pesiz[ abs(m[i][j][next(k)]-angolonuovo)];
    expmenodeltaE*= pesiz[ abs(m[i][j][previous(k)]-angolonuovo)];



      expmenodeltaE*= pesim[ abs(m[next(i)][j][k]-m[i][j][k]) ];
      expmenodeltaE*= pesim[ abs(m[previous(i)][j][k]-m[i][j][k]) ];
      expmenodeltaE*= pesim[ abs(m[i][next(j)][k]-m[i][j][k]) ];
      expmenodeltaE*= pesim[ abs(m[i][previous(j)][k]-m[i][j][k]) ];
      expmenodeltaE*= pesimz[ abs(m[i][j][next(k)]-m[i][j][k]) ];
      expmenodeltaE*= pesimz[ abs(m[i][j][previous(k)]-m[i][j][k]) ];







  if(randomreal(& rng)<expmenodeltaE) {m[i][j][k]=angolonuovo;}
  // printpiano(m,k); }//cout << "spin flippato in "<< angolonuovo << endl;}


  //    cout << "si è flippato" <<peso<< endl;
}

int flipangolo(int a,int r)    // prende numero da 1 a p e dà numero da 1 a p
{
 return ( (r-a+1+p/2+p)%p ) +1;
}

void settapesi(double be)
{
  // vettore con possibili valori coseni
  for(int jj=0;jj<p;jj++)
    {
      v[jj]=cos(2.*pi*double(jj)/p); if( abs(v[jj])<1E-7) v[jj]=0;
      vseni[jj]=sin(2.*pi*double(jj)/p);
    }

  // vettore con possibili pesi boltzmann

  for(int jj=0;jj<p;jj++)
    {
      pesi[jj]=exp(be*v[jj]);
      pesim[jj]=1/pesi[jj];
    }

    //stessa cosa ma su z

      for(int jj=0;jj<p;jj++)
        {
          pesiz[jj]=exp(Delta*be*v[jj]);
          pesimz[jj]=1/pesiz[jj];
        }

  for(int jj=0;jj<2*p;jj++)
    {
      for(int kk=0;kk<2*p;kk++)
        {
    pbond[jj][kk]=1-exp(-2* be *cos(pi*double(jj)/p)* cos(pi*double(kk)/p) );
        if(pbond[jj][kk]<0) pbond[jj][kk]=0;
        }
    }


}




void wolff() // matrice Lx*Lx*Lx cluster serve per controllare che sito non sia gia' stato preso. vettore cluster serve per flippare tutti assieme alla fine. vettore sitidafare e' quello che si svuota gradualmente
{
  for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Lx; j++)
      for (int k = 0; k < Lx; k++)
      cluster[i][j][k] = 0;


  int i0 = randomx();
  int j0 = randomy();
  int k0 = randomy();

  int r= randomspin()-1;
 // per evitare di fare roba ricorsiva si fa un while su siti da fare. ogni terna che indica un sito viene convertita in un numero solo, che è scritto in una base in cui una cifra è in base Lx e le altre due in base Lx. poi le tre funzioni ida,jdan,kdan riconvertono questo numero in i,j,k


  vector<int> sitidafare;
  vector<int> siticluster;

     cluster[i0][j0][k0]=1;
     sitidafare.push_back(ndaijk(i0, j0,k0));
     siticluster.push_back(ndaijk(i0, j0,k0));


  while(!sitidafare.empty())
  {   //prima si determina chi va in cluster, poi si flippa tutti
    int i=idan(sitidafare[0]);
    int j=jdan(sitidafare[0]);
    int k=kdan(sitidafare[0]);





        if(!cluster[previous(i)][j][k])
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[previous(i)][j][k]-1)-r) ] )
          {
          cluster[previous(i)][j][k]=1;
          sitidafare.push_back(ndaijk(previous(i), j,k));
          siticluster.push_back(ndaijk(previous(i), j,k));
          }
        if(!cluster[next(i)][j][k]&&m[next(i)][j][k]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[next(i)][j][k]-1)-r) ] )
          {
          cluster[next(i)][j][k]=1;
          //if(i!=Lx-1||cluster[0][j][k]==0)  // questo e' per evitare che si raggiunga bordo contemporaneamente da due
          sitidafare.push_back(ndaijk(next(i), j,k));
          siticluster.push_back(ndaijk(next(i), j,k));
          }
        if(!cluster[i][previous(j)][k]&&m[i][previous(j)][k]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i][previous(j)][k]-1)-r) ] )
          {
          cluster[i][previous(j)][k]=1;
          sitidafare.push_back(ndaijk(i, previous(j),k));
          siticluster.push_back(ndaijk(i, previous(j),k));
          }

        if(!cluster[i][next(j)][k]&&m[i][next(j)][k]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i][next(j)][k]-1)-r) ] )
          {
          cluster[i][next(j)][k]=1;
          sitidafare.push_back(ndaijk(i, next(j),k));
          siticluster.push_back(ndaijk(i, next(j),k));
          }

        if(!cluster[i][j][next(k)]&&m[i][j][next(k)]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i][j][next(k)]-1)-r) ] )
          {
          cluster[i][j][next(k)]=1;
          sitidafare.push_back(ndaijk(i, j , next(k) ) );
          siticluster.push_back(ndaijk(i, j , next(k) ) );
          }

        if(!cluster[i][j][previous(k)]&&m[i][j][previous(k)]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i][j][previous(k)]-1)-r) ] )
          {
          cluster[i][j][previous(k)]=1;
          sitidafare.push_back(ndaijk(i, j,previous(k) ));
          siticluster.push_back(ndaijk(i, j,previous(k) ));
          }


    sitidafare.erase(sitidafare.begin());
  }// fine del while
//cout<<siticluster.size()<<endl;
  for(vector<int>::iterator q=siticluster.begin();q!=siticluster.end();q++)
  {
    int is=idan(*q);
    int js=jdan(*q);
    int ks=kdan(*q);
   //cout<<is<<" "<<js<<" "<<ks<<endl;
    m[is][js][ks]=flipangolo(m[is][js][ks],r);
  }//cout<<"taglia "<<siticluster.size()<<endl;

}



void SweepEWolff()
{
  int i=0; int j=0; int k=0;

    for(int rr=0;rr<6;rr++)
      {
        for(int sweep=0;sweep<N;sweep++)
          {
            spinflip(i,j,k);
            sitodopo(i,j,k);
          }
         // wolff(); tolto perché bisogna controllare come farlo con coupling su z diverso
      }
}

//int rifletti(int x){return 16-x;}

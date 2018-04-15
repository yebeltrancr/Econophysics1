#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>

//Constantes
const int N = 1000; //Numero de agentes
const int M = 1000; //Dinero total
const int P = pow(10, 4); //Numero de transacciones, pasos de tiempo
const int S = pow(10, 4); //Numero de simulaciones
const double Mprom = double(M)/N; //Dinero promedio por agente
const int m = 20; //Estados calculo entropia
const double lambda = 0.3; //Coeficiente Saving propency

//Vectores
std::vector<double> V (N, Mprom); //dinero*agente inicializado con Mprom
std::vector<double> H (P, 0); //vector entropia
std::vector<double> Vprom (N, 0); //dinero*agente varias simulaciones

//Archivos de Datos, ¡Cambiar nombres para cada simulacion!
//Los archivos se reescriben si no se cambia
std::string Distribution = "datosCyC.dat";
std::ofstream DatosD(Distribution.c_str());

std::string Entropy = "HCyC.dat";
std::ofstream DatosH(Entropy.c_str());

//Funciones que se usan en main
void Intercambio (std::vector<double> & V,double lambda, int i, int j, double r);
double Entropia (std::vector<double> V,int N, int m, double Mx);
double Gini (std::vector<double> V);

int main(){
  //Generadores de numeros aleatorios
  //3 diferentes semillas una para cada punto
  std::random_device r0, r1, r2;
  std::mt19937_64 gen1, gen2, gen3;
  gen1.seed(r0());
  gen2.seed(r1());
  gen3.seed(r2());

  //Distribucione numero entero para elegir agentes i y j
  std::uniform_int_distribution<> disi(0, N-1);
  //Distribución numero real para delta en el intercambio
  std::uniform_real_distribution<> disr(0, 1);

  //Varias simulaciones
  double Mx = 0.0;
  for(int jj = 0; jj < S; ++jj ){
    //Simulacion
    for(int ii = 0; ii < P; ++ii){
      Intercambio(V, lambda, disi(gen1), disi(gen2), disr(gen3));
      if(jj != 0 ) H[ii] += Entropia(V, N, m, Mx); //Entropia*P
    }
    std::sort(V.begin(), V.end()); //Organizacion descendente
    if(jj == 0) Mx = V[N-1] + 0.5;
    for(int kk = 0; kk < N; ++kk){
      Vprom[kk] += V[kk]; //sumar cada componente de los agentes con el mismo dinero
      V[kk] = Mprom; //volver a poner V con la condicion inicial
    }
  }

  //Dividir la suma entre el numero de simulaciones hechas para encontrar
  //una distribucion promedio
  for(int ii = 0; ii < N; ++ii) Vprom[ii] = Vprom[ii]/S;
  for(int ii = 0; ii < P; ++ii) H[ii] = H[ii]/(S-1.0);
  Mx = Vprom[N-1];

  //Generar el archivo con los datos distribución
  for(int ii = 0; ii < N; ++ii) DatosD << Vprom[ii] << "\n";
  for(int ii = 0; ii < P; ++ii) DatosH << ii << "\t" << H[ii] << "\n";

  //Ver variables en la consola Dinero total, entropia y gini
  double Dt = 0;
  for (int i = 0; i < N; ++i) Dt += Vprom[i];
  std::cout << "DineroT: " << Dt << "\n";
  std::cout << "Entropia: " << Entropia(Vprom, N, m, Mx) << "\n";
  std::cout << "Gini: " << Gini(Vprom) << "\n";
  //std::cout << Mx << "\n";

  return 0;
}

void Intercambio (std::vector<double> & V,double lambda , int i, int j, double r){
  //Regla de intercambio modelo CyC
  double delta = (1 - lambda)*(V[i] - r*(V[i]+V[j]));
  if(i!=j){
    V[i] -= delta;
    V[j] += delta;
  }
  if(V[i] <= 0){
    V[i] += delta;
    V[j] -= delta;
  }
}

double Gini (std::vector<double> Vp){

  std::vector<double> V = Vp;

  //Organizar vector
  std::sort(V.begin(), V.end());

  //Población acumulada Xi+1 - Xi
  int Pt = V.size(); //tamaño de la población

  std::vector<double> Xt ;
  for (int i = 0; i < Pt; ++i)
    Xt.push_back(1.0/Pt);

  //Dinero acumulado

  double Dt = 0.0; //Dinero total
  for (int i = 0; i < Pt; ++i) Dt += V[i]; //suma componentes

  std::vector<double> Y = V; //normalización
  for (int i = 0; i < Pt; ++i) Y[i] /= Dt; 

  std::vector<double> AY = Y; //acumulado
  for (int i = 0; i < Pt; ++i){
    if (i == 0) AY[i] = Y[i];
    else AY[i]= AY[i-1] + Y[i];
  }

  std::vector<double> Yt = AY; // Yi+1 + Yi
  for (int i = 0; i < Pt; ++i){
    if (i==1) Yt[i] = AY[i];
    else Yt[i] = AY[i-1] + AY[i];
  }

  double L = 0.0;
  for (int i = 0; i < Pt; ++i) L += Yt[i]*Xt[i];

  return 1 - L;
}

double Entropia (std::vector<double> V, int N, int m, double Mx){
  //Vector frecuencias
  std::vector<double> F (m, 0);

  //delta
  double Delta = Mx/m;

  //Contando agentes en el intervalo
  for(int i = 0; i < N; ++i){
    for(int j = 0; j < m; ++j){
      if (V[i] > Delta*j && V[i] <= (j+1)*Delta ){
        F[j] = F[j] + 1;
        break;
      }
    }
  } 

  //Probabilidad
  for(int i = 0; i < m; ++i) F[i] = F[i]/N;

  //Entropia
  double H = 0.0;
  for (int i = 0; i < m ; ++i){
    if (F[i] != 0){
      H += F[i]*log(F[i]);
    }
  }


  return -H;
}

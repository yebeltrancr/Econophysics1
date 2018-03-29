#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>

//Constantes
const int N = 1000; //Numero de agentes
const int M = 1000; //Dinero total
const int P = 10E5; //Numero de transacciones, pasos de tiempo
const double Mprom = double(M)/N; //Dinero promedio por agente

//Archivos de Datos
std::string Data = "datos.dat";
std::string Entropy = "entropia.dat";

std::ofstream Datos(Data.c_str());
std::ofstream Entropia(Entropy.c_str());

//Funciones que se usan en main
void Intercambio (std::vector<double> & V, int i, int j, double r);

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

  //Vector de dinero por agente inicializado con el dinero promedio
  std::vector<double> V (N, Mprom);

  //Simulacion
  for(int ii = 0; ii < P; ++ii){
    Intercambio(V, disi(gen1), disi(gen2), disr(gen3));
  }

  for(int i = 0; i < N; ++i) Datos << V[i] << "\n";

  return 0;
}

void Intercambio (std::vector<double> & V, int i, int j, double r){
  double delta = r*Mprom;
  if(i!=j){
    V[i] -= delta;
    V[j] += delta;
  }
  if(V[i] <= 0){
    V[i] += delta;
    V[j] -= delta;
  }
}

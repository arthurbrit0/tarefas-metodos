#include <iostream>
#include <iomanip>
#include <cmath>
#include "funcoes.h"

int main() {
    std::cout << std::fixed << std::setprecision(4);

    // ------------------------------------------------------------------
    // q1: f(a) = -(e^a)/2 + 2 cos(a)
    //    (a) isolamento
    //    (b) bisseção
    //    (c) posição falsa
    //    (d) verificar resultados
    // ------------------------------------------------------------------

    double eps = 1.0e-4;
    int criterioParada = 50;

    double aEsq = 0.0;
    double aDir = 1.0;

    std::cout << "[Q1] f(" << aEsq << ") = " << fPendulo(aEsq) << "\n";
    std::cout << "[Q1] f(" << aDir << ") = " << fPendulo(aDir) << "\n";

    // bisseção
    double raizB = bissecao(aEsq, aDir, eps, criterioParada, fPendulo);
    std::cout << "[Q1-b] Bissecao: raiz = " << raizB 
              << "; f(raiz) = " << fPendulo(raizB) << "\n";

    // posição falsa
    double raizF = posicaoFalsa(aEsq, aDir, eps, criterioParada, fPendulo);
    std::cout << "[Q1-c] Posicao Falsa: raiz = " << raizF
              << "; f(raiz) = " << fPendulo(raizF) << "\n";

    // ------------------------------------------------------------------
    // q2
    //   (a) newton raphson
    //   (b) secante (a0=0.5, a1=1.0)
    //   (c) checar se ultrapassa pi/4
    // ------------------------------------------------------------------

    criterioParada = 50; 
    double a0 = 0.5;
    double raizNR = newtonRaphson(a0, eps, criterioParada, fPendulo, dfPendulo);
    std::cout << "\n[Q2-a] Newton raphson: raiz = " << raizNR
              << "; f(raiz) = " << fPendulo(raizNR) << "\n";

    double raizSec = secante(0.5, 1.0, eps, criterioParada, fPendulo);
    std::cout << "[Q2-b] secante: raiz = " << raizSec
              << "; f(raiz) = " << fPendulo(raizSec) << "\n";

    double limite = M_PI / 4.0;
    if (std::fabs(raizNR) > limite) {
        std::cout << "[Q2-d] Newton raphson rompeu limite de pi/4\n";
    }
    if (std::fabs(raizSec) > limite) {
        std::cout << "[Q2-d] secante rompeu limite de pi/4\n";
    }

    // ------------------------------------------------------------------
    // q3: f(a) = a^3 - 9a + 3
    //   (a) polinomio
    //   (b) ponto fixo
    // ------------------------------------------------------------------

    std::cout << "\n[Q3] polinomio: f(a)= a^3 - 9a + 3\n";

    double xE = 0.0;
    double xD = 1.0;

    std::cout << "f(0) = " << fPolinomio(0) << "; f(1) = " << fPolinomio(1) << "\n";
    double raizPoli = bissecao(0.0, 1.0, eps, 50, fPolinomio);
    std::cout << "[Q3-a] bissecao no polinômio: raiz = " << raizPoli
              << "; f(raiz) = " << fPolinomio(raizPoli) << "\n";

    auto gPoli = [](double a){
        double alpha = 0.1;
        return a - alpha*fPolinomio(a);
    };

    static double alpha = 0.1;
    static auto gPoliStatic = [](double a){
        return a - alpha*fPolinomio(a);
    };

    double raizPF = pontoFixo(1.0, eps, 50, gPoliStatic);
    std::cout << "[Q3-b] ponto fixo no polinômio: raiz = " << raizPF
              << "; f(raiz) = " << fPolinomio(raizPF) << "\n";

    // ------------------------------------------------------------------
    // q4: f(x)= x - xln(x) = 0
    //   (a) newton raphson
    //   (b) ponto fixo
    // ------------------------------------------------------------------

    std::cout << std::setprecision(6);
    eps = 1.0e-5;
    criterioParada = 100;
    double x0 = 3.0; 

    double raizN4 = newtonRaphson(x0, eps, criterioParada, fLog, dfLog);
    std::cout << "\n[Q4-a] newton raphson: raiz = " << raizN4
              << "; f(raiz) = " << fLog(raizN4) << "\n";

    double raizP4 = pontoFixo(x0, eps, criterioParada, gLog);
    std::cout << "[Q4-b] ponto fixo:     raiz = " << raizP4
              << "; f(raiz) = " << fLog(raizP4) << "\n";
    return 0; // explode
}



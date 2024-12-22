#include "funcoes.h"
#include <iostream>
#include <cmath>
#include <limits>

//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//      funcoes pra q1 e q2
//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

double fPendulo(double a) {
    // f(a) = -(e^a / 2) + 2*cos(a)
    return -(std::exp(a)/2.0) + 2.0*std::cos(a);
}
double dfPendulo(double a) {
    // derivada: d/d(a)[-(exp(a)/2) + 2 cos(a)] = -(exp(a)/2) - 2 sin(a)
    return -(std::exp(a)/2.0) - 2.0*std::sin(a);
}

//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//      2) funcao pra q3
//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

double fPolinomio(double a) {
    // f(a) = a^3 - 9a + 3
    return a*a*a - 9.0*a + 3.0;
}
double dfPolinomio(double a) {
    // d/d(a)[a^3 - 9a + 3] = 3a^2 - 9
    return 3.0*a*a - 9.0;
}

//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//      3) funcoes pra q4
//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

double fLog(double x) {
    // f(x) = x - x ln(x) = x(1-ln(x))
    return x - x*std::log(x);
}
double dfLog(double x) {
    // d/dx[x - xln(x)] = 1 - [ln(x) + 1] = -ln(x)
    if (x <= 0.0) {
        std::cerr << "erro: x<=0 (ln indefinido)\n";
        return std::numeric_limits<double>::quiet_NaN();
    }
    return -std::log(x);
}
double gLog(double x) {
double alfa = 0.1;
    return x - alfa*fLog(x);
}

//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// isolamento
//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

bool isolamentoAnalitico(double &aEsq, double &aDir, double passo, double (*func)(double)) {
    double f1 = func(aEsq);
    for (double x = aEsq+passo; x <= aDir; x += passo) {
        double f2 = func(x);
        if (f1*f2 < 0.0) {
            aEsq = x - passo;
            aDir = x;
            return true;
        }
        f1 = f2;
    }
    return false;
}

//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//      metodos numericos
//=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

// bissecao
double bissecao(double aEsq, double aDir, double eps, int criterioParada, double (*func)(double)) {
    if (func(aEsq)*func(aDir) > 0.0) {
        std::cerr << "[Bissecao] erro: f(aEsq)*f(aDir) > 0\n";
        return std::numeric_limits<double>::quiet_NaN();
    }

    double xCentral = 0.0;
    for (int i=0; i < criterioParada; i++) {
        xCentral = 0.5*(aEsq + aDir);
        double fCentral = func(xCentral);

        if (std::fabs(fCentral) < eps)
            break;

        if (func(aEsq)*fCentral < 0.0)
            aDir = xCentral;
        else
            aEsq = xCentral;
    }
    return xCentral;
}

// posicao falsa
double posicaoFalsa(double aEsq, double aDir, double eps, int criterioParada, double (*func)(double)) {
    if (func(aEsq)*func(aDir) > 0.0) {
        std::cerr << "[Posicao falsa] Erro: f(aEsq)*f(aDir) > 0\n";
        return std::numeric_limits<double>::quiet_NaN();
    }

    double xRaiz = 0.0;
    for (int i = 0; i < criterioParada; i++) {
        double fEsq = func(aEsq);
        double fDir = func(aDir);
        xRaiz = aDir - fDir*(aEsq - aDir)/(fEsq - fDir);

        double fRaiz = func(xRaiz);
        if (std::fabs(fRaiz) < eps)
            break;

        if (fEsq * fRaiz < 0.0)
            aDir = xRaiz;
        else
            aEsq = xRaiz;
    }
    return xRaiz;
}

// newton raphson
double newtonRaphson(double x0, double eps, int criterioParada, double (*func)(double), double (*dfunc)(double)){
    double x = x0;
    for (int i = 0; i < criterioParada; i++) {
        double fx  = func(x);
        double dfx = dfunc(x);

        if (std::fabs(dfx) < 1.0e-14) {
            std::cerr << "[Newton] Erro: derivada ~0.\n";
            return x;
        }
        double x1 = x - fx/dfx;

        if (std::fabs(func(x1)) < eps)
            return x1;

        x = x1;
    }
    return x;
}

// secante
double secante(double x0, double x1, double eps, int criterioParada, double (*func)(double)){
    double f0 = func(x0);
    double f1 = func(x1);

    for (int i = 0; i < criterioParada; i++) {
        if (std::fabs(f1 - f0) < 1.0e-14) {
            std::cerr << "[Secante] Erro: denom ~0.\n";
            return x1;
        }
        double x2 = x1 - f1*(x1 - x0)/(f1 - f0);
        double f2 = func(x2);

        if (std::fabs(f2) < eps)
            return x2;

        x0 = x1;  
        f0 = f1;
        x1 = x2;  
        f1 = f2;
    }
    return x1;
}

// ponto fixo
double pontoFixo(double x0, double eps, int criterioParada, double (*g)(double)){
    double x = x0;
    for (int i = 0; i < criterioParada; i++) {
        double x1 = g(x);
        if (std::fabs(x1 - x) < eps)
            return x1;
        x = x1;
    }
    return x;
}

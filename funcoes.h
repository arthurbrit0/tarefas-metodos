#ifndef FUNCOES_H
#define FUNCOES_H

double fPendulo(double a);

double fPolinomio(double a);

double fLog(double x);

double dfPendulo(double a);  

double dfPolinomio(double a);  

double dfLog(double x);   

double bissecao(double aEsq, double aDir, double eps, int criterioParada, double (*func)(double));

double posicaoFalsa(double aEsq, double aDir, double eps, int criterioParada, double (*func)(double));

double newtonRaphson(double x0, double eps, int criterioParada, double (*func)(double), double (*dfunc)(double));

double secante(double x0, double x1, double eps, int criterioParada, double (*func)(double));

double pontoFixo(double x0, double eps, int criterioParada, double (*g)(double));

double gLog(double x);

bool isolamentoAnalitico(double &aEsq, double &aDir, double step, double (*func)(double));

#endif



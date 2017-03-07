#ifndef __EIGENSOLVER_HPP__
#define __EIGENSOLVER_HPP__

#include <stdlib.h>
#include <assert.h>
#include <error.h>

#include <complex.h>
#include <math.h>

#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>


class Eigensolver
{
private:
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> hermitian_solver;
    Eigen::RowVectorXcd* eengs;
    Eigen::MatrixXcd* evecs;
public:
    // consturctors
    Eigensolver();
    Eigensolver(const unsigned n);
    
    // destructors
    ~Eigensolver();
    
    
    // methods
    void eigensystem(double __complex__* H, double __complex__** E, double __complex__** V, const unsigned n, double* t);
    void eigensystem(double __complex__* H, Eigen::RowVectorXcd & E, Eigen::MatrixXcd & V, const unsigned n, double* t);
};

Eigensolver::Eigensolver()
{
    eengs = NULL;
    evecs = NULL;
    Eigen::initParallel();
}

Eigensolver::Eigensolver(const unsigned n)
{
    eengs = new Eigen::RowVectorXcd(n);
    evecs = new Eigen::MatrixXcd(n,n);
    
    Eigen::initParallel();
}

Eigensolver::~Eigensolver()
{
    if (eengs) delete eengs;
    if (evecs) delete evecs;
}


void Eigensolver::eigensystem(double __complex__* H, double __complex__** E, double __complex__** V, const unsigned n, double* t = NULL)
{
    // map H to Eigen matrix
    Eigen::Map<Eigen::MatrixXcd> h_mat((std::complex<double>*)H,n,n); // https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
    
    
    
    clock_t start = clock();
    hermitian_solver.compute(h_mat);
    clock_t end = clock();
    if (t) *t += (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
    
    *eengs = hermitian_solver.eigenvalues().transpose();
    *evecs = hermitian_solver.eigenvectors();
    
    *E = (double __complex__*) &(*eengs)(0);
    *V = (double __complex__*) &(*evecs)(0,0);
}

void Eigensolver::eigensystem(double __complex__* H, Eigen::RowVectorXcd & E, Eigen::MatrixXcd & V, const unsigned n, double* t=NULL) 
{
    // map H to Eigen matrix
    Eigen::Map<Eigen::MatrixXcd> h_mat((std::complex<double>*)H,n,n); // https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
    
    clock_t start = clock();
    hermitian_solver.compute(h_mat);
    clock_t end = clock();
    if (t) *t += (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
    
    E = hermitian_solver.eigenvalues().transpose();
    V = hermitian_solver.eigenvectors();
}

#endif
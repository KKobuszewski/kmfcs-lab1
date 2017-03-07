#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <unistd.h>

#include <Eigensolver.hpp>

#ifndef M_PI
#define M_PI
#endif


// ==================================================== BASE =============================================================

#ifdef COSSINBASE

inline double __complex__ even_functions(const double x, const double A, const unsigned n)
{
    if (n==0)
        return 1./sqrt(2.*A) + 0.*I;
    else
        return cos(x*((double)n)*M_PI/A)/sqrt(A) + 0.*I;
}

inline double __complex__ odd_functions(const double x, const double A, const unsigned n)
{
    return sin(x*((double)n)*M_PI/A)/sqrt(A) + 0.*I;
}


inline double __complex__ even_funcs_2derivative(const double x, const double A, const unsigned n)
{
    return n*n*M_PI*M_PI*cos(x*((double)n)*M_PI/A)/sqrt(A)/A/A + 0.*I;
}

inline double __complex__ odd_funcs_2derivative(const double x, const double A, const unsigned n)
{
    return n*n*M_PI*M_PI*sin(x*((double)n)*M_PI/A)/sqrt(A)/A/A + 0.*I;
}

#endif


// ==================================================== POTENTIAL =============================================================
#ifdef WELL_POTENTIAL

inline double potential(const double x, const double V0, const double a)
{
    if ( (x > a) || (x < -a) )
        return 0.;
    else
        return -V0;
}

#endif


#ifdef WOODS_SAXON

#ifndef R
#define R 1.
#endif

inline double potential(const double x, const double V0, const double a)
{
    return  - V0 / (1 + exp( (fabs(x)-a)/R ));
}

#endif


// ========================================== HAMILTONIAN MATRIX EVALUATION ====================================================

#define EVEN 0
#define ODD 1

template <int type>
inline double __complex__ hij_on_grid(const double x, const double A, const double V0, const double a, const unsigned i, const unsigned j)
{
    switch(type)
    {
        case EVEN:
        {
            //     conjugate psi_i                   kinetic term psi_j             potential term psi_j
            return conj(even_functions(x,A,i)) * ( even_funcs_2derivative(x,A,j) + potential(x,V0,a)*even_functions(x,A,j) ); // hij(x)
        }
        case ODD:
        {
            //     conjugate psi_i                   kinetic term psi_j             potential term psi_j
            return conj(odd_functions(x,A,i)) * ( odd_funcs_2derivative(x,A,j) + potential(x,V0,a)*odd_functions(x,A,j) ); // hij(x)
        }
    }
}

// TODO: test if really quicker with inline
template <int type>
double __complex__ integrate_hij_on_grid(const unsigned nx, const double xmin, const double xmax, 
                                const double A, const double V0, const double a, const unsigned i, const unsigned j)
{
    const double dx = (xmax - xmin)/nx;
    double __complex__ val = 0. + 0.*I;
    for (unsigned ii =0; ii < nx; ii++)
    {
        // x_i = x_min + i*dx, Euler integration
        val += hij_on_grid<type>(xmin + ii*dx,A,V0,a,i,j);
    }
    
#ifdef DEBUG
    //printf("dx: %lf\n",dx);
    switch(type)
    {
    case EVEN:
    {
        //printf( "(%u,%u)\tx:%lf\t%lf\n",i,j,0., creal(even_functions(0.,A,i)) );
        printf( "even (%u,%u)\tx:%lf\t%lf\n",i,j,0., creal(val*dx) );
        break;
    }
    case ODD:
    {
        printf( "odd (%u,%u)\tx:%lf\t%lf\n",i,j,0., creal(val*dx) );
        break;
    }
    }
#endif
    
    return val*dx;
}


template <int type>
void construct_hamiltonian(double __complex__* hij, unsigned nlevels, const unsigned nx, const double xmin, const double xmax, 
                           const double A, const double V0, const double a)
{
    for (unsigned ii=0; ii<nlevels; ii++)
    for (unsigned jj=0; jj<nlevels; jj++)
    {
        hij[ii + nlevels*jj] = integrate_hij_on_grid<type>(nx, xmin, xmax, A, V0, a, ii, jj);// 2*ii+type, 2*jj+type);
    }
}


void print_hamiltonian(double __complex__* h_sym, double __complex__* h_asym, const unsigned nlevels)
{
    printf("------------------------------------------------------------------------------\n");
    for (unsigned ii=0; ii < nlevels; ii++)
    {
        printf("|  ");
        for (unsigned jj=0; jj < nlevels; jj++)
        {
            double __complex__ vij = 0. + 0.*I;
            if ((ii%2 == 0) && (jj%2 == 0))
            {
                 //printf("\nsym %u, %u, %u\n;",ii,jj,(ii/2)*(nlevels/2) + jj/2);
                 vij = h_sym[(ii/2)*(nlevels/2) + jj/2];
            }
            else if ((ii%2 == 1) && (jj%2 == 1))
            {
                 //printf("\nanym %u, %u, %u\n;",ii,jj,((ii-1)/2)*(nlevels/2) + (jj-1)/2);
                 vij = h_asym[((ii-1)/2)*(nlevels/2) + (jj-1)/2];
            }
            if (creal(vij) >= 0) printf(" ");
            if (cimag(vij) >= 0) printf( "%.2lf+%.1lfj ", creal(vij), fabs(cimag(vij)) );
            else                 printf( "%.2lf-%.1lfj ", creal(vij), -1*(cimag(vij))  );
            //usleep(1000);
        }
        printf("\n");
    }
    printf("------------------------------------------------------------------------------\n");
}


int main(int argc, char* argv[])
{
#ifdef DEBUG
#ifndef NLEVELS
    if (argc < 2) { fprintf(stderr,"Error! no number of levels given!\n"); exit(EXIT_FAILURE); }
    unsigned nlevels = (size_t) atoi(argv[1]);  // set number of levels by cmd line
#else
    const unsigned nlevels = NLEVELS;    // 10-20
#endif
#endif
    const double A  = 15.;
    const double V0 = 2.5;
    const double a  = 10./3.;
    const unsigned nx = 1024;
    const double xmin = -A;
    const double xmax =  A; 
    
#ifndef DEBUG
    for (unsigned nlevels = 8; nlevels <= 100; nlevels += 2)
    {
#endif
        // declare two matrices for 
        //printf( "%u , %u\n", (nlevels/2), (nlevels/2)*(nlevels/2) );
        double __complex__ *h_sym  = (double __complex__*) malloc( (nlevels/2)*(nlevels/2) * sizeof(double __complex__) );
        double __complex__ *h_asym = (double __complex__*) malloc( (nlevels/2)*(nlevels/2) * sizeof(double __complex__) );
        
        // count matrix elements
        clock_t start = clock();
        construct_hamiltonian<EVEN>(h_sym, nlevels/2, nx, xmin, xmax, A, V0, a);
        construct_hamiltonian<ODD>(h_asym, nlevels/2, nx, xmin, xmax, A, V0, a);
        clock_t end = clock();
        double time_matrix = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
#ifdef DEBUG
        // print matrices
        printf("hamiltonian matrix elements:\n");
        print_hamiltonian(h_sym, h_asym, nlevels);
        printf("\n");
        printf("\n");
#endif
        double time_eigs = 0.;
        Eigensolver s;
        double __complex__ *E_sym, *V_sym;
        s.eigensystem(h_sym,&E_sym,&V_sym,nlevels/2,&time_eigs);
    //     double __complex__ *E_asym, *V_asym;
    //     s_asym.eigensystem(h_asym,&E_asym,&V_asym,nlevels/2);
        
#ifdef DEBUG
        printf("Eigenvalues:\n");
#else
        printf("%u\t",nlevels);
#endif
        for (unsigned ii=0; ii < 4; ii++)
        {
            printf("%.3lf\t",creal(E_sym[ii]));
        }
        
        double __complex__ *E_asym, *V_asym;
        s.eigensystem(h_asym,&E_asym,&V_asym,nlevels/2,&time_eigs);
        for (unsigned ii=0; ii < 4; ii++)
        {
            printf("%.3lf\t",creal(E_asym[ii]));
        }
        printf("%lf\t%lf",time_matrix,time_eigs);
        printf("\n");
        
        free(h_sym);
        free(h_asym);
#ifndef DEBUG
    }
#endif
    
    
    return EXIT_SUCCESS;
}

// 2pcf: least squares two-point correlation function code
// ---
// author:  Nicolas Tessore <nicolas.tessore@manchester.ac.uk>
// date:    19 Aug 2018

#define _XOPEN_SOURCE 600

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// LAPACK
extern void dptsv_(int*, int*, double*, double*, double*, int*, int*);

static inline double complex phase(unsigned n, double x, double y)
{
    double complex w = (x + I*y)/(x - I*y);
    double complex z = 1;
    for(; n > 1; n -= 2)
        z *= w;
    if(n == 1)
        z *= (x + I*y)/hypot(x, y);
    return z;
}

static const int DW = 8;

int mapsort(const void* a, const void* b)
{
    const double* x = a;
    const double* y = b;
    
    if(x[7] < y[7])
        return -1;
    if(x[7] > y[7])
        return +1;
    
    if(x[2] < y[2])
        return -1;
    if(x[2] > y[2])
        return +1;
    if(x[1] < y[1])
        return -1;
    if(x[1] > y[1])
        return +1;
    if(x[0] < y[0])
        return -1;
    if(x[0] > y[0])
        return +1;
    
    return 0;
}

static inline int index(double x, double y, double z, double s, int w, int h)
{
    return (int)(z/s)*(w*h) + (int)(y/s)*w + (int)(x/s);
}

static inline int query(int q, const int ma[], int gx, int gy, int gz,
                                                    int gr, int* c, int v[])
{
    int i, il, ih, j, jl, jh, k, kl, kh, l, m, n, p;
    
    i = q/(gx*gy);
    j = (q/gx)%gy;
    k = q%gx;
    
    il = i > gr ? i-gr : 0;
    ih = i+gr < gz ? i+gr+1 : gz;
    
    jl = j > gr ? j-gr : 0;
    jh = j+gr < gy ? j+gr+1 : gy;
    
    kl = k > gr ? k-gr : 0;
    kh = k+gr < gx ? k+gr+1 : gx;
    
    n = 0;
    p = -1;
    
    for(i = il; i < ih; ++i)
    {
        for(j = jl; j < jh; ++j)
        {
            k = (i*gy + j)*gx;
            l = ma[k + kl];
            m = ma[k + kh];
            
            if(l == p)
                p = (v[2*n-1] = m);
            else
                p = (v[2*n+0] = l, v[2*n+1] = m), ++n;
        }
    }
    
    *c = n;
    
    return q;
}

#include "io.c" // yes, really

static const char* ANIM[] = {
    "\033[34m\xe2\xa0\xb7\033[0m", "\033[34m\xe2\xa0\xaf\033[0m",
    "\033[34m\xe2\xa0\x9f\033[0m", "\033[34m\xe2\xa0\xbb\033[0m",
    "\033[34m\xe2\xa0\xbd\033[0m", "\033[34m\xe2\xa0\xbe\033[0m"
};
static const int NANIM = sizeof(ANIM)/sizeof(*ANIM);

volatile sig_atomic_t AL;
volatile sig_atomic_t QQ;

void handler(int s)
{
    AL = (s == SIGALRM);
    QQ = (s == SIGQUIT);
    signal(s, handler);
}

int main(int argc, char* argv[])
{
    char* cfgfile;
    struct config cfg;
    
    bool xc, ls, sc, tc;
    int nd;
    double dl, dh, d0, dm, Dl, Dh;
    double ui, uo;
    int S1, S2;
    
    int n1, n2;
    double* c1;
    double* c2;
    
    double gs;
    int gr;
    double xl, xh, yl, yh, zl, zh;
    int gx, gy, gz, ng;
    int* ma;
    
    double* N;
    double* W;
    double* Y;
    double* A;
    double* X;
    
    time_t st;
    int dt;
    double nn;
    
    int i, j;
    
    char* bf, *nf, *sv, *sx;
    
    if(isatty(fileno(stdout)))
    {
        bf = "\033[1m";
        nf = "\033[0m";
        sv = "\033[32m\xe2\x9c\x94\033[0m";
        sx = "\033[31m\xe2\x9c\x98\033[0m";
    }
    else
    {
        bf = nf = "";
        sv = sx = ">";
    }
    
    cfgfile = NULL;
    memset(&cfg, 0, sizeof(cfg));
    
    if(argc > 5)
        goto err_usage;
    
    if(argc > 1)
    {
        char** posarg[] = {
            NULL,
            &cfgfile,
            &cfg.catalog1,
            &cfg.catalog2,
            &cfg.output
        };
        
        for(int c = 1; c < argc; ++c)
            if(strcmp(argv[c], "--") != 0)
                *posarg[c] = strdup(argv[c]);
    }
    
    if(!cfgfile)
        cfgfile = strdup("2pcf.cfg");
    
    readcfg(cfgfile, &cfg);
    
    printf("\n");
    printf("%sconfiguration file %s%s\n", bf, cfgfile, nf);
    printf("\n");
    printcfg(&cfg);
    printf("\n");
    
    xc = cfg.catalog2 != NULL;
    sc = cfg.coords >= COORDS_LONLAT;
    ls = cfg.spacing == SPACING_LOG;
    
    ui = UCONV[cfg.dunit];
    uo = UCONV[cfg.thunit];
    
    nd = cfg.nth;
    dl = cfg.thmin*uo;
    dh = cfg.thmax*uo;
    
    S1 = cfg.spin1;
    S2 = cfg.spin2;
    
    if(cfg.coords == COORDS_3D && (S1 != 0 || S2 != 0))
    {
        fprintf(stderr, "error: 3D fields must be spin 0\n");
        return EXIT_FAILURE;
    }
    
#ifdef _OPENMP
    if(cfg.num_threads)
        omp_set_num_threads(cfg.num_threads);
    tc = cfg.thread_data == TDATA_COPY;
#else
    tc = false;
#endif
    
    if(sc)
    {
        dl = 2*sin(0.5*dl);
        dh = 2*sin(0.5*dh);
    }
    
    if(ls)
    {
        d0 = log(dl);
        dm = (nd - 1)/(log(dh) - d0);
    }
    else
    {
        d0 = dl;
        dm = (nd - 1)/(dh - d0);
    }
    
    Dl = dl*dl;
    Dh = dh*dh;
    
    printf("%sread catalog%s%s\n", bf, xc ? " 1" : "", nf);
    fflush(stdout);
    
    c1 = readc(cfg.catalog1, cfg.coords, ui, cfg.field1, cfg.signs1, &n1);
    
    printf("%s done with %d points\n", sv, n1);
    printf("\n");
    
    if(xc)
    {
        printf("%sread catalog 2%s\n", bf, nf);
        fflush(stdout);
        
        c2 = readc(cfg.catalog2, cfg.coords, ui, cfg.field2, cfg.signs2, &n2);
        
        printf("%s done with %d points\n", sv, n2);
        printf("\n");
    }
    else
    {
        c2 = c1;
        n2 = n1;
        S2 = S1;
    }
    
    printf("%sbuild index%s\n", bf, nf);
    fflush(stdout);
    
    gs = 0.25*dh;
    gr = ceil(dh/gs);
    
    xl = xh = c1[0];
    yl = yh = c1[1];
    zl = zh = c1[2];
    for(i = 1; i < n1; ++i)
    {
        if(c1[i*DW+0] < xl) xl = c1[i*DW+0];
        if(c1[i*DW+0] > xh) xh = c1[i*DW+0];
        if(c1[i*DW+1] < yl) yl = c1[i*DW+1];
        if(c1[i*DW+1] > yh) yh = c1[i*DW+1];
        if(c1[i*DW+2] < zl) zl = c1[i*DW+2];
        if(c1[i*DW+2] > zh) zh = c1[i*DW+2];
    }
    if(xc)
    {
        for(i = 0; i < n2; ++i)
        {
            if(c2[i*DW+0] < xl) xl = c2[i*DW+0];
            if(c2[i*DW+0] > xh) xh = c2[i*DW+0];
            if(c2[i*DW+1] < yl) yl = c2[i*DW+1];
            if(c2[i*DW+1] > yh) yh = c2[i*DW+1];
            if(c2[i*DW+2] < zl) zl = c2[i*DW+2];
            if(c2[i*DW+2] > zh) zh = c2[i*DW+2];
        }
    }
    
    gx = floor((xh - xl)/gs) + 1;
    gy = floor((yh - yl)/gs) + 1;
    gz = floor((zh - zl)/gs) + 1;
    
    ng = gx*gy*gz;
    
    for(i = 0; i < n1; ++i)
        c1[i*DW+7] =
            index(c1[i*DW+0]-xl, c1[i*DW+1]-yl, c1[i*DW+2]-zl, gs, gx, gy);
    qsort(c1, n1, DW*sizeof(double), mapsort);
    
    if(xc)
    {
        for(i = 0; i < n2; ++i)
            c2[i*DW+7] =
                index(c2[i*DW+0]-xl, c2[i*DW+1]-yl, c2[i*DW+2]-zl, gs, gx, gy);
        qsort(c2, n2, DW*sizeof(double), mapsort);
    }
    
    ma = malloc((ng+1)*sizeof(int));
    if(!ma)
        goto err_alloc;
    
    for(i = 0, j = 0; i < ng; ++i)
    {
        while(j < n2 && c2[j*DW+7] < i)
            j += 1;
        ma[i] = j;
    }
    ma[ng] = n2;
    
    printf("%s done with %d x %d x %d grid cells\n", sv, gx, gy, gz);
    printf("\n");
    
    N = calloc(nd, sizeof(double));
    W = calloc(nd*2, sizeof(double));
    Y = calloc(nd*4, sizeof(double));
    if(!N || !W || !Y)
        goto err_alloc;
    
    signal(SIGALRM, handler);
    signal(SIGQUIT, handler);
    AL = QQ = 0;
    
    printf("%scorrelations%s\n", bf, nf);
    fflush(stdout);
    
    st = time(NULL);
    dt = 0;
    
    #pragma omp parallel default(none) shared(st, dt, N, W, Y, AL, QQ) \
        private(i, j) firstprivate(xc, ls, sc, tc, nd, d0, dm, Dl, Dh, gr, \
            gx, gy, gz, ng, n1, n2, c1, c2, ma, S1, S2, ANIM, NANIM, stdout)
    {
        int qc, jh;
        int q, nq;
        int* qr;
        
        double* c1_;
        double* c2_;
        int* ma_;
        
        double* N_;
        double* W_;
        double* Y_;
        
        bool fb;
        
        nq = 0;
        qr = malloc((2*gr+1)*(2*gr+1)*2*sizeof(int));
        if(!qr)
            perror(NULL), abort();
        
        if(tc)
        {
            c1_ = malloc(n1*DW*sizeof(double));
            if(xc)
                c2_ = malloc(n2*DW*sizeof(double));
            else
                c2_ = c1_;
            ma_ = malloc((ng+1)*sizeof(int));
            if(!c1_ || !c2_ || !ma_)
                perror(NULL), abort();
            
            memcpy(c1_, c1, n1*DW*sizeof(double));
            if(xc)
                memcpy(c2_, c2, n2*DW*sizeof(double));
            memcpy(ma_, ma, (ng+1)*sizeof(int));
        }
        else
        {
            c1_ = c1;
            c2_ = c2;
            ma_ = ma;
        }
        
        N_ = calloc(nd, sizeof(double));
        W_ = calloc(nd*2, sizeof(double));
        Y_ = calloc(nd*4, sizeof(double));
        if(!N_ || !W_ || !Y_)
            perror(NULL), abort();
        
        fb = false;
        
        #pragma omp master
        if(isatty(fileno(stdout)))
        {
            fb = true;
            AL = false;
            alarm(1);
            
#ifdef _OPENMP
            printf("\r%s %d thread(s) ", ANIM[0], omp_get_num_threads());
            fflush(stdout);
#endif
        }
        
        qc = -1;
        
        #pragma omp for schedule(static, 1) nowait
        for(i = 0; i < n1; ++i)
        {
            const double xi = c1_[i*DW+0];
            const double yi = c1_[i*DW+1];
            const double zi = c1_[i*DW+2];
            const double ui = c1_[i*DW+3];
            const double vi = c1_[i*DW+4];
            const double wi = c1_[i*DW+5];
            const int    qi = c1_[i*DW+7];
            
            if(QQ)
                continue;
            
            if(AL && fb)
            {
                dt = difftime(time(NULL), st);
                
                printf("\r%s %.2f%%", ANIM[dt%NANIM], 100.*i/n1);
                printf(" in %02d:%02d:%02d ", dt/3600, (dt/60)%60, dt%60);
                fflush(stdout);
                
                AL = false;
                alarm(1);
            }
            
            if(qi != qc)
                qc = query(qi, ma, gx, gy, gz, gr, &nq, qr);
            
            for(q = 0; q < nq; ++q)
            {
                j = qr[2*q+0];
                jh = qr[2*q+1];
                
                if(!xc && j < i+1)
                    j = i+1;
                
                for(; j < jh; ++j)
                {
                    const double xj = c2_[j*DW+0];
                    const double yj = c2_[j*DW+1];
                    const double zj = c2_[j*DW+2];
                    const double uj = c2_[j*DW+3];
                    const double vj = c2_[j*DW+4];
                    const double wj = c2_[j*DW+5];
                    
                    const double dx = xi - xj;
                    const double dy = yi - yj;
                    const double dz = zi - zj;
                    
                    const double D = dx*dx + dy*dy + dz*dz;
                    
                    if(D >= Dl && D < Dh)
                    {
                        double complex gi, gj, xip, xim;
                        
                        double fl, fh, ww;
                        int nl, nh;
                        
                        const double dp = xi*xj + sc*(yi*yj+zi*zj);
                        const double cp = xj*yi - yj*xi;
                        
                        fl = dm*((ls ? 0.5*log(D) : sqrt(D)) - d0);
                        nl = floor(fl);
                        nh = nl + 1;
                        fl = nh - fl;
                        fh = 1 - fl;
                        
                        ww = wi*wj;
                        
                        gi = (ui + I*vi)*phase(S1, zi*dp - zj, cp);
                        gj = (uj + I*vj)*phase(S2, zi - zj*dp, cp);
                        
                        xip = gi*conj(gj);
                        xim = gi*gj;
                        
                        N_[nl] += fl;
                        N_[nh] += fh;
                        
                        W_[0*nd+nl] += ww*fl*fl;
                        W_[0*nd+nh] += ww*fh*fh;
                        W_[1*nd+nl] += ww*fl*fh;
                        
                        Y_[0*nd+nl] += ww*fl*creal(xip);
                        Y_[0*nd+nh] += ww*fh*creal(xip);
                        Y_[1*nd+nl] += ww*fl*creal(xim);
                        Y_[1*nd+nh] += ww*fh*creal(xim);
                        Y_[2*nd+nl] += ww*fl*cimag(xip);
                        Y_[2*nd+nh] += ww*fh*cimag(xip);
                        Y_[3*nd+nl] += ww*fl*cimag(xim);
                        Y_[3*nd+nh] += ww*fh*cimag(xim);
                    }
                }
            }
        }
        
        #pragma omp critical
        {
            for(int n = 0; n < nd; ++n)
            {
                N[n] += N_[n];
                
                W[0*nd+n] += W_[0*nd+n];
                W[1*nd+n] += W_[1*nd+n];
                
                Y[0*nd+n] += Y_[0*nd+n];
                Y[1*nd+n] += Y_[1*nd+n];
                Y[2*nd+n] += Y_[2*nd+n];
                Y[3*nd+n] += Y_[3*nd+n];
            }
        }
        
        free(qr);
        free(N_);
        free(W_);
        free(Y_);
        
        if(tc)
        {
            free(c1_);
            if(xc)
                free(c2_);
            free(ma_);
        }
    }
    
    nn = 0;
    for(int n = 0; n < nd; ++n)
        nn += N[n];
    
    dt = difftime(time(NULL), st);
    
    if(isatty(fileno(stdin)))
        printf("\r");
    printf("%s done with %.0f pairs", sv, nn);
    printf(" in %02d:%02d:%02d  \n", dt/3600, (dt/60)%60, dt%60);
    printf("\n");
    
    A = calloc(nd*2, sizeof(double));
    X = calloc(nd*4, sizeof(double));
    if(!A || !X)
        goto err_alloc;
    
    memcpy(A, W, nd*2*sizeof(double));
    memcpy(X, Y, nd*4*sizeof(double));
    
    // solve A.X = Y
    {
        int n = nd, m = 4, err;
        
        printf("%sleast squares%s\n", bf, nf);
        fflush(stdout);
        
        dptsv_(&n, &m, A, A+n, X, &n, &err);
        
        if(!err)
            printf("%s success\n", sv);
        else
        {
            printf("%s error: ", sx);
            if(err > 0)
                printf("the %dx%d submatrix is not pos. def.", err, err);
            else
                printf("illegal argument");
            printf("\n");
        }
        printf("\n");
    }
    
    output(cfg.output, nd, dl, dh, uo, sc, ls, N, X, W, Y);
    
    free(N);
    free(W);
    free(Y);
    free(A);
    free(X);
    free(ma);
    free(c1);
    if(xc)
        free(c2);
    
    free(cfgfile);
    freecfg(&cfg);
    
    return EXIT_SUCCESS;
    
err_usage:
    fprintf(stderr, "usage: 2pcf [config] [catalog] [catalog2] [output]\n");
    return EXIT_FAILURE;
    
err_alloc:
    perror(NULL);
    return EXIT_FAILURE;
}

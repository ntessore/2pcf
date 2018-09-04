// 2pcf: least squares two-point correlation function code
// ---
// author:  Nicolas Tessore <nicolas.tessore@manchester.ac.uk>
// date:    19 Aug 2018

#define _XOPEN_SOURCE 600

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// LAPACK
extern void dposv_(char*, int*, int*, double*, int*, double*, int*, int*);

static const double PI_HALF = 1.5707963267948966192;
static const double TWO_PI = 6.2831853071795864769;

static inline void cmul(double x, double y, double u, double v,
                                                    double* re, double* im)
{
    const double k1 = x*(u + v);
    const double k2 = v*(x + y);
    const double k3 = u*(y - x);
    *re = k1 - k2;
    *im = k1 + k3;
}

static inline void sincos(double x, double* s, double* c)
{
    *s = sin(x);
    *c = cos(x);
}

static inline void nsincos(int n, double x, double y, double* s, double* c)
{
    double h;
    
    switch(n)
    {
    case 0:
        *c = 1;
        *s = 0;
        return;
        
    case 1:
        h = hypot(x, y);
        *c = x/h;
        *s = y/h;
        return;
        
    case 2:
        h = x*x + y*y;
        *c = (x*x - y*y)/h;
        *s = (2*x*y)/h;
        return;
        
    default:
        h = hypot(x, y);
        x = x/h;
        y = y/h;
        *c = x*x - y*y;
        *s = 2*x*y;
        for(int i = 2; i < n; ++i)
            cmul(x, y, *c, *s, c, s);
        return;
    }
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
    
    if(x[0] < y[0])
        return -1;
    if(x[0] > y[0])
        return +1;
    
    return 0;
}

static inline int index(double x, double y, double dx, double dy, int w)
{
    return (int)(y/dy)*w + (int)(x/dx);
}

static inline void query(int k, int w, int h, int dy, const int dx[],
                                                        int* qc, int qv[])
{
    const int i0 = k/w;
    const int il = i0 > dy ? i0-dy : 0;
    const int ih = i0+dy < h ? i0+dy+1 : h;
    
    const int j = k%w;
    
    int n = 0, qq = -1;
    
    for(int i = il; i < ih; ++i)
    {
        const int di = dx[i];
        if(di < 0)
        {
            const int q0 = i*w;
            const int qw = q0 + w;
            const int ql = q0 + (j+w+di)%w;
            const int qh = q0 + (j-di+1)%w;
            
            if(ql < qh)
            {
                if(ql == qq)
                    qq = (qv[2*n-1] = qh);
                else
                    qq = (qv[2*n+0] = ql, qv[2*n+1] = qh), ++n;
            }
            else
            {
                if(q0 == qq)
                    qq = (qv[2*n-1] = qh);
                else
                    qq = (qv[2*n+0] = q0, qv[2*n+1] = qh), ++n;
                
                if(ql == qq)
                    qq = (qv[2*n-1] = qw);
                else
                    qq = (qv[2*n+0] = ql, qv[2*n+1] = qw), ++n;
            }
        }
        else
        {
            const int q0 = i*w;
            const int ql = q0 + (j > di ? j-di : 0);
            const int qh = q0 + (j+di < w ? j+di+1 : w);
            
            if(ql == qq)
                qq = (qv[2*n-1] = qh);
            else
                qq = (qv[2*n+0] = ql, qv[2*n+1] = qh), ++n;
        }
    }
    
    *qc = n;
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
    const char* cfgfile;
    struct config cfg;
    
    bool xc, ls, sc, tc;
    int nd;
    double dl, dh, sdh, Dl, Dh, D0, Dm;
    double ui, uo;
    int S1, S2;
    
    size_t n1, n2;
    double* c1;
    double* c2;
    
    double xl, xh, yl, yh;
    size_t gw, gh, ng;
    double gx, gy;
    int dy;
    int* dx;
    int* ma;
    
    double* N;
    double* W;
    double* X;
    
    time_t st;
    int dt;
    double nn;
    
    size_t i, j;
    
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
    
    if(argc > 4)
        goto err_usage;
    
    cfgfile = argc > 1 ? argv[1] : "2pcf.cfg";
    readcfg(cfgfile, &cfg);
    
    if(argc > 2)
        free(cfg.catalog1), cfg.catalog1 = copystr(argv[2]);
    if(argc > 3)
        free(cfg.catalog2), cfg.catalog2 = copystr(argv[3]);
    
    xc = cfg.catalog2 != NULL;
    sc = cfg.coords != COORDS_FLAT;
    ls = cfg.spacing == SPACING_LOG;
    
#ifdef _OPENMP
    if(cfg.num_threads)
        omp_set_num_threads(cfg.num_threads);
    tc = cfg.thread_data == TDATA_COPY;
#else
    tc = false;
#endif
    
    printf("\n");
    printf("%sconfiguration file %s%s\n", bf, cfgfile, nf);
    printf("\n");
    printcfg(&cfg);
    printf("\n");
    
    ui = UCONV[cfg.dunit];
    uo = UCONV[cfg.thunit];
    
    nd = cfg.nth;
    dl = cfg.thmin*uo;
    dh = cfg.thmax*uo;
    
    sdh = sin(dh);
    
    Dl = dl;
    Dh = dh;
    if(sc)
    {
        Dl = 2*sin(0.5*Dl);
        Dh = 2*sin(0.5*Dh);
    }
    Dl = Dl*Dl;
    Dh = Dh*Dh;
    
    if(ls)
    {
        D0 = log(Dl);
        Dm = (nd - 1)/(log(Dh) - D0);
    }
    else
    {
        D0 = Dl;
        Dm = (nd - 1)/(Dh - D0);
    }
    
    S1 = cfg.spin1;
    S2 = cfg.spin2;
    
    printf("%sreading catalog%s%s\n", bf, xc ? " 1" : "", nf);
    fflush(stdout);
    
    c1 = readc(cfg.catalog1, ui, cfg.field1, cfg.signs1, &n1);
    
    printf("%s done with %zu points\n", sv, n1);
    printf("\n");
    
    if(xc)
    {
        printf("%sreading catalog 2%s\n", bf, nf);
        fflush(stdout);
        
        c2 = readc(cfg.catalog2, ui, cfg.field2, cfg.signs2, &n2);
        
        printf("%s done with %zu points\n", sv, n2);
        printf("\n");
    }
    else
    {
        c2 = c1;
        n2 = n1;
        S2 = S1;
    }
    
    printf("%sbuilding index%s\n", bf, nf);
    fflush(stdout);
    
    gx = cfg.gridx*uo;
    gy = cfg.gridy*uo;
    
    if(sc)
    {
        xl = 0;
        xh = TWO_PI;
        yl = -PI_HALF;
        yh = +PI_HALF;
        
        gw = (int)(fmax(1, floor((xh - xl)/gx))) | 1;
        gh = fmax(1, floor((yh - yl)/gy));
        
        gx = (xh - xl)/gw;
        gy = (yh - yl)/gh;
    }
    else
    {
        xl = xh = c1[0];
        yl = yh = c1[1];
        for(i = 1; i < n1; ++i)
        {
            if(c1[i*DW+0] < xl) xl = c1[i*DW+0];
            if(c1[i*DW+0] > xh) xh = c1[i*DW+0];
            if(c1[i*DW+1] < yl) yl = c1[i*DW+1];
            if(c1[i*DW+1] > yh) yh = c1[i*DW+1];
        }
        if(xc)
        {
            for(i = 0; i < n2; ++i)
            {
                if(c2[i*DW+0] < xl) xl = c2[i*DW+0];
                if(c2[i*DW+0] > xh) xh = c2[i*DW+0];
                if(c2[i*DW+1] < yl) yl = c2[i*DW+1];
                if(c2[i*DW+1] > yh) yh = c2[i*DW+1];
            }
        }
        
        gw = floor((xh - xl)/gx) + 1;
        gh = floor((yh - yl)/gy) + 1;
    }
    
    ng = gw*gh;
    
    dy = ceil(dh/gy);
    dx = malloc(gh*sizeof(int));
    if(!dx)
        goto err_alloc;
    
    if(sc)
    {
        for(i = 0; i < gh; ++i)
        {
            const double cy = fmin(cos(yl + i*gy), cos(yl + (i+1)*gy));
            const int di = sdh > cy ? gw/2 : ceil(asin(sdh/cy)/gx);
            dx[i] = -di;
        }
    }
    else
    {
        for(i = 0; i < gh; ++i)
            dx[i] = ceil(dh/gx);
    }
    
    for(i = 0; i < n1; ++i)
        c1[i*DW+7] = index(c1[i*DW+0] - xl, c1[i*DW+1] - yl, gx, gy, gw);
    qsort(c1, n1, DW*sizeof(double), mapsort);
    
    if(xc)
    {
        for(i = 0; i < n2; ++i)
            c2[i*DW+7] = index(c2[i*DW+0] - xl, c2[i*DW+1] - yl, gx, gy, gw);
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
    
    printf("%s done with %zu x %zu grid cells\n", sv, gw, gh);
    printf("\n");
    
    if(sc)
    {
        for(i = 0; i < n1; ++i)
        {
            sincos(c1[i*DW+0], &c1[i*DW+0], &c1[i*DW+2]);
            sincos(c1[i*DW+1], &c1[i*DW+1], &c1[i*DW+3]);
        }
        if(xc)
        {
            for(i = 0; i < n2; ++i)
            {
                sincos(c2[i*DW+0], &c2[i*DW+0], &c2[i*DW+2]);
                sincos(c2[i*DW+1], &c2[i*DW+1], &c2[i*DW+3]);
            }
        }
    }
    
    N = calloc(nd, sizeof(double));
    W = calloc(nd*nd, sizeof(double));
    X = calloc(nd*4, sizeof(double));
    if(!N || !W || !X)
        goto err_alloc;
    
    signal(SIGALRM, handler);
    signal(SIGQUIT, handler);
    AL = QQ = 0;
    
    printf("%scalculating correlations%s\n", bf, nf);
    fflush(stdout);
    
    st = time(NULL);
    dt = 0;
    
    #pragma omp parallel default(none) shared(st, dt, N, W, X, AL, QQ) \
        private(i, j) firstprivate(xc, ls, sc, tc, nd, Dl, Dh, D0, Dm, ng, \
            gw, gh, dx, dy, n1, n2, c1, c2, ma, S1, S2, ANIM, NANIM, stdout)
    {
        size_t qc, jh;
        int q, nq;
        int* qr;
        
        double* c1_;
        double* c2_;
        int* ma_;
        
        double* N_;
        double* W_;
        double* X_;
        
        bool fb;
        
        nq = 0;
        qr = malloc((2*dy+1)*4*sizeof(int));
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
        W_ = calloc(nd*nd, sizeof(double));
        X_ = calloc(nd*4, sizeof(double));
        if(!N_ || !W_ || !X_)
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
            const double sxi = c1_[i*DW+0];
            const double syi = c1_[i*DW+1];
            const double cxi = c1_[i*DW+2];
            const double cyi = c1_[i*DW+3];
            const double ui  = c1_[i*DW+4];
            const double vi  = c1_[i*DW+5];
            const double wi  = c1_[i*DW+6];
            const size_t qi  = c1_[i*DW+7];
            
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
            {
                qc = qi;
                query(qi, gw, gh, dy, dx, &nq, qr);
                for(q = 0; q < 2*nq; ++q)
                    qr[q] = ma_[qr[q]];
            }
            
            for(q = 0; q < nq; ++q)
            {
                j = qr[2*q+0];
                jh = qr[2*q+1];
                
                if(!xc && j < i+1)
                    j = i+1;
                
                for(; j < jh; ++j)
                {
                    const double sxj = c2_[j*DW+0];
                    const double syj = c2_[j*DW+1];
                    const double cxj = c2_[j*DW+2];
                    const double cyj = c2_[j*DW+3];
                    const double uj  = c2_[j*DW+4];
                    const double vj  = c2_[j*DW+5];
                    const double wj  = c2_[j*DW+6];
                    
                    const double sdx = cxi*sxj - sxi*cxj;
                    const double cdx = cxi*cxj + sc*sxi*sxj;
                    
                    const double D1 = cdx*cyi - cyj;
                    const double D2 = sdx*cyi;
                    const double D3 = syj - syi;
                    
                    const double D = D1*D1 + D2*D2 + D3*D3;
                    
                    if(D >= Dl && D < Dh)
                    {
                        double ai, bi, aj, bj;
                        double xip_re, xip_im, xim_re, xim_im;
                        
                        double fl, fh, ww;
                        int nl, nh;
                        
                        fl = Dm*((ls ? log(D) : D) - D0);
                        nl = floor(fl);
                        nh = nl + 1;
                        fl = nh - fl;
                        fh = 1 - fl;
                        
                        ww = wi*wj;
                        
                        // e^{I phi_ij} unnormalised
                        ai = cyi*syj - syi*cyj*cdx;
                        bi = -cyj*sdx;
                        // e^{I S1 phi_ij}
                        nsincos(S1, ai, bi, &bi, &ai);
                        // ai + I bi = (ui + I vi) e^{-I S1 phi_ij}
                        cmul(ui, vi, ai, -bi, &ai, &bi);
                        
                        // e^{I phi_ji} unnormalised
                        aj = cyj*syi - syj*cyi*cdx;
                        bj = cyi*sdx;
                        // e^{I S2 phi_ji}
                        nsincos(S2, aj, bj, &bj, &aj);
                        // aj + I bj = (uj + I vj) e^{-I S2 phi_ji}
                        cmul(uj, vj, aj, -bj, &aj, &bj);
                        
                        // xip = (ai + I bi) (aj - I bj)
                        cmul(ai, bi, aj, -bj, &xip_re, &xip_im);
                        // xim = (ai + I bi) (aj + I bj)
                        cmul(ai, bi, aj, bj, &xim_re, &xim_im);
                        
                        N_[nl] += fl;
                        N_[nh] += fh;
                        
                        W_[nl*nd+nl] += ww*fl*fl;
                        W_[nl*nd+nh] += ww*fl*fh;
                        W_[nh*nd+nh] += ww*fh*fh;
                        
                        X_[0*nd+nl] += ww*fl*xip_re;
                        X_[0*nd+nh] += ww*fh*xip_re;
                        
                        X_[1*nd+nl] += ww*fl*xim_re;
                        X_[1*nd+nh] += ww*fh*xim_re;
                        
                        X_[2*nd+nl] += ww*fl*xip_im;
                        X_[2*nd+nh] += ww*fh*xip_im;
                        
                        X_[3*nd+nl] += ww*fl*xim_im;
                        X_[3*nd+nh] += ww*fh*xim_im;
                    }
                }
            }
        }
        
        #pragma omp critical
        {
            for(int n = 0; n < nd; ++n)
            {
                N[n] += N_[n];
                
                for(int m = n; m < nd; ++m)
                    W[n*nd+m] += W_[n*nd+m];
                
                X[0*nd+n] += X_[0*nd+n];
                X[1*nd+n] += X_[1*nd+n];
                X[2*nd+n] += X_[2*nd+n];
                X[3*nd+n] += X_[3*nd+n];
            }
        }
        
        free(qr);
        free(N_);
        free(W_);
        free(X_);
        
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
    
    for(int n = 0; n < nd; ++n)
        for(int m = 0; m < n; ++m)
            W[n*nd+m] = W[m*nd+n];
    
    dt = difftime(time(NULL), st);
    
    if(isatty(fileno(stdin)))
        printf("\r");
    printf("%s done with %.0f pairs", sv, nn);
    printf(" in %02d:%02d:%02d  \n", dt/3600, (dt/60)%60, dt%60);
    printf("\n");
    
    if(cfg.matrix)
        writetxt(cfg.matrix, nd, nd, W);
    if(cfg.rhs)
        writetxt(cfg.rhs, nd, 4, X);
    
    // solve A.X = B
    {
        int n = nd, m = 4, err;
        
        printf("%ssolving normal equations%s\n", bf, nf);
        fflush(stdout);
        
        dposv_("L", &n, &m, W, &n, X, &n, &err);
        
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
    
    writexi(cfg.output, nd, Dl, Dh, uo, sc, ls, N, X);
    
    free(N);
    free(W);
    free(X);
    free(ma);
    free(c1);
    if(xc)
        free(c2);
    free(dx);
    
    freecfg(&cfg);
    
    return EXIT_SUCCESS;
    
err_usage:
    fprintf(stderr, "usage: 2pcf [config] [catalog] [catalog2]\n");
    return EXIT_FAILURE;
    
err_alloc:
    perror(NULL);
    return EXIT_FAILURE;
}

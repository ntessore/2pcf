// 2pcf: two-point correlation function code
// ---
// author:  Nicolas Tessore <nicolas.tessore@manchester.ac.uk>
// date:    08 Aug 2018

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

volatile sig_atomic_t fb;
volatile sig_atomic_t qQ;

#pragma omp threadprivate(fb)

void handler(int s)
{
    fb = (s == SIGALRM);
    qQ = (s == SIGQUIT);
}

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

static inline size_t upper_bound(double x, const double v[], size_t n)
{
    size_t i = 0;
    while(n > 0)
    {
        const size_t m = n/2;
        const size_t j = i+n-m;
        const double u = v[i+m];
        n = m;
        i = x < u ? i : j;
    }
    return i;
}

#include "io.c" // yes, really

int main(int argc, char* argv[])
{
    const char* cfgfile;
    struct config cfg;
    
    bool pt, xc, ls, sc, tc;
    size_t nd;
    double dl, dh, sdh;
    double* db;
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
    
    int* m1;
    int* m2;
    
    double* N;
    double* W;
    double* X;
    
    int p, np;
    time_t st;
    int dt;
    size_t ni, nj;
    double* ci;
    double* cj;
    int* mj;
    int Si, Sj;
    double nn;
    
    size_t i, j;
    
    if(argc > 2)
        goto err_usage;
    
    cfgfile = argc > 1 ? argv[1] : "2pcf.cfg";
    readcfg(cfgfile, &cfg);
    printcfg(cfgfile, &cfg);
    
    pt = cfg.mode == MODE_POINTS;
    xc = cfg.catalog2 != NULL;
    sc = cfg.coords != COORDS_FLAT;
    ls = cfg.spacing == SPACING_LOG;
    
#ifdef _OPENMP
    tc = true;
#else
    tc = false;
#endif
    
    if(pt && !xc)
    {
        fprintf(stderr, "error: point mode requires `catalog2`\n");
        return EXIT_FAILURE;
    }
    
    ui = UCONV[cfg.dunit];
    uo = UCONV[cfg.thunit];
    
    nd = cfg.nth;
    dl = cfg.thmin*uo;
    dh = cfg.thmax*uo;
    
    sdh = sin(dh);
    
    db = malloc((nd+1)*sizeof(double));
    if(!db)
        goto err_alloc;
    
    for(i = 0; i <= nd; ++i)
    {
        if(ls)
            db[i] = exp(log(dl) + i*(log(dh) - log(dl))/nd);
        else
            db[i] = dl + i*(dh - dl)/nd;
        if(sc)
            db[i] = 2*sin(0.5*db[i]);
        db[i] = db[i]*db[i];
    }
    
    S1 = cfg.spin1;
    S2 = cfg.spin2;
    
    printf("reading %s\n", pt ? "data catalog" : xc ? "catalog 1" : "catalog");
    
    c1 = readc(cfg.catalog1, ui, pt, cfg.field1, cfg.signs1, &n1);
    
    printf("> done with %zu points\n", n1);
    printf("\n");
    
    if(xc)
    {
        printf("reading %s\n", pt ? "random catalog" : "catalog 2");
        
        c2 = readc(cfg.catalog2, ui, pt, cfg.field2, cfg.signs2, &n2);
        
        printf("> done with %zu points\n", n2);
        printf("\n");
    }
    else
    {
        c2 = NULL;
        n2 = 0;
    }
    
    printf("building index\n");
    
    if(sc)
    {
        xl = 0;
        xh = TWO_PI;
        yl = -PI_HALF;
        yh = +PI_HALF;
        
        gx = cfg.gridx*UCONV[UNIT_DEG];
        gy = cfg.gridy*UCONV[UNIT_DEG];
        
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
        for(i = 0; i < n2; ++i)
        {
            if(c2[i*DW+0] < xl) xl = c2[i*DW+0];
            if(c2[i*DW+0] > xh) xh = c2[i*DW+0];
            if(c2[i*DW+1] < yl) yl = c2[i*DW+1];
            if(c2[i*DW+1] > yh) yh = c2[i*DW+1];
        }
        
        gw = floor((xh - xl)/dh) + 1;
        gh = floor((yh - yl)/dh) + 1;
        
        gx = dh;
        gy = dh;
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
    
    m1 = malloc((ng+1)*sizeof(int));
    if(!m1)
        goto err_alloc;
    
    for(i = 0; i < n1; ++i)
        c1[i*DW+7] = index(c1[i*DW+0] - xl, c1[i*DW+1] - yl, gx, gy, gw);
    
    qsort(c1, n1, DW*sizeof(double), mapsort);
    
    for(i = 0, j = 0; i < ng; ++i)
    {
        while(j < n1 && c1[j*DW+7] < i)
            j += 1;
        m1[i] = j;
    }
    m1[ng] = n1;
    
    if(xc)
    {
        m2 = malloc((ng+1)*sizeof(int));
        if(!m2)
            goto err_alloc;
        
        for(i = 0; i < n2; ++i)
            c2[i*DW+7] = index(c2[i*DW+0] - xl, c2[i*DW+1] - yl, gx, gy, gw);
        
        qsort(c2, n2, DW*sizeof(double), mapsort);
        
        for(i = 0, j = 0; i < ng; ++i)
        {
            while(j < n2 && c2[j*DW+7] < i)
                j += 1;
            m2[i] = j;
        }
        m2[ng] = n2;
    }
    else
        m2 = NULL;
    
    printf("> done with %zu x %zu grid cells\n", gw, gh);
    printf("\n");
    
    if(sc)
    {
        for(i = 0; i < n1; ++i)
        {
            __sincos(c1[i*DW+0], &c1[i*DW+0], &c1[i*DW+2]);
            __sincos(c1[i*DW+1], &c1[i*DW+1], &c1[i*DW+3]);
        }
        for(i = 0; i < n2; ++i)
        {
            __sincos(c2[i*DW+0], &c2[i*DW+0], &c2[i*DW+2]);
            __sincos(c2[i*DW+1], &c2[i*DW+1], &c2[i*DW+3]);
        }
    }
    
    if(pt)
    {
        np = 3;
        N = calloc(nd*np, sizeof(double));
        W = calloc(nd*np, sizeof(double));
        X = NULL;
        if(!N || !W)
            goto err_alloc;
    }
    else
    {
        np = 1;
        N = calloc(nd, sizeof(double));
        W = calloc(nd, sizeof(double));
        X = calloc(nd*4, sizeof(double));
        if(!N || !W || !X)
            goto err_alloc;
    }
    
    signal(SIGQUIT, handler);
    qQ = 0;
    
    for(p = 0; p < np && !qQ; ++p)
    {
        if(pt)
        {
            Si = Sj = 0;
            
            switch(p)
            {
            case 0:
                printf("calculating DD correlations\n");
                
                xc = false;
                ci = c1, ni = n1;
                cj = c1, nj = n1, mj = m1;
                
                break;
                
            case 1:
                printf("calculating DR correlations\n");
                
                xc = true;
                ci = c1, ni = n1;
                cj = c2, nj = n2, mj = m2;
                
                break;
                
            case 2:
                printf("calculating RR correlations\n");
                
                xc = false;
                ci = c2, ni = n2;
                cj = c2, nj = n2, mj = m2;
                
                break;
            }
        }
        else
        {
            printf("calculating correlations\n");
            
            ci = c1, ni = n1, Si = S1;
            if(xc)
                cj = c2, nj = n2, mj = m2, Sj = S2;
            else
                cj = c1, nj = n1, mj = m1, Sj = S1;
        }
        
        st = time(NULL);
        dt = 0;
        fb = 0;
        
        #pragma omp parallel default(none) shared(st, dt, nn, N, W, X, qQ) \
            private(i, j) firstprivate(pt, xc, sc, tc, nd, db, ng, gw, gh, \
                dx, dy, p, ni, nj, ci, cj, mj, Si, Sj, stdout)
        {
            size_t qc;
            int q, nq;
            int* qr;
            
            double* db_;
            double* ci_;
            double* cj_;
            int* mj_;
            
            double* N_;
            double* W_;
            double* X_;
            
            nq = 0;
            qr = malloc((2*dy+1)*4*sizeof(int));
            if(!qr)
                perror(NULL), abort();
            
            if(tc)
            {
                db_ = malloc((nd+1)*sizeof(double));
                ci_ = malloc(ni*DW*sizeof(double));
                if(cj != ci)
                    cj_ = malloc(nj*DW*sizeof(double));
                else
                    cj_ = ci_;
                mj_ = malloc((ng+1)*sizeof(int));
                if(!db_ || !ci_ || !cj_ || !mj_)
                    perror(NULL), abort();
                
                memcpy(db_, db, (nd+1)*sizeof(double));
                memcpy(ci_, ci, ni*DW*sizeof(double));
                if(cj != ci)
                    memcpy(cj_, cj, nj*DW*sizeof(double));
                memcpy(mj_, mj, (ng+1)*sizeof(int));
            }
            else
            {
                db_ = db;
                ci_ = ci;
                cj_ = cj;
                mj_ = mj;
            }
            
            N_ = calloc(nd, sizeof(double));
            W_ = calloc(nd, sizeof(double));
            X_ = X ? calloc(nd*4, sizeof(double)) : NULL;
            if(!N_ || !W_ || (X && !X_))
                perror(NULL), abort();
            
            #pragma omp master
            {
                signal(SIGALRM, handler);
                alarm(1);
            }
            
            qc = -1;
            
            #pragma omp for schedule(static, 1) nowait
            for(i = 0; i < ni; ++i)
            {
                if(qQ)
                    continue;
                
                if(fb)
                {
                    dt = difftime(time(NULL), st);
                    
                    printf("\r> %.2f%%", 100.*i/ni);
                    printf(" - %02d:%02d:%02d ", dt/3600, (dt/60)%60, dt%60);
                    fflush(stdout);
                    
                    fb = 0;
                    alarm(1);
                }
                
                const double sxi = ci_[i*DW+0];
                const double syi = ci_[i*DW+1];
                const double cxi = ci_[i*DW+2];
                const double cyi = ci_[i*DW+3];
                const double ui  = ci_[i*DW+4];
                const double vi  = ci_[i*DW+5];
                const double wi  = ci_[i*DW+6];
                const size_t qi  = ci_[i*DW+7];
                
                if(qi != qc)
                {
                    qc = qi;
                    query(qi, gw, gh, dy, dx, &nq, qr);
                    for(q = 0; q < 2*nq; ++q)
                        qr[q] = mj_[qr[q]];
                }
                
                for(q = 0; q < nq; ++q)
                {
                    j = qr[2*q+0];
                    nj = qr[2*q+1];
                    
                    if(!xc && j < i+1)
                        j = i+1;
                    
                    for(; j < nj; ++j)
                    {
                        const double sxj = cj_[j*DW+0];
                        const double syj = cj_[j*DW+1];
                        const double cxj = cj_[j*DW+2];
                        const double cyj = cj_[j*DW+3];
                        const double uj  = cj_[j*DW+4];
                        const double vj  = cj_[j*DW+5];
                        const double wj  = cj_[j*DW+6];
                        
                        const double sdx = cxi*sxj - sxi*cxj;
                        const double cdx = cxi*cxj + sc*sxi*sxj;
                        
                        const double d1 = cdx*cyi - cyj;
                        const double d2 = sdx*cyi;
                        const double d3 = syj - syi;
                        
                        const double d = d1*d1 + d2*d2 + d3*d3;
                        
                        if(d >= db_[0] && d < db_[nd])
                        {
                            const size_t n = upper_bound(d, db+1, nd-1);
                            
                            const double ww = wi*wj;
                            
                            N_[n] += 1;
                            W_[n] += ww;
                            
                            if(!pt)
                            {
                                double sij, cij, sji, cji;
                                double ai, bi, aj, bj;
                                double xip_re, xip_im, xim_re, xim_im;
                                
                                const double wn = ww/W_[n];
                                
                                // e^{I phi_ij} unnormalised
                                cij = cyi*syj - syi*cyj*cdx;
                                sij = -cyj*sdx;
                                // cij + I sij = e^{I Si phi_ij}
                                nsincos(Si, cij, sij, &sij, &cij);
                                // ai + I bi = (ui + I vi) e^{-I Si phi_i}
                                cmul(ui, vi, cij, -sij, &ai, &bi);
                                
                                // e^{I phi_ji} unnormalised
                                cji = cyj*syi - syj*cyi*cdx;
                                sji = cyi*sdx;
                                // cji + I cji = e^{I Sj phi_ji}
                                nsincos(Sj, cji, sji, &sji, &cji);
                                // aj + I bj = (uj + I vj) e^{-I Sj phi_ji}
                                cmul(uj, vj, cji, -sji, &aj, &bj);
                                
                                // xip = (ai + I bi) (aj - I bj)
                                cmul(ai, bi, aj, -bj, &xip_re, &xip_im);
                                // xim = (ai + I bi) (aj + I bj)
                                cmul(ai, bi, aj, bj, &xim_re, &xim_im);
                                
                                X_[0*nd+n] += wn*(xip_re - X_[0*nd+n]);
                                X_[1*nd+n] += wn*(xim_re - X_[1*nd+n]);
                                X_[2*nd+n] += wn*(xip_im - X_[2*nd+n]);
                                X_[3*nd+n] += wn*(xim_im - X_[3*nd+n]);
                            }
                        }
                    }
                }
            }
            
            #pragma omp master
            {
                signal(SIGALRM, SIG_IGN);
            }
            
            #pragma omp critical
            {
                for(i = 0; i < nd; ++i)
                {
                    N[p*nd+i] += N_[i];
                    W[p*nd+i] += W_[i];
                    
                    if(X)
                    {
                        const double w = W_[i] > 0 ? W_[i]/W[p*nd+i] : 0;
                        
                        X[0*nd+i] += w*(X_[0*nd+i] - X[0*nd+i]);
                        X[1*nd+i] += w*(X_[1*nd+i] - X[1*nd+i]);
                        X[2*nd+i] += w*(X_[2*nd+i] - X[2*nd+i]);
                        X[3*nd+i] += w*(X_[3*nd+i] - X[3*nd+i]);
                    }
                }
            }
            
            free(qr);
            free(N_);
            free(W_);
            free(X_);
            
            if(tc)
            {
                free(db_);
                free(ci_);
                if(cj != ci)
                    free(cj_);
                free(mj_);
            }
        }
        
        nn = 0;
        for(i = 0; i < nd; ++i)
            nn += N[p*nd+i];
        
        dt = difftime(time(NULL), st);
        
        printf("\r> done with %.0f pairs", nn);
        printf(" in %02d:%02d:%02d  \n", dt/3600, (dt/60)%60, dt%60);
        printf("\n");
    }
    
    if(pt)
    {
        const size_t ndd = n1*(n1-1)/2;
        const size_t ndr = n1*n2;
        const size_t nrr = n2*(n2-1)/2;
        
        for(i = 0; i < nd; ++i)
        {
            W[0*nd+i] /= ndd;
            W[1*nd+i] /= ndr;
            W[2*nd+i] /= nrr;
        }
    }
    
    writexi(cfg.output, nd, dl, dh, ls, uo, pt, N, W, X);
    
    free(W);
    free(X);
    free(m1);
    free(m2);
    free(c1);
    free(c2);
    free(dx);
    
    freecfg(&cfg);
    
    return EXIT_SUCCESS;
    
err_usage:
    fprintf(stderr, "usage: 2pcf [FILE]\n");
    return EXIT_FAILURE;
    
err_alloc:
    perror(NULL);
    return EXIT_FAILURE;
}

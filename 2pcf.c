// 2pcf: two-point correlation function semi-tree code
// ---------------------------------------------------
// author:  Nicolas Tessore <nicolas.tessore@manchester.ac.uk>
// date:    08 Aug 2018

#define _XOPEN_SOURCE 600

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef LINELEN
#define LINELEN 1024
#endif

enum {
    CORR_DD,
    CORR_DR,
    CORR_RR,
    /*-------*/
    CORR_DONE
};

enum {
    UNIT_RAD,
    UNIT_DEG,
    UNIT_ARCMIN,
    UNIT_ARCSEC,
    /*-------*/
    NUM_UNITS
};

enum {
    SPACING_LIN,
    SPACING_LOG
};

enum {
    COORDS_FLAT,
    COORDS_RADEC
};

const double UCONV[NUM_UNITS] = {
    1.0,
    0.017453292519943295769,
    0.00029088820866572159615,
    0.0000048481368110953599359
};

const char* UNAME[NUM_UNITS] = {
    "rad",
    "deg",
    "arcmin",
    "arcsec"
};

static const double TWO_PI = 6.2831853071795864769;

int cmp0(const void* a, const void* b)
{
    const double* x = a;
    const double* y = b;
    return (x[0] > y[0]) - (x[0] < y[0]);
}

int cmp1(const void* a, const void* b)
{
    const double* x = a;
    const double* y = b;
    return (x[1] > y[1]) - (x[1] < y[1]);
}

typedef struct node {
    size_t i;
    size_t n;
    double x[2];
    double y[2];
    struct node* l;
    struct node* r;
} node;

node* tree(double* p, size_t i, size_t n, size_t m, int d, size_t s, size_t* c)
{
    node* t;
    
    t = malloc(sizeof(node));
    if(!t)
    {
        perror("tree()");
        exit(EXIT_FAILURE);
    }
    
    t->i = i;
    t->n = n;
    
    if(d == 0)
    {
        if(c)
            *c = 0;
        qsort(p+i*m, n-i, m*sizeof(double), cmp1);
    }
    
    if(d & 1)
    {
        t->x[0] = p[i*m];
        t->x[1] = p[(n-1)*m];
    }
    else
    {
        t->y[0] = p[i*m+1];
        t->y[1] = p[(n-1)*m+1];
    }
    
    qsort(p+i*m, n-i, m*sizeof(double), d%2 ? cmp1 : cmp0);
    
    if(d & 1)
    {
        t->y[0] = p[i*m+1];
        t->y[1] = p[(n-1)*m+1];
    }
    else
    {
        t->x[0] = p[i*m];
        t->x[1] = p[(n-1)*m];
    }
    
    if(n-i > s)
    {
        t->l = tree(p, i, (i+n)/2, m, d+1, s, c);
        t->r = tree(p, (i+n)/2, n, m, d+1, s, c);
    }
    else
    {
        t->l = NULL;
        t->r = NULL;
    }
    
    if(c)
        *c += 1;
    
    return t;
}

void tree_free(node* t)
{
    if(t)
    {
        tree_free(t->l);
        tree_free(t->r);
        free(t);
    }
}

static double* readc(const char* f, size_t* n, double u, bool rd)
{
    FILE* fp;
    int l;
    char buf[LINELEN];
    size_t i, a;
    double* d;
    char* sx;
    char* sy;
    char* sw;
    double x, y, w;
    
    fp = fopen(f, "r");
    if(!fp)
    {
        perror(f);
        exit(EXIT_FAILURE);
    }
    
    i = 0;
    a = 1;
    d = malloc(a*7*sizeof(double));
    if(!d)
    {
        perror(NULL);
        abort();
    }
    
    for(l = 1; fgets(buf, sizeof buf, fp); ++l)
    {
        sx = strtok(buf, " \t\r\n");
        sy = strtok(NULL, " \t\r\n");
        sw = strtok(NULL, " \t\r\n");
        
        if(!sx || *sx == '#')
            continue;
        
        if(!sy || *sy == '#')
        {
            fprintf(stderr, "error: %s:%d: missing `y` value\n", f, l);
            exit(EXIT_FAILURE);
        }
        
        x = atof(sx)*u;
        y = atof(sy)*u;
        
        if(!sw || *sw == '#')
            w = 1;
        else
            w = atof(sw);
        
        d[i*7+0] = x;
        d[i*7+1] = y;
        d[i*7+2] = w;
        
        if(rd)
        {
            d[i*7+3] = sin(x);
            d[i*7+4] = cos(x);
            d[i*7+5] = sin(y);
            d[i*7+6] = cos(y);
        }
        
        i += 1;
        
        if(i == a)
        {
            a *= 2;
            d = realloc(d, a*7*sizeof(double));
            if(!d)
            {
                perror(NULL);
                abort();
            }
        }
    }
    
    fclose(fp);
    
    d = realloc(d, i*7*sizeof(double));
    if(!d)
    {
        perror(NULL);
        abort();
    }
    
    *n = i;
    
    return d;
}

volatile sig_atomic_t flag_fb = 0;

#pragma omp threadprivate(flag_fb)

void fbhandler(int s)
{
    flag_fb = 1;
}

int main(int argc, char* argv[])
{
    FILE* fp;
    
    const char* infile;
    char buf[LINELEN];
    size_t line;
    char* key;
    char* val;
    
    struct {
        char output[LINELEN];
        char data[LINELEN];
        char random[LINELEN];
        int dunit;
        int coords;
        int nth;
        double thmin;
        double thmax;
        int thunit;
        int spacing;
        int leafpts;
    } cfg;
    
    size_t ndat, nran;
    double* cdat;
    double* cran;
    
    size_t nnod;
    node* tdat;
    node* tran;
    
    bool ls, rd;
    int nd;
    double dl, dh, sdh, dm, d0;
    double ui, uo;
    
    double* W;
    
    int act;
    size_t ni, nj, ii;
    double* ci;
    double* cj;
    node* tj;
    bool xc;
    
    size_t NDD, NDR, NRR;
    double DD, DR, RR;
    
    if(argc > 2)
        goto err_usage;
    
    memset(&cfg, 0, sizeof cfg);
    
    cfg.leafpts = 8;
    
    infile = argc > 1 ? argv[1] : "2pcf.cfg";
    
    fp = fopen(infile, "r");
    if(!fp)
    {
        perror(infile);
        return EXIT_FAILURE;
    }
    
    for(line = 1; fgets(buf, sizeof(buf), fp); ++line)
    {
        key = strtok(buf, " \t\r\n");
        val = strtok(NULL, " \t\r\n");
        
        if(!key || *key == '#')
            continue;
        if(!val || *val == '#')
            goto err_cfg_no_value;
        
        if(strcmp(key, "output") == 0)
            strncpy(cfg.output, val, sizeof cfg.output);
        else if(strcmp(key, "data") == 0)
            strncpy(cfg.data, val, sizeof cfg.data);
        else if(strcmp(key, "random") == 0)
            strncpy(cfg.random, val, sizeof cfg.random);
        else if(strcmp(key, "dunit") == 0)
        {
            for(cfg.dunit = 0; cfg.dunit < NUM_UNITS; ++cfg.dunit)
                if(strcmp(val, UNAME[cfg.dunit]) == 0)
                    break;
            if(cfg.dunit == NUM_UNITS)
                goto err_cfg_bad_value;
        }
        else if(strcmp(key, "coords") == 0)
        {
            if(strcmp(val, "flat") == 0)
                cfg.coords = COORDS_FLAT;
            else if(strcmp(val, "radec") == 0)
                cfg.coords = COORDS_RADEC;
            else
                goto err_cfg_bad_value;
        }
        else if(strcmp(key, "nth") == 0)
            cfg.nth = atoi(val);
        else if(strcmp(key, "thmin") == 0)
            cfg.thmin = atof(val);
        else if(strcmp(key, "thmax") == 0)
            cfg.thmax = atof(val);
        else if(strcmp(key, "thunit") == 0)
        {
            for(cfg.thunit = 0; cfg.thunit < NUM_UNITS; ++cfg.thunit)
                if(strcmp(val, UNAME[cfg.thunit]) == 0)
                    break;
            if(cfg.thunit == NUM_UNITS)
                goto err_cfg_bad_value;
        }
        else if(strcmp(key, "spacing") == 0)
        {
            if(strcmp(val, "lin") == 0)
                cfg.spacing = SPACING_LIN;
            else if(strcmp(val, "log") == 0)
                cfg.spacing = SPACING_LOG;
            else
                goto err_cfg_bad_value;
        }
        else if(strcmp(key, "leafpts") == 0)
            cfg.leafpts = atoi(val);
        else
            goto err_cfg_bad_key;
    }
    
    fclose(fp);
    
    if(!strlen(cfg.output))
        { key = "output"; goto err_cfg_missing_key; }
    if(!strlen(cfg.data))
        { key = "data"; goto err_cfg_missing_key; }
    if(!strlen(cfg.random))
        { key = "random"; goto err_cfg_missing_key; }
    if(!cfg.nth)
        { key = "nth"; goto err_cfg_missing_key; }
    if(!cfg.thmin)
        { key = "thmin"; goto err_cfg_missing_key; }
    if(!cfg.thmax)
        { key = "thmax"; goto err_cfg_missing_key; }
    
    ui = UCONV[cfg.dunit];
    uo = UCONV[cfg.thunit];
    
    rd = cfg.coords == COORDS_RADEC;
    ls = cfg.spacing == SPACING_LOG;
    
    nd = cfg.nth;
    dl = cfg.thmin*uo;
    dh = cfg.thmax*uo;
    
    sdh = sin(dh);
    
    if(ls)
    {
        d0 = log(dl);
        dm = nd/(log(dh) - d0);
    }
    else
    {
        d0 = dl;
        dm = nd/(dh - d0);
    }
    
    printf("\n");
    printf("configuration ... %s\n", infile);
    printf("data catalog .... %s\n", cfg.data);
    printf("random catalog .. %s\n", cfg.random);
    printf("data units ...... %s\n", UNAME[cfg.dunit]);
    printf("coordinates ..... %s\n", rd ? "spherical" : "flat");
    printf("\n");
    printf("output file ..... %s\n", cfg.output);
    printf("bin count ....... %u\n", cfg.nth);
    printf("bin range ....... %lg to %lg %s\n",
                                    cfg.thmin, cfg.thmax, UNAME[cfg.thunit]);
    printf("bin spacing ..... %s\n", ls ? "logarithmic" : "linear");
    printf("\n");
    printf("leaf points ..... %d\n", cfg.leafpts);
    printf("\n");
    
    printf("reading data catalog\n");
    
    cdat = readc(cfg.data, &ndat, ui, rd);
    if(!cdat)
    {
        fprintf(stderr, "error: %s: could not read catalog\n", cfg.data);
        return EXIT_FAILURE;
    }
    
    printf("> done with %zu points\n", ndat);
    printf("\n");
    
    printf("building data tree\n");
    
    tdat = tree(cdat, 0, ndat, 7, 0, cfg.leafpts, &nnod);
    
    printf("> done with %zu nodes\n", nnod);
    printf("\n");
    
    printf("reading random catalog\n");
    
    cran = readc(cfg.random, &nran, ui, rd);
    if(!cran)
    {
        fprintf(stderr, "error: %s: could not read catalog\n", cfg.random);
        return EXIT_FAILURE;
    }
    
    printf("> done with %zu points\n", nran);
    printf("\n");
    
    printf("building random tree\n");
    
    tran = tree(cran, 0, nran, 7, 0, cfg.leafpts, &nnod);
    
    printf("> done with %zu nodes\n", nnod);
    printf("\n");
    
    W = calloc(3*nd, sizeof(double));
    if(!W)
    {
        perror(NULL);
        abort();
    }
    
    for(act = CORR_DD; act != CORR_DONE; ++act)
    {
        printf("calculating %s correlations\n",
                    act == CORR_DD ? "DD" : act == CORR_DR ? "DR" : "RR");
        
        switch(act)
        {
        case CORR_DD:
            ci = cdat;
            ni = ndat;
            
            cj = cdat;
            tj = tdat;
            
            xc = false;
            
            break;
            
        case CORR_DR:
            ci = cdat;
            ni = ndat;
            
            cj = cran;
            tj = tran;
            
            xc = true;
            
            break;
            
        default:
            ci = cran;
            ni = nran;
            
            cj = cran;
            tj = tran;
            
            xc = false;
            
            break;
        }
        
        ii = 0;
        
        #pragma omp parallel
        {
            time_t Tini, Tnow;
            int dT;
            
            node** tl;
            size_t ta;
            
            size_t i, j, k;
            double xi, yi, wi, xj, yj, wj, xk, yk, dx, xl, xh;
            double sxi, cxi, syi, cyi, sxj, cxj, syj, cyj, d;
            int n;
            
            double* Wi;
            
            ta = 1;
            tl = malloc(ta*sizeof(node*));
            if(!tl)
            {
                perror(NULL);
                abort();
            }
            
            Wi = calloc(3*nd, sizeof(double));
            if(!Wi)
            {
                perror(NULL);
                abort();
            }
            
            #pragma omp master
            {
                time(&Tini);
                signal(SIGALRM, fbhandler);
                alarm(1);
            }
            
            sxi = cxi = syi = cyi = 0;
            
            #pragma omp for schedule(static, 1) nowait
            for(i = 0; i < ni; ++i)
            {
                #pragma omp atomic
                ii += 1;
                
                if(flag_fb)
                {
                    time(&Tnow);
                    dT = difftime(Tnow, Tini);
                    printf("\r> %.2f%% - %02d:%02d:%02d ",
                                100.*ii/ni, dT/3600, (dT/60)%60, dT%60);
                    fflush(stdout);
                    flag_fb = 0;
                    alarm(1);
                }
                
                xi = ci[i*7+0];
                yi = ci[i*7+1];
                wi = ci[i*7+2];
                
                if(rd)
                {
                    sxi = ci[i*7+3];
                    cxi = ci[i*7+4];
                    syi = ci[i*7+5];
                    cyi = ci[i*7+6];
                }
                
                tl[0] = tj;
                
                for(k = 1; k > 0; --k)
                {
                    if(!xc && tl[k-1]->n <= i+1)
                        continue;
                    
                    xj = tl[k-1]->x[0];
                    yj = tl[k-1]->y[0];
                    xk = tl[k-1]->x[1];
                    yk = tl[k-1]->y[1];
                    
                    if(yj - yi >= dh || yi - yk >= dh)
                        continue;
                    
                    if(rd)
                    {
                        if(cyi > sdh)
                        {
                            dx = asin(sdh/cyi);
                            xl = xi - dx;
                            xh = xi + dx;
                            if(xl < 0)
                            {
                                if(xl + TWO_PI >= xk && xh <= xj)
                                    continue;
                            }
                            else if(xh > TWO_PI)
                            {
                                if(xl >= xk && xh - TWO_PI <= xj)
                                    continue;
                            }
                            else
                            {
                                if(xl >= xk || xh <= xj)
                                    continue;
                            }
                        }
                    }
                    else
                    {
                        if(xj - xi >= dh || xi - xk >= dh)
                            continue;
                    }
                    
                    if(tl[k-1]->l)
                    {
                        if(k == ta)
                        {
                            ta *= 2;
                            tl = realloc(tl, ta*sizeof(node*));
                            if(!tl)
                            {
                                perror(NULL);
                                abort();
                            }
                        }
                        
                        tl[k] = tl[k-1]->l;
                        tl[k-1] = tl[k-1]->r;
                        
                        k += 2;
                        
                        continue;
                    }
                    
                    j = tl[k-1]->i;
                    nj = tl[k-1]->n;
                    
                    if(!xc && j < i+1)
                        j = i+1;
                    
                    for(; j < nj; ++j)
                    {
                        xj = cj[j*7+0];
                        yj = cj[j*7+1];
                        wj = cj[j*7+2];
                        
                        if(rd)
                        {
                            sxj = cj[j*7+3];
                            cxj = cj[j*7+4];
                            syj = cj[j*7+5];
                            cyj = cj[j*7+6];
                            
                            d = acos(syi*syj + cyi*cyj*(sxi*sxj + cxi*cxj));
                        }
                        else
                            d = hypot(xj - xi, yj - yi);
                        
                        if(d < dl || d >= dh)
                            continue;
                        
                        if(ls)
                            d = log(d);
                        
                        n = dm*(d - d0);
                        
                        Wi[act*nd+n] += wi*wj;
                    }
                }
            }
            
            #pragma omp master
            {
                signal(SIGALRM, SIG_IGN);
                
                time(&Tnow);
                dT = difftime(Tnow, Tini);
                printf("\r> done in %02d:%02d:%02d  \n",
                                            dT/3600, (dT/60)%60, dT%60);
                printf("\n");
            }
            
            #pragma omp critical
            for(n = 0; n < 3*nd; ++n)
                W[n] += Wi[n];
            
            free(tl);
            free(Wi);
        }
    }
    
    NDD = ndat*(ndat-1)/2;
    NDR = ndat*nran;
    NRR = nran*(nran-1)/2;
    
    fp = fopen(cfg.output, "w");
    if(!fp)
    {
        perror(cfg.output);
        return EXIT_FAILURE;
    }
    
    fprintf(fp, "# theta w\n");
    
    for(int n = 0; n < nd; ++n)
    {
        double d = d0 + (n + 0.5)/dm;
        
        if(ls)
            d = exp(d);
        
        DD = W[0*nd+n]/NDD;
        DR = W[1*nd+n]/NDR;
        RR = W[2*nd+n]/NRR;
        
        fprintf(fp, "%.18e %.18e\n", d/uo, (DD - 2*DR + RR)/RR);
    }
    
    fclose(fp);
    
    free(W);
    tree_free(tdat);
    tree_free(tran);
    free(cdat);
    free(cran);
    
    return EXIT_SUCCESS;
    
err_usage:
    fprintf(stderr, "usage: 2pcf [FILE]\n");
    return EXIT_FAILURE;
    
err_cfg_bad_key:
    fprintf(stderr, "error: %s:%zu: invalid key `%s`\n", infile, line, key);
    return EXIT_FAILURE;
    
err_cfg_missing_key:
    fprintf(stderr, "error: %s: missing required `%s` key\n", infile, key);
    return EXIT_FAILURE;
    
err_cfg_no_value:
    fprintf(stderr, "error: %s:%zu: missing value\n", infile, line);
    return EXIT_FAILURE;
    
err_cfg_bad_value:
    fprintf(stderr, "error: %s:%zu: invalid value `%s`\n", infile, line, val);
    return EXIT_FAILURE;
}

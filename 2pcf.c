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
        d[i*7+3] = sin(x);
        d[i*7+4] = cos(x);
        d[i*7+5] = sin(y);
        d[i*7+6] = cos(y);
        
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
        char catalog1[LINELEN];
        char catalog2[LINELEN];
        int dunit;
        int coords;
        int nth;
        double thmin;
        double thmax;
        int thunit;
        int spacing;
        int leafpts;
    } cfg;
    
    bool xc, ls, rd;
    int nd;
    double dl, dh, sdh, dm, d0;
    double ui, uo;
    
    size_t nc1, nc2;
    double* c1;
    double* c2;
    
    size_t nn;
    node* t1;
    node* t2;
    
    double* W;
    
    int p, np;
    size_t ni, nj, ii;
    double* ci;
    double* cj;
    node* tj;
    
    size_t NDD, NDR, NRR;
    double DD, DR, RR;
    
    if(argc > 2)
        goto err_usage;
    
    infile = argc > 1 ? argv[1] : "2pcf.cfg";
    
    fp = fopen(infile, "r");
    if(!fp)
    {
        perror(infile);
        return EXIT_FAILURE;
    }
    
    memset(&cfg, 0, sizeof cfg);
    
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
        else if(strcmp(key, "catalog") == 0)
            strncpy(cfg.catalog1, val, sizeof cfg.catalog1);
        else if(strcmp(key, "catalog1") == 0)
            strncpy(cfg.catalog1, val, sizeof cfg.catalog1);
        else if(strcmp(key, "catalog2") == 0)
            strncpy(cfg.catalog2, val, sizeof cfg.catalog2);
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
    if(!strlen(cfg.catalog1))
        { key = "catalog"; goto err_cfg_missing_key; }
    if(!cfg.nth)
        { key = "nth"; goto err_cfg_missing_key; }
    if(!cfg.thmin)
        { key = "thmin"; goto err_cfg_missing_key; }
    if(!cfg.thmax)
        { key = "thmax"; goto err_cfg_missing_key; }
    
    if(!cfg.leafpts)
        cfg.leafpts = 8;
    
    xc = !!strlen(cfg.catalog2);
    rd = cfg.coords == COORDS_RADEC;
    ls = cfg.spacing == SPACING_LOG;
    
    if(!xc)
    {
        fprintf(stderr, "error: requires `catalog2` for random points\n");
        return EXIT_FAILURE;
    }
    
    ui = UCONV[cfg.dunit];
    uo = UCONV[cfg.thunit];
    
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
    if(!xc)
        printf("catalog ......... %s\n", cfg.catalog1);
    else
    {
        printf("catalog 1 ....... %s\n", cfg.catalog1);
        printf("catalog 2 ....... %s\n", cfg.catalog2);
    }
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
    
    printf("reading catalog%s\n", xc ? " 1" : "");
    
    c1 = readc(cfg.catalog1, &nc1, ui, rd);
    if(!c1)
    {
        fprintf(stderr, "error: could not read %s\n", cfg.catalog1);
        return EXIT_FAILURE;
    }
    
    printf("> done with %zu points\n", nc1);
    printf("\n");
    
    printf("building tree%s\n", xc ? " 1" : "");
    
    t1 = tree(c1, 0, nc1, 7, 0, cfg.leafpts, &nn);
    
    printf("> done with %zu nodes\n", nn);
    printf("\n");
    
    if(xc)
    {
        printf("reading catalog 2\n");
        
        c2 = readc(cfg.catalog2, &nc2, ui, rd);
        if(!c2)
        {
            fprintf(stderr, "error: could not read %s\n", cfg.catalog2);
            return EXIT_FAILURE;
        }
        
        printf("> done with %zu points\n", nc2);
        printf("\n");
        
        printf("building tree 2\n");
        
        t2 = tree(c2, 0, nc2, 7, 0, cfg.leafpts, &nn);
        
        printf("> done with %zu nodes\n", nn);
        printf("\n");
    }
    
    W = calloc(3*nd, sizeof(double));
    if(!W)
    {
        perror(NULL);
        abort();
    }
    
    for(p = 0, np = 3; p < np; ++p)
    {
        switch(p)
        {
        case 0:
            printf("calculating DD correlations\n");
            
            xc = false;
            
            ci = c1;
            ni = nc1;
            
            cj = c1;
            tj = t1;
            
            break;
            
        case 1:
            printf("calculating DR correlations\n");
            
            xc = true;
            
            ci = c1;
            ni = nc1;
            
            cj = c2;
            tj = t2;
            
            break;
            
        case 2:
            printf("calculating RR correlations\n");
            
            xc = false;
            
            ci = c2;
            ni = nc2;
            
            cj = c2;
            tj = t2;
            
            break;
        }
        
        ii = 0;
        
        #pragma omp parallel
        {
            time_t Tini, Tnow;
            int dT;
            
            size_t tn, ta;
            node** tl;
            
            size_t i, j;
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
                
                xi  = ci[i*7+0];
                yi  = ci[i*7+1];
                wi  = ci[i*7+2];
                sxi = ci[i*7+3];
                cxi = ci[i*7+4];
                syi = ci[i*7+5];
                cyi = ci[i*7+6];
                
                tl[0] = tj;
                
                for(tn = 1; tn > 0; --tn)
                {
                    if(!xc && tl[tn-1]->n <= i+1)
                        continue;
                    
                    xj = tl[tn-1]->x[0];
                    yj = tl[tn-1]->y[0];
                    xk = tl[tn-1]->x[1];
                    yk = tl[tn-1]->y[1];
                    
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
                    
                    if(tl[tn-1]->l)
                    {
                        if(tn == ta)
                        {
                            ta *= 2;
                            tl = realloc(tl, ta*sizeof(node*));
                            if(!tl)
                            {
                                perror(NULL);
                                abort();
                            }
                        }
                        
                        tl[tn] = tl[tn-1]->l;
                        tl[tn-1] = tl[tn-1]->r;
                        
                        tn += 2;
                        
                        continue;
                    }
                    
                    j = tl[tn-1]->i;
                    nj = tl[tn-1]->n;
                    
                    if(!xc && j < i+1)
                        j = i+1;
                    
                    for(; j < nj; ++j)
                    {
                        xj  = cj[j*7+0];
                        yj  = cj[j*7+1];
                        wj  = cj[j*7+2];
                        sxj = cj[j*7+3];
                        cxj = cj[j*7+4];
                        syj = cj[j*7+5];
                        cyj = cj[j*7+6];
                        
                        if(rd)
                            d = acos(syi*syj + cyi*cyj*(sxi*sxj + cxi*cxj));
                        else
                            d = hypot(xj - xi, yj - yi);
                        
                        if(d < dl || d >= dh)
                            continue;
                        
                        if(ls)
                            d = log(d);
                        
                        n = dm*(d - d0);
                        
                        Wi[p*nd+n] += wi*wj;
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
    
    NDD = nc1*(nc1-1)/2;
    NDR = nc1*nc2;
    NRR = nc2*(nc2-1)/2;
    
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
    tree_free(t1);
    tree_free(t2);
    free(c1);
    free(c2);
    
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

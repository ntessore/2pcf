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

enum { MODE_NONE, MODE_POINTS, MODE_FIELD };
enum { SPACING_LIN, SPACING_LOG };
enum { COORDS_FLAT, COORDS_RADEC };
enum { FIELD_REAL, FIELD_COMPLEX };

enum {
    UNIT_RAD,
    UNIT_DEG,
    UNIT_ARCMIN,
    UNIT_ARCSEC,
    /*-------*/
    NUM_UNITS
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

static const int DW = 9;

static double clamp(double x, double xl, double xh)
{
  const double y = x < xl ? xl : x;
  return y > xh ? xh : y;
}

static void nsincos(int n, double x, double y, double* s, double* c)
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
        {
            const double k1 = x*(*c + *s);
            const double k2 = (x + y)*(*s);
            const double k3 = (y - x)*(*c);
            *c = k1 - k2;
            *s = k1 + k3;
        }
        return;
    }
}

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
    size_t i, n;
    double xl, xh;
    double yl, yh;
    struct node* l;
    struct node* r;
} node;

node* tree(double* p, size_t i, size_t n, int d, size_t s, size_t* c)
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
        qsort(p+i*DW, n-i, DW*sizeof(double), cmp1);
    }
    
    if(d & 1)
    {
        t->xl = p[i*DW+0];
        t->xh = p[(n-1)*DW+0];
    }
    else
    {
        t->yl = p[i*DW+1];
        t->yh = p[(n-1)*DW+1];
    }
    
    qsort(p+i*DW, n-i, DW*sizeof(double), d & 1 ? cmp1 : cmp0);
    
    if(d & 1)
    {
        t->yl = p[i*DW+1];
        t->yh = p[(n-1)*DW+1];
    }
    else
    {
        t->xl = p[i*DW+0];
        t->xh = p[(n-1)*DW+0];
    }
    
    if(n-i > s)
    {
        t->l = tree(p, i, (i+n)/2, d+1, s, c);
        t->r = tree(p, (i+n)/2, n, d+1, s, c);
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

double* readc(const char* f, double ui, bool rd, bool pt, bool cf, size_t* n)
{
    FILE* fp;
    int l;
    char buf[LINELEN];
    size_t i, a;
    double* d;
    char* sx;
    char* sy;
    char* su;
    char* sv;
    char* sw;
    double x, y, u, v, w;
    
    fp = fopen(f, "r");
    if(!fp)
    {
        perror(f);
        exit(EXIT_FAILURE);
    }
    
    i = 0;
    a = 1;
    
    d = malloc(a*DW*sizeof(double));
    if(!d)
    {
        perror(NULL);
        abort();
    }
    
    for(l = 1; fgets(buf, sizeof buf, fp); ++l)
    {
        sx = strtok(buf, " \t\r\n");
        sy = strtok(NULL, " \t\r\n");
        su = !pt ? strtok(NULL, " \t\r\n") : NULL;
        sv = !pt && cf ? strtok(NULL, " \t\r\n") : NULL;
        sw = strtok(NULL, " \t\r\n");
        
        if(!sx || *sx == '#')
            continue;
        
        if(!sy || *sy == '#')
        {
            fprintf(stderr, "error: %s:%d: missing `y` value\n", f, l);
            exit(EXIT_FAILURE);
        }
        
        x = atof(sx)*ui;
        y = atof(sy)*ui;
        
        if(!pt)
        {
            if(!su || *su == '#')
            {
                fprintf(stderr, "error: %s:%d: missing `u` value\n", f, l);
                exit(EXIT_FAILURE);
            }
            
            u = atof(su);
            
            if(cf)
            {
                if(!sv || *sv == '#')
                {
                    fprintf(stderr, "error: %s:%d: missing `v` value\n", f, l);
                    exit(EXIT_FAILURE);
                }
                
                v = atof(sv);
            }
            else
                v = 0;
        }
        else
        {
            u = 0;
            v = 0;
        }
        
        if(!sw || *sw == '#')
            w = 1;
        else
            w = atof(sw);
        
        d[i*DW+0] = x;
        d[i*DW+1] = y;
        d[i*DW+2] = u;
        d[i*DW+3] = v;
        d[i*DW+4] = w;
        d[i*DW+5] = rd ? sin(x) : 0;
        d[i*DW+6] = rd ? cos(x) : 0;
        d[i*DW+7] = rd ? sin(y) : 0;
        d[i*DW+8] = rd ? cos(y) : 0;
        
        i += 1;
        
        if(i == a)
        {
            a *= 2;
            d = realloc(d, a*DW*sizeof(double));
            if(!d)
            {
                perror(NULL);
                abort();
            }
        }
    }
    
    fclose(fp);
    
    d = realloc(d, i*DW*sizeof(double));
    if(!d)
    {
        perror(NULL);
        abort();
    }
    
    *n = i;
    
    return d;
}

volatile sig_atomic_t fb;

#pragma omp threadprivate(fb)

void fbhandler(int s)
{
    fb = s;
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
        int mode;
        char catalog1[LINELEN];
        char catalog2[LINELEN];
        int dunit;
        int coords;
        int field1;
        int field2;
        int spin1;
        int spin2;
        char output[LINELEN];
        int nth;
        double thmin;
        double thmax;
        int thunit;
        int spacing;
        int leafpts;
    } cfg;
    
    bool pt, xc, ls, rd;
    size_t nd;
    double dl, dh, sdh, dm, d0;
    double ui, uo;
    int S1, S2;
    
    size_t n1, n2;
    double* c1;
    double* c2;
    
    size_t nn;
    node* t1;
    node* t2;
    
    double* W;
    double* X;
    
    int p, np;
    size_t ni, nj, ii;
    double* ci;
    double* cj;
    node* tj;
    int Si, Sj;
    
    size_t i, j;
    
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
        
        if(strcmp(key, "mode") == 0)
        {
            if(strcmp(val, "points") == 0)
                cfg.mode = MODE_POINTS;
            else if(strcmp(val, "field") == 0)
                cfg.mode = MODE_FIELD;
            else
                goto err_cfg_bad_value;
        }
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
        else if(strcmp(key, "field") == 0)
        {
            if(strcmp(val, "real") == 0)
                cfg.field1 = FIELD_REAL;
            else if(strcmp(val, "complex") == 0)
                cfg.field1 = FIELD_COMPLEX;
            else
                goto err_cfg_bad_value;
        }
        else if(strcmp(key, "field1") == 0)
        {
            if(strcmp(val, "real") == 0)
                cfg.field1 = FIELD_REAL;
            else if(strcmp(val, "complex") == 0)
                cfg.field1 = FIELD_COMPLEX;
            else
                goto err_cfg_bad_value;
        }
        else if(strcmp(key, "field2") == 0)
        {
            if(strcmp(val, "real") == 0)
                cfg.field2 = FIELD_REAL;
            else if(strcmp(val, "complex") == 0)
                cfg.field2 = FIELD_COMPLEX;
            else
                goto err_cfg_bad_value;
        }
        else if(strcmp(key, "spin") == 0)
        {
            cfg.spin1 = atoi(val) + 1;
            if(cfg.spin1 < 1)
                goto err_cfg_bad_value;
        }
        else if(strcmp(key, "spin1") == 0)
        {
            cfg.spin1 = atoi(val) + 1;
            if(cfg.spin1 < 1)
                goto err_cfg_bad_value;
        }
        else if(strcmp(key, "spin2") == 0)
        {
            cfg.spin2 = atoi(val) + 1;
            if(cfg.spin2 < 1)
                goto err_cfg_bad_value;
        }
        else if(strcmp(key, "output") == 0)
            strncpy(cfg.output, val, sizeof cfg.output);
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
    
    if(cfg.mode == MODE_NONE)
        { key = "mode"; goto err_cfg_missing_key; }
    if(!strlen(cfg.catalog1))
        { key = "catalog"; goto err_cfg_missing_key; }
    if(!strlen(cfg.output))
        { key = "output"; goto err_cfg_missing_key; }
    if(!cfg.nth)
        { key = "nth"; goto err_cfg_missing_key; }
    if(!cfg.thmin)
        { key = "thmin"; goto err_cfg_missing_key; }
    if(!cfg.thmax)
        { key = "thmax"; goto err_cfg_missing_key; }
    
    if(!cfg.leafpts)
        cfg.leafpts = 8;
    
    if(cfg.spin1)
        cfg.field1 = FIELD_COMPLEX;
    if(cfg.spin2 || cfg.spin1)
        cfg.field2 = FIELD_COMPLEX;
    
    pt = cfg.mode == MODE_POINTS;
    xc = !!strlen(cfg.catalog2);
    rd = cfg.coords == COORDS_RADEC;
    ls = cfg.spacing == SPACING_LOG;
    
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
    
    S1 = cfg.spin1 ? cfg.spin1 - 1 : 0;
    S2 = cfg.spin2 ? cfg.spin2 - 1 : S1;
    
    printf("\n");
    printf("configuration ... %s\n", infile);
    printf("\n");
    printf("input type ...... %s\n", pt ? "points" : "field");
    if(pt)
    {
        printf("data catalog .... %s\n", cfg.catalog1);
        printf("random catalog .. %s\n", cfg.catalog2);
    }
    else if(!xc)
        printf("catalog ......... %s\n", cfg.catalog1);
    else
    {
        printf("catalog 1 ....... %s\n", cfg.catalog1);
        printf("catalog 2 ....... %s\n", cfg.catalog2);
    }
    printf("data units ...... %s\n", UNAME[cfg.dunit]);
    printf("coordinates ..... %s\n", rd ? "spherical" : "flat");
    if(!pt)
    {
        if(!xc)
        {
            if(cfg.field1)
                printf("field type ...... complex spin-%d\n", S1);
            else
                printf("field type ...... real\n");
        }
        else
        {
            if(cfg.field1)
                printf("field 1 type .... complex spin-%d\n", S1);
            else
                printf("field 1 type .... real\n");
            if(cfg.field2)
                printf("field 2 type .... complex spin-%d\n", S2);
            else
                printf("field 2 type .... real\n");
        }
    }
    printf("\n");
    printf("output file ..... %s\n", cfg.output);
    printf("bin count ....... %u\n", cfg.nth);
    printf("bin range ....... %lg to %lg %s\n",
                                    cfg.thmin, cfg.thmax, UNAME[cfg.thunit]);
    printf("bin spacing ..... %s\n", ls ? "logarithmic" : "linear");
    printf("\n");
    printf("leaf points ..... %d\n", cfg.leafpts);
    printf("\n");
    
    printf("reading %s\n", pt ? "data catalog" : xc ? "catalog 1" : "catalog");
    
    c1 = readc(cfg.catalog1, ui, rd, pt, cfg.field1, &n1);
    if(!c1)
    {
        fprintf(stderr, "error: could not read %s\n", cfg.catalog1);
        return EXIT_FAILURE;
    }
    
    printf("> done with %zu points\n", n1);
    printf("\n");
    
    printf("building %s\n", pt ? "data tree" : xc ? "tree 1" : "tree");
    
    t1 = tree(c1, 0, n1, 0, cfg.leafpts, &nn);
    
    printf("> done with %zu nodes\n", nn);
    printf("\n");
    
    if(xc)
    {
        printf("reading %s\n", pt ? "random catalog" : "catalog 2");
        
        c2 = readc(cfg.catalog2, ui, rd, pt, cfg.field2, &n2);
        if(!c2)
        {
            fprintf(stderr, "error: could not read %s\n", cfg.catalog2);
            return EXIT_FAILURE;
        }
        
        printf("> done with %zu points\n", n2);
        printf("\n");
        
        printf("building %s\n", pt ? "random tree" : "tree 2");
        
        t2 = tree(c2, 0, n2, 0, cfg.leafpts, &nn);
        
        printf("> done with %zu nodes\n", nn);
        printf("\n");
    }
    else
    {
        c2 = NULL;
        t2 = NULL;
    }
    
    W = calloc(3*nd, sizeof(double));
    X = calloc(3*nd, sizeof(double));
    if(!W || !X)
    {
        perror(NULL);
        abort();
    }
    
    for(p = 0, np = pt ? 3 : 1; p < np; ++p)
    {
        if(pt)
        {
            Si = Sj = 0;
            
            switch(p)
            {
            case 0:
                printf("calculating DD correlations\n");
                
                xc = false;
                
                ci = c1;
                ni = n1;
                
                cj = c1;
                tj = t1;
                
                break;
                
            case 1:
                printf("calculating DR correlations\n");
                
                xc = true;
                
                ci = c1;
                ni = n1;
                
                cj = c2;
                tj = t2;
                
                break;
                
            case 2:
                printf("calculating RR correlations\n");
                
                xc = false;
                
                ci = c2;
                ni = n2;
                
                cj = c2;
                tj = t2;
                
                break;
            }
        }
        else
        {
            printf("calculating correlations\n");
            
            ci = c1;
            ni = n1;
            Si = S1;
            
            cj = xc ? c2 : c1;
            tj = xc ? t2 : t1;
            Sj = xc ? S2 : S1;
        }
        
        ii = 0;
        nn = 0;
        
        #pragma omp parallel default(none) shared(ii, nn, W, X) \
            private(i, j, nj) firstprivate(pt, xc, rd, ls, nd, dl, dh, d0, \
                dm, sdh, p, ni, ci, cj, tj, Si, Sj, stdout)
        {
            time_t st;
            int dt;
            
            size_t tn, ta;
            node** tl;
            
            double xi, yi, ui, vi, wi, sxi, cxi, syi, cyi;
            double xj, yj, uj, vj, wj, sxj, cxj, syj, cyj, sdx, cdx;
            double xl, xh, yl, yh, dx;
            
            double d;
            int n;
            
            double ww, uu, uv, vu, vv;
            double sij, cij, sji, cji;
            double sp, cp, sm, cm;
            double xip, xim, xix;
            
            double* Wi;
            double* Xi;
            
            ta = 1;
            tl = malloc(ta*sizeof(node*));
            if(!tl)
            {
                perror(NULL);
                abort();
            }
            
            Wi = calloc(3*nd, sizeof(double));
            Xi = calloc(3*nd, sizeof(double));
            if(!Wi || !Xi)
            {
                perror(NULL);
                abort();
            }
            
            fb = 0;
            st = time(NULL);
            dt = 0;
            
            #pragma omp master
            {
                signal(SIGALRM, fbhandler);
                alarm(1);
            }
            
            #pragma omp for schedule(static, 1) nowait
            for(i = 0; i < ni; ++i)
            {
                #pragma omp atomic
                ii += 1;
                
                if(fb)
                {
                    dt = difftime(time(NULL), st);
                    
                    printf("\r> %.2f%%", 100.*ii/ni);
                    printf(" - %02d:%02d:%02d ", dt/3600, (dt/60)%60, dt%60);
                    fflush(stdout);
                    
                    fb = 0;
                    alarm(1);
                }
                
                xi  = ci[i*DW+0];
                yi  = ci[i*DW+1];
                ui  = ci[i*DW+2];
                vi  = ci[i*DW+3];
                wi  = ci[i*DW+4];
                sxi = ci[i*DW+5];
                cxi = ci[i*DW+6];
                syi = ci[i*DW+7];
                cyi = ci[i*DW+8];
                
                tl[0] = tj;
                
                for(tn = 1; tn > 0; --tn)
                {
                    if(!xc && tl[tn-1]->n <= i+1)
                        continue;
                    
                    xl = tl[tn-1]->xl;
                    xh = tl[tn-1]->xh;
                    yl = tl[tn-1]->yl;
                    yh = tl[tn-1]->yh;
                    
                    if(yi - dh >= yh || yi + dh <= yl)
                        continue;
                    
                    if(rd)
                    {
                        if(cyi > sdh)
                        {
                            dx = asin(sdh/cyi);
                            if((xi - dx >= xh && xi + dx - TWO_PI <= xl) ||
                                    (xi + dx <= xl && xi - dx + TWO_PI >= xh))
                                continue;
                        }
                    }
                    else
                    {
                        if(xi - dh >= xh || xi + dh <= xl)
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
                        xj  = cj[j*DW+0];
                        yj  = cj[j*DW+1];
                        uj  = cj[j*DW+2];
                        vj  = cj[j*DW+3];
                        wj  = cj[j*DW+4];
                        sxj = cj[j*DW+5];
                        cxj = cj[j*DW+6];
                        syj = cj[j*DW+7];
                        cyj = cj[j*DW+8];
                        
                        sdx = cxi*sxj - sxi*cxj;
                        cdx = cxi*cxj + sxi*sxj;
                        
                        if(rd)
                            d = acos(clamp(syi*syj + cyi*cyj*cdx, -1, 1));
                        else
                            d = hypot(xj - xi, yj - yi);
                        
                        if(d < dl || d >= dh)
                            continue;
                        
                        if(ls)
                            d = log(d);
                        
                        n = dm*(d - d0);
                        
                        ww = wi*wj;
                        
                        if(!pt)
                        {
                            uu = ui*uj;
                            uv = ui*vj;
                            vu = vi*uj;
                            vv = vi*vj;
                            
                            if(rd)
                            {
                                cij = cyi*syj - syi*cyj*cdx;
                                sij = cyj*sdx;
                                cji = cyj*syi - syj*cyi*cdx;
                                sji = -cyi*sdx;
                            }
                            else
                            {
                                cij = xj - xi;
                                sij = yj - yi;
                                cji = xi - xj;
                                sji = yi - yj;
                            }
                            
                            nsincos(Si, cij, sij, &sij, &cij);
                            nsincos(Sj, cji, sji, &sji, &cji);
                            
                            cp = cij*cji + sij*sji;
                            sp = sij*cji - cij*sji;
                            cm = cij*cji - sij*sji;
                            sm = sij*cji + cij*sji;
                            
                            xip = (uu + vv)*cp - (uv - vu)*sp;
                            xim = (uu - vv)*cm - (uv + vu)*sm;
                            xix = (uv + vu)*cm + (uu - vv)*sm;
                            
                            Xi[0*nd+n] += ww*xip;
                            Xi[1*nd+n] += ww*xim;
                            Xi[2*nd+n] += ww*xix;
                        }
                        
                        Wi[p*nd+n] += ww;
                        
                        #pragma omp atomic
                        nn += 1;
                    }
                }
            }
            
            #pragma omp critical
            for(i = 0; i < 3*nd; ++i)
            {
                W[i] += Wi[i];
                X[i] += Xi[i];
            }
            
            free(tl);
            free(Wi);
            free(Xi);
            
            #pragma omp barrier
            #pragma omp master
            {
                signal(SIGALRM, SIG_IGN);
                alarm(0);
                
                dt = difftime(time(NULL), st);
                
                printf("\r> done with %zu pairs", nn);
                printf(" in %02d:%02d:%02d  \n", dt/3600, (dt/60)%60, dt%60);
                printf("\n");
            }
        }
    }
    
    fp = fopen(cfg.output, "w");
    if(!fp)
    {
        perror(cfg.output);
        return EXIT_FAILURE;
    }
    
    if(pt)
        fprintf(fp, "# theta w\n");
    else
        fprintf(fp, "# theta xip xim xix\n");
    
    for(i = 0; i < nd; ++i)
    {
        double d;
        
        d = d0 + (i + 0.5)/dm;
        if(ls)
            d = exp(d);
        d /= uo;
        
        if(pt)
        {
            size_t ndd, ndr, nrr;
            double dd, dr, rr;
            
            ndd = n1*(n1-1)/2;
            ndr = n1*n2;
            nrr = n2*(n2-1)/2;
            
            dd = W[0*nd+i]/ndd;
            dr = W[1*nd+i]/ndr;
            rr = W[2*nd+i]/nrr;
            
            fprintf(fp, "%.18e %.18e\n", d, (dd - 2*dr + rr)/rr);
        }
        else
        {
            double nor, xip, xim, xix;
            
            nor = W[i];
            xip = X[0*nd+i]/nor;
            xim = X[1*nd+i]/nor;
            xix = X[2*nd+i]/nor;
            
            fprintf(fp, "%.18e %.18e %.18e %.18e\n", d, xip, xim, xix);
        }
    }
    
    fclose(fp);
    
    free(W);
    free(X);
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

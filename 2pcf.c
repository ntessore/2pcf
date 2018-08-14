// 2pcf: two-point correlation function code
// ---
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

static const double PI_HALF = 1.5707963267948966192;
static const double TWO_PI = 6.2831853071795864769;

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

static const int DW = 8;

double* readc(const char* f, double ui, bool pt, bool cf, size_t* n)
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
        d[i*DW+2] = 1;
        d[i*DW+3] = 1;
        d[i*DW+4] = u;
        d[i*DW+5] = v;
        d[i*DW+6] = w;
        
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
        double gridx;
        double gridy;
    } cfg;
    
    bool pt, xc, ls, rd;
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
    
    double* W;
    double* X;
    
    int p, np;
    time_t st;
    int dt;
    size_t ni, nj, ii, nn;
    double* ci;
    double* cj;
    int* mj;
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
        else if(strcmp(key, "gridx") == 0)
            cfg.gridx = atof(val);
        else if(strcmp(key, "gridy") == 0)
            cfg.gridy = atof(val);
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
    
    if(cfg.gridx <= 0)
        cfg.gridx = 0.1;
    if(cfg.gridy <= 0)
        cfg.gridy = 0.1;
    
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
    
    db = malloc((nd+1)*sizeof(double));
    if(!db)
    {
        perror(NULL);
        abort();
    }
    
    for(i = 0; i <= nd; ++i)
    {
        if(ls)
            db[i] = exp(log(dl) + i*(log(dh) - log(dl))/nd);
        else
            db[i] = dl + i*(dh - dl)/nd;
        if(rd)
            db[i] = 2*sin(0.5*db[i]);
        db[i] = db[i]*db[i];
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
    if(rd)
    {
        printf("grid size ....... %g x %g deg^2\n", cfg.gridx, cfg.gridy);
        printf("\n");
    }
    
    printf("reading %s\n", pt ? "data catalog" : xc ? "catalog 1" : "catalog");
    
    c1 = readc(cfg.catalog1, ui, pt, cfg.field1, &n1);
    if(!c1)
    {
        fprintf(stderr, "error: could not read %s\n", cfg.catalog1);
        return EXIT_FAILURE;
    }
    
    printf("> done with %zu points\n", n1);
    printf("\n");
    
    if(xc)
    {
        printf("reading %s\n", pt ? "random catalog" : "catalog 2");
        
        c2 = readc(cfg.catalog2, ui, pt, cfg.field2, &n2);
        if(!c2)
        {
            fprintf(stderr, "error: could not read %s\n", cfg.catalog2);
            return EXIT_FAILURE;
        }
        
        printf("> done with %zu points\n", n2);
        printf("\n");
    }
    else
    {
        c2 = NULL;
        n2 = 0;
    }
    
    printf("building index\n");
    
    if(rd)
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
    {
        perror(NULL);
        abort();
    }
    
    if(rd)
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
    {
        perror(NULL);
        abort();
    }
    
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
        {
            perror(NULL);
            abort();
        }
        
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
    
    if(rd)
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
    
    W = calloc(3*nd, sizeof(double));
    X = calloc(4*nd, sizeof(double));
    if(!W || !X)
    {
        perror(NULL);
        abort();
    }
    
    signal(SIGQUIT, handler);
    qQ = 0;
    
    for(p = 0, np = pt ? 3 : 1; p < np && !qQ; ++p)
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
                mj = m1;
                
                break;
                
            case 1:
                printf("calculating DR correlations\n");
                
                xc = true;
                
                ci = c1;
                ni = n1;
                
                cj = c2;
                mj = m2;
                
                break;
                
            case 2:
                printf("calculating RR correlations\n");
                
                xc = false;
                
                ci = c2;
                ni = n2;
                
                cj = c2;
                mj = m2;
                
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
            mj = xc ? m2 : m1;
            Sj = xc ? S2 : S1;
        }
        
        st = time(NULL);
        dt = 0;
        fb = 0;
        
        ii = 0;
        nn = 0;
        
        #pragma omp parallel default(none) shared(st, dt, ii, nn, W, X, qQ) \
            private(i, j, nj) firstprivate(pt, xc, rd, ls, nd, db, gw, gh, \
                dx, dy, p, ni, ci, cj, mj, Si, Sj, stdout)
        {
            size_t qi;
            int q, nq;
            int* qr;
            
            double* Wi;
            double* Xi;
            
            nq = 0;
            qr = malloc((2*dy+1)*4*sizeof(int));
            if(!qr)
            {
                perror(NULL);
                abort();
            }
            
            Wi = calloc(nd, sizeof(double));
            Xi = calloc(4*nd, sizeof(double));
            if(!Wi || !Xi)
            {
                perror(NULL);
                abort();
            }
            
            #pragma omp master
            {
                signal(SIGALRM, handler);
                alarm(1);
            }
            
            qi = -1;
            
            #pragma omp for schedule(static, 1) nowait
            for(i = 0; i < ni; ++i)
            {
                if(qQ)
                    continue;
                
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
                
                const double sxi = ci[i*DW+0];
                const double syi = ci[i*DW+1];
                const double cxi = ci[i*DW+2];
                const double cyi = ci[i*DW+3];
                const double ui  = ci[i*DW+4];
                const double vi  = ci[i*DW+5];
                const double wi  = ci[i*DW+6];
                
                if(ci[i*DW+7] != qi)
                {
                    qi = ci[i*DW+7];
                    query(qi, gw, gh, dy, dx, &nq, qr);
                    for(q = 0; q < 2*nq; ++q)
                        qr[q] = mj[qr[q]];
                }
                
                for(q = 0; q < nq; ++q)
                {
                    j = qr[2*q+0];
                    nj = qr[2*q+1];
                    
                    if(!xc && j < i+1)
                        j = i+1;
                    
                    for(; j < nj; ++j)
                    {
                        const double sxj = cj[j*DW+0];
                        const double syj = cj[j*DW+1];
                        const double cxj = cj[j*DW+2];
                        const double cyj = cj[j*DW+3];
                        const double uj  = cj[j*DW+4];
                        const double vj  = cj[j*DW+5];
                        const double wj  = cj[j*DW+6];
                        
                        const double sdx = cxi*sxj - sxi*cxj;
                        const double cdx = cxi*cxj + rd*sxi*sxj;
                        
                        const double d1 = cdx*cyi - cyj;
                        const double d2 = sdx*cyi;
                        const double d3 = syj - syi;
                        
                        double d = d1*d1 + d2*d2 + d3*d3;
                        
                        if(d < db[0] || d >= db[nd])
                            continue;
                        
                        size_t n = 0;
                        while(db[n+1] <= d)
                            ++n;
                        
                        const double ww = wi*wj;
                        
                        Wi[n] += ww;
                        
                        if(!pt)
                        {
                            const double uu = ui*uj;
                            const double uv = ui*vj;
                            const double vu = vi*uj;
                            const double vv = vi*vj;
                            
                            const double xij = cyi*syj - syi*cyj*cdx;
                            const double yij = cyj*sdx;
                            const double xji = cyj*syi - syj*cyi*cdx;
                            const double yji = -cyi*sdx;
                            
                            double sij, cij, sji, cji;
                            nsincos(Si, xij, yij, &sij, &cij);
                            nsincos(Sj, xji, yji, &sji, &cji);
                            
                            const double cp = cij*cji + sij*sji;
                            const double sp = sij*cji - cij*sji;
                            const double cm = cij*cji - sij*sji;
                            const double sm = sij*cji + cij*sji;
                            
                            const double xip_re = (uu + vv)*cp - (uv - vu)*sp;
                            const double xip_im = (uv - vu)*cp + (uu + vv)*sp;
                            const double xim_re = (uu - vv)*cm - (uv + vu)*sm;
                            const double xim_im = (uv + vu)*cm + (uu - vv)*sm;
                            
                            Xi[0*nd+n] += (ww/Wi[n])*(xip_re - Xi[0*nd+n]);
                            Xi[1*nd+n] += (ww/Wi[n])*(xim_re - Xi[1*nd+n]);
                            Xi[2*nd+n] += (ww/Wi[n])*(xip_im - Xi[2*nd+n]);
                            Xi[3*nd+n] += (ww/Wi[n])*(xim_im - Xi[3*nd+n]);
                        }
                        
                        #pragma omp atomic
                        nn += 1;
                    }
                }
            }
            
            #pragma omp master
            {
                signal(SIGALRM, SIG_IGN);
            }
            
            #pragma omp critical
            for(i = 0; i < nd; ++i)
            {
                W[p*nd+i] += Wi[i];
                X[0*nd+i] += (Wi[i]/W[p*nd+i])*(Xi[0*nd+i] - X[0*nd+i]);
                X[1*nd+i] += (Wi[i]/W[p*nd+i])*(Xi[1*nd+i] - X[1*nd+i]);
                X[2*nd+i] += (Wi[i]/W[p*nd+i])*(Xi[2*nd+i] - X[2*nd+i]);
                X[3*nd+i] += (Wi[i]/W[p*nd+i])*(Xi[3*nd+i] - X[3*nd+i]);
            }
            
            free(qr);
            free(Wi);
            free(Xi);
        }
        
        dt = difftime(time(NULL), st);
        
        printf("\r> done with %zu pairs", nn);
        printf(" in %02d:%02d:%02d  \n", dt/3600, (dt/60)%60, dt%60);
        printf("\n");
    }
    
    fp = fopen(cfg.output, "w");
    if(!fp)
    {
        perror(cfg.output);
        return EXIT_FAILURE;
    }
    
    if(pt)
        fprintf(fp, "%-25s %-25s %-25s %-25s %-25s\n",
                                        "# theta", "xi", "DD", "DR", "RR");
    else
        fprintf(fp, "%-25s %-25s %-25s %-25s %-25s\n",
                                "# theta", "xip", "xim", "xip_im", "xim_im");
    
    for(i = 0; i < nd; ++i)
    {
        double d;
        
        if(ls)
            d = exp(log(dl) + (i + 0.5)*(log(dh) - log(dl))/nd);
        else
            d = dl + (i + 0.5)*(dh - dl)/nd;
        
        d /= uo;
        
        if(pt)
        {
            size_t ndd, ndr, nrr;
            double dd, dr, rr, xi;
            
            ndd = n1*(n1-1)/2;
            ndr = n1*n2;
            nrr = n2*(n2-1)/2;
            
            dd = W[0*nd+i]/ndd;
            dr = W[1*nd+i]/ndr;
            rr = W[2*nd+i]/nrr;
            
            xi = (dd - 2*dr + rr)/rr;
            
            fprintf(fp, "% .18e % .18e % .18e % .18e % .18e \n",
                                                        d, xi, dd, dr, rr);
        }
        else
        {
            double xip_re, xim_re, xip_im, xim_im;
            
            xip_re = X[0*nd+i];
            xim_re = X[1*nd+i];
            xip_im = X[2*nd+i];
            xim_im = X[3*nd+i];
            
            fprintf(fp, "% .18e % .18e % .18e % .18e % .18e\n",
                                        d, xip_re, xim_re, xip_im, xim_im);
        }
    }
    
    fclose(fp);
    
    free(W);
    free(X);
    free(m1);
    free(m2);
    free(c1);
    free(c2);
    free(dx);
    
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

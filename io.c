#include <string.h>

enum { SPACING_LIN, SPACING_LOG, NUM_SPACING };
char* CFG_SPACING[] = { "lin", "log" };
char* PRN_SPACING[] = { "linear", "logarithmic" };

enum { COORDS_FLAT, COORDS_LONLAT, NUM_COORDS };
char* CFG_COORDS[] = { "flat", "lonlat" };
char* PRN_COORDS[] = { "flat", "lon, lat" };

enum { FIELD_REAL, FIELD_COMPLEX, NUM_FIELD };
char* CFG_FIELD[] = { "real", "complex" };
char* PRN_FIELD[] = { "real", "complex" };

enum { SIGNS_PP, SIGNS_PM, SIGNS_MP, SIGNS_MM, NUM_SIGNS };
char* CFG_SIGNS[] = { "++", "+-", "-+", "--" };
char* PRN_SIGNS[] = { "u + i v", "u - i v", "-u + i v", "-u - i v" };

enum { TDATA_SHARE, TDATA_COPY, NUM_TDATA };
char* CFG_TDATA[] = { "share", "copy" };
char* PRN_TDATA[] = { "share", "copy" };

enum { UNIT_RAD, UNIT_DEG, UNIT_ARCMIN, UNIT_ARCSEC, NUM_UNIT };
char* CFG_UNIT[] = { "rad", "deg", "arcmin", "arcsec" };
char* PRN_UNIT[] = { "radian", "degree", "arcmin", "arcsec" };

const double UCONV[] = {
    1.0,
    0.017453292519943295769,
    0.00029088820866572159615,
    0.0000048481368110953599359
};

#ifndef LINELEN
#define LINELEN 1024
#endif

struct config {
    char* catalog1;
    char* catalog2;
    int dunit;
    int coords;
    int field1;
    int field2;
    int spin1;
    int spin2;
    int signs1;
    int signs2;
    char* output;
    int nth;
    double thmin;
    double thmax;
    int thunit;
    int spacing;
    int num_threads;
    int thread_data;
};

static int findstr(const char* s, int n, char* v[])
{
    int i;
    for(i = 0; i < n; ++i)
        if(strcmp(s, v[i]) == 0)
            break;
    return i;
}

void readcfg(const char* f, struct config* cfg)
{
    FILE* fp;
    char buf[LINELEN];
    size_t l;
    char* key;
    char* val;
    
    fp = fopen(f, "r");
    if(!fp)
        goto err_fopen;
    
    for(l = 1; fgets(buf, sizeof(buf), fp); ++l)
    {
        key = strtok(buf, " \t\r\n");
        val = strtok(NULL, " \t\r\n");
        
        if(!key || *key == '#')
            continue;
        if(!val || *val == '#')
            goto err_no_value;
        
        if(strcmp(key, "catalog") == 0)
            cfg->catalog1 = strdup(val);
        else if(strcmp(key, "catalog1") == 0)
            cfg->catalog1 = strdup(val);
        else if(strcmp(key, "catalog2") == 0)
            cfg->catalog2 = strdup(val);
        else if(strcmp(key, "dunit") == 0)
        {
            cfg->dunit = findstr(val, NUM_UNIT, CFG_UNIT);
            if(cfg->dunit == NUM_UNIT)
                goto err_bad_value;
        }
        else if(strcmp(key, "coords") == 0)
        {
            cfg->coords = findstr(val, NUM_COORDS, CFG_COORDS);
            if(cfg->coords == NUM_COORDS)
                goto err_bad_value;
        }
        else if(strcmp(key, "field") == 0)
        {
            cfg->field1 = cfg->field2 = findstr(val, NUM_FIELD, CFG_FIELD);
            if(cfg->field1 == NUM_FIELD)
                goto err_bad_value;
        }
        else if(strcmp(key, "field1") == 0)
        {
            cfg->field1 = findstr(val, NUM_FIELD, CFG_FIELD);
            if(cfg->field1 == NUM_FIELD)
                goto err_bad_value;
        }
        else if(strcmp(key, "field2") == 0)
        {
            cfg->field2 = findstr(val, NUM_FIELD, CFG_FIELD);
            if(cfg->field2 == NUM_FIELD)
                goto err_bad_value;
        }
        else if(strcmp(key, "spin") == 0)
        {
            cfg->spin1 = cfg->spin2 = atoi(val);
            if(cfg->spin1 < 0)
                goto err_bad_value;
            cfg->field1 = cfg->field2 = FIELD_COMPLEX;
        }
        else if(strcmp(key, "spin1") == 0)
        {
            cfg->spin1 = atoi(val);
            if(cfg->spin1 < 0)
                goto err_bad_value;
            cfg->field1 = FIELD_COMPLEX;
        }
        else if(strcmp(key, "spin2") == 0)
        {
            cfg->spin2 = atoi(val);
            if(cfg->spin2 < 0)
                goto err_bad_value;
            cfg->field2 = FIELD_COMPLEX;
        }
        else if(strcmp(key, "signs") == 0)
        {
            cfg->signs1 = cfg->signs2 = findstr(val, NUM_SIGNS, CFG_SIGNS);
            if(cfg->signs1 == NUM_SIGNS)
            {
                if(strcmp(val, "+") == 0)
                    cfg->signs1 = cfg->signs2 = SIGNS_PP;
                else if(strcmp(val, "-") == 0)
                    cfg->signs1 = cfg->signs2 = SIGNS_MM;
                else
                    goto err_bad_value;
            }
        }
        else if(strcmp(key, "signs1") == 0)
        {
            cfg->signs1 = findstr(val, NUM_SIGNS, CFG_SIGNS);
            if(cfg->signs1 == NUM_SIGNS)
            {
                if(strcmp(val, "+") == 0)
                    cfg->signs1 = SIGNS_PP;
                else if(strcmp(val, "-") == 0)
                    cfg->signs1 = SIGNS_MM;
                else
                    goto err_bad_value;
            }
        }
        else if(strcmp(key, "signs2") == 0)
        {
            cfg->signs2 = findstr(val, NUM_SIGNS, CFG_SIGNS);
            if(cfg->signs2 == NUM_SIGNS)
            {
                if(strcmp(val, "+") == 0)
                    cfg->signs2 = SIGNS_PP;
                else if(strcmp(val, "-") == 0)
                    cfg->signs2 = SIGNS_MM;
                else
                    goto err_bad_value;
            }
        }
        else if(strcmp(key, "output") == 0)
            cfg->output = strdup(val);
        else if(strcmp(key, "nth") == 0)
        {
            cfg->nth = atoi(val);
            if(cfg->nth <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "thmin") == 0)
        {
            cfg->thmin = atof(val);
            if(cfg->thmin <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "thmax") == 0)
        {
            cfg->thmax = atof(val);
            if(cfg->thmax <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "thunit") == 0)
        {
            cfg->thunit = findstr(val, NUM_UNIT, CFG_UNIT);
            if(cfg->thunit == NUM_UNIT)
                goto err_bad_value;
        }
        else if(strcmp(key, "spacing") == 0)
        {
            cfg->spacing = findstr(val, NUM_SPACING, CFG_SPACING);
            if(cfg->spacing == NUM_SPACING)
                goto err_bad_value;
        }
        else if(strcmp(key, "num_threads") == 0)
        {
            cfg->num_threads = atoi(val);
            if(cfg->num_threads <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "thread_data") == 0)
        {
            cfg->thread_data = findstr(val, NUM_TDATA, CFG_TDATA);
            if(cfg->thread_data == NUM_TDATA)
                goto err_bad_value;
        }
        else
            goto err_bad_key;
    }
    
    fclose(fp);
    
    if(!cfg->catalog1)
        { key = "catalog"; goto err_missing_key; }
    if(!cfg->output)
        { key = "output"; goto err_missing_key; }
    if(!cfg->nth)
        { key = "nth"; goto err_missing_key; }
    if(!cfg->thmin)
        { key = "thmin"; goto err_missing_key; }
    if(!cfg->thmax)
        { key = "thmax"; goto err_missing_key; }
    
    return;
    
err_fopen:
    perror(f);
    exit(EXIT_FAILURE);
    
err_bad_key:
    fprintf(stderr, "error: %s:%zu: invalid key `%s`\n", f, l, key);
    exit(EXIT_FAILURE);
    
err_missing_key:
    fprintf(stderr, "error: %s: missing required `%s` key\n", f, key);
    exit(EXIT_FAILURE);
    
err_no_value:
    fprintf(stderr, "error: %s:%zu: missing value\n", f, l);
    exit(EXIT_FAILURE);
    
err_bad_value:
    fprintf(stderr, "error: %s:%zu: invalid value `%s`\n", f, l, val);
    exit(EXIT_FAILURE);
}

void freecfg(struct config* cfg)
{
    free(cfg->catalog1);
    free(cfg->catalog2);
    free(cfg->output);
}

void printcfg(const struct config* cfg)
{
    if(!cfg->catalog2)
        printf("catalog ......... %s\n", cfg->catalog1);
    else
    {
        printf("catalog 1 ....... %s\n", cfg->catalog1);
        printf("catalog 2 ....... %s\n", cfg->catalog2);
    }
    printf("data units ...... %s\n", PRN_UNIT[cfg->dunit]);
    printf("coordinates ..... %s\n", PRN_COORDS[cfg->coords]);
    if(!cfg->catalog2)
    {
        if(cfg->field1)
        {
            printf("field type ...... complex spin-%d\n", cfg->spin1);
            printf("signature ....... %s\n", PRN_SIGNS[cfg->signs1]);
        }
        else
        {
            printf("field type ...... real\n");
            printf("signature ....... %su\n", cfg->signs1&2 ? "-" : "+");
        }
    }
    else
    {
        if(cfg->field1)
        {
            printf("field 1 type .... complex spin-%d\n", cfg->spin1);
            printf("signature 1 ..... %s\n", PRN_SIGNS[cfg->signs1]);
        }
        else
        {
            printf("field 1 type .... real\n");
            printf("signature 1 ..... %su\n", cfg->signs1&2 ? "-" : "+");
        }
        if(cfg->field2)
        {
            printf("field 2 type .... complex spin-%d\n", cfg->spin2);
            printf("signature 2 ..... %s\n", PRN_SIGNS[cfg->signs2]);
        }
        else
        {
            printf("field 2 type .... real\n");
            printf("signature 2 ..... %su\n", cfg->signs2&2 ? "-" : "+");
        }
    }
    printf("\n");
    printf("output file ..... %s\n", cfg->output);
    printf("num. points ..... %u\n", cfg->nth);
    printf("point range ..... %lg to %lg %s\n", cfg->thmin, cfg->thmax,
                                                    PRN_UNIT[cfg->thunit]);
    printf("point spacing ... %s\n", PRN_SPACING[cfg->spacing]);
#ifdef _OPENMP
    printf("\n");
    printf("num. threads .... %d\n", omp_get_max_threads());
    printf("thread data ..... %s\n", PRN_TDATA[cfg->thread_data]);
#endif
}

double* readc(const char* f, int co, double ui, bool cf, int sg, int* n)
{
    FILE* fp;
    int l;
    char buf[LINELEN];
    int i, a;
    double* d;
    char* sx;
    char* sy;
    char* su;
    char* sv;
    char* sw;
    double x, y, u, v, w;
    
    if(strcmp(f, "-") == 0)
        fp = stdin;
    else
    {
        fp = fopen(f, "r");
        if(!fp)
            goto err_fopen;
    }
    
    i = 0;
    a = 1;
    
    d = malloc(a*DW*sizeof(double));
    if(!d)
        goto err_alloc;
    
    for(l = 1; fgets(buf, sizeof buf, fp); ++l)
    {
        sx = strtok(buf, " \t\r\n");
        sy = strtok(NULL, " \t\r\n");
        su = strtok(NULL, " \t\r\n");
        sv = cf ? strtok(NULL, " \t\r\n") : NULL;
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
        
        if(!sw || *sw == '#')
            w = 1;
        else
            w = atof(sw);
        
        switch(co)
        {
        case COORDS_FLAT:
            d[i*DW+0] = 1;
            d[i*DW+1] = x;
            d[i*DW+2] = y;
            break;
        case COORDS_LONLAT:
            d[i*DW+0] = cos(x)*cos(y);
            d[i*DW+1] = sin(x)*cos(y);
            d[i*DW+2] = sin(y);
            break;
        }
        
        d[i*DW+3] = sg&2 ? -u : u;
        d[i*DW+4] = sg&1 ? -v : v;
        d[i*DW+5] = w;
        
        i += 1;
        
        if(i == a)
        {
            a *= 2;
            d = realloc(d, a*DW*sizeof(double));
            if(!d)
                goto err_alloc;
        }
    }
    
    if(fp != stdin)
        fclose(fp);
    
    d = realloc(d, i*DW*sizeof(double));
    if(!d)
        goto err_alloc;
    
    if(n)
        *n = i;
    
    return d;
    
err_fopen:
    perror(f);
    exit(EXIT_FAILURE);
    
err_alloc:
    perror(NULL);
    exit(EXIT_FAILURE);
}

void output(const char* f, size_t n, double a, double b, double u,
                bool sc, bool ls, double* N, double* X, double* W, double* Y)
{
    FILE* fp;
    size_t i;
    
    double d;
    
    if(strcmp(f, "-") == 0)
        fp = stdout;
    else
    {
        fp = fopen(f, "w");
        if(!fp)
            goto err_fopen;
    }
    
    fprintf(fp, "# %-24s", "theta");
    fprintf(fp, "  %-24s", "xip");
    fprintf(fp, "  %-24s", "xim");
    fprintf(fp, "  %-24s", "xip_im");
    fprintf(fp, "  %-24s", "xim_im");
    fprintf(fp, "  %-24s", "np");
    fprintf(fp, "  %-24s", "mat_dia");
    fprintf(fp, "  %-24s", "mat_sub");
    fprintf(fp, "  %-24s", "rhs_p");
    fprintf(fp, "  %-24s", "rhs_m");
    fprintf(fp, "  %-24s", "rhs_p_im");
    fprintf(fp, "  %-24s", "rhs_m_im");
    fprintf(fp, "\n");
    
    for(i = 0; i < n; ++i)
    {
        if(ls)
            d = exp(log(a) + i*(log(b) - log(a))/(n - 1));
        else
            d = a + i*(b - a)/(n - 1);
        
        if(sc)
            d = 2*asin(0.5*d);
        
        fprintf(fp, " % .18e", d/u);
        
        fprintf(fp, " % .18e", X[0*n+i]);
        fprintf(fp, " % .18e", X[1*n+i]);
        fprintf(fp, " % .18e", X[2*n+i]);
        fprintf(fp, " % .18e", X[3*n+i]);
        
        fprintf(fp, " % .18e", N[i]);
        
        fprintf(fp, " % .18e", W[0*n+i]);
        fprintf(fp, " % .18e", W[1*n+i]);
        
        fprintf(fp, " % .18e", Y[0*n+i]);
        fprintf(fp, " % .18e", Y[1*n+i]);
        fprintf(fp, " % .18e", Y[2*n+i]);
        fprintf(fp, " % .18e", Y[3*n+i]);
        
        fprintf(fp, "\n");
    }
    
    if(fp != stdout)
        fclose(fp);
    
    return;
    
err_fopen:
    perror(f);
    exit(EXIT_FAILURE);
}

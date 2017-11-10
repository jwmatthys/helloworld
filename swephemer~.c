#include "m_pd.h"
#include "m_imp.h"
#include "swephexp.h"   /* this includes  "sweodef.h" */
#include <regex.h>

#define UNUSED(x) (void)(x)
typedef enum { false, true } bool;

static t_class *swephemer_tilde_class;

typedef struct _swephemer_tilde {
        t_object x_obj;
        t_sample x_sample_rate;
        double x_jul_date;
        double x_loop_start;
        double x_loop_end;
        double x_loop_val;
        t_float x_step;
        long x_body;
        bool x_dsp_flag;
        bool x_loop;
        long x_iflag;
        t_outlet* x_longitude_out, *x_latitude_out, *x_distance_out, *x_date_out;
        t_outlet* x_signal_out;
} t_swephemer_tilde;

void *swephemer_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
        UNUSED(s);
        t_swephemer_tilde *x = (t_swephemer_tilde *)pd_new(swephemer_tilde_class);
        swe_set_ephe_path(NULL); // tell Pd to look for ephemera files in same dir as this external
        char *svers = malloc(1024);
        swe_version(svers);
        post ("swephemer~ v0.5-alpha2");
        post ("Swiss Ephemeris version %s", svers);
        unsigned short i;
        char snam[40];
        for (i = 0; i < argc; i++)
        {
                switch (argv[i].a_type)
                {
                case A_FLOAT:
                        x->x_body = (long)argv[i].a_w.w_float;
                        swe_get_planet_name(x->x_body, snam);
                        post ("Celestial body: %s", snam);
                        break;
                case A_SYMBOL:
                        //error("Ignoring %s", argv[i].a_w.w_symbol->s_name);
                        break;
                default:
                        error("Got something unexpected...");
                }
        }
        x->x_dsp_flag = false;
        x->x_sample_rate = sys_getsr();
        outlet_new(&x->x_obj, &s_signal); //gensym("signal"));
        x->x_longitude_out = outlet_new(&x->x_obj, &s_float);
        x->x_latitude_out = outlet_new(&x->x_obj, &s_float);
        x->x_distance_out = outlet_new(&x->x_obj, &s_float);
        x->x_date_out = outlet_new(&x->x_obj, &s_float);
        x->x_body = 0;
        x->x_iflag = 0;
        x->x_step = 1.0; // one day
        x->x_jul_date = 0;
        x->x_loop_start = 0;
        x->x_loop_end = 0;
        x->x_loop_val = 0;
        x->x_loop = false;
        return (void *)x;
}

void swephemer_tilde_free(t_swephemer_tilde *x)
{
        outlet_free(x->x_longitude_out);
        outlet_free(x->x_latitude_out);
        outlet_free(x->x_distance_out);
        outlet_free(x->x_date_out);
        //outlet_free(x->x_signal_out);
        swe_close();
}

static t_int *swephemer_tilde_perform(t_int *w)
{
        t_swephemer_tilde *x = (t_swephemer_tilde *)(w[1]);
        t_sample *out = (t_sample*)(w[2]);
        double x2[6];
        char serr[256];
        int n = (int)(w[3]);
        long iflgret = 0;
        t_sample outval = 0;
        while (n--)
        {
                if (x->x_dsp_flag == false) outval = 0;
                else
                {
                        iflgret = swe_calc(x->x_loop_val, x->x_body, x->x_iflag, x2, serr);
                        if (iflgret < 0) x->x_dsp_flag = false;
                        t_float val = (t_float)sin (3.14159265 * x2[0] / 180);
                        //post ("val: %f", val);
                        outval = (t_sample)val;
                        x->x_loop_val += x->x_step;
                        if (x->x_loop == true && x->x_loop_val > x->x_loop_end) x->x_loop_val = x->x_loop_start;
                }
                *out++ = outval;
        }
        if (iflgret < 0) error("%s", serr);
        return (w+4);
}

static void swephemer_tilde_dsp(t_swephemer_tilde *x, t_signal **sp)
{
        dsp_add(swephemer_tilde_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}

double atom2double (t_atom *argv, t_int index)
{
        unsigned short len;
        char *subbuf;
        double result = 0;
        switch (argv[index].a_type)
        {
        case A_FLOAT:
                result = (double)argv[index].a_w.w_float;
                break;
        case A_SYMBOL:
                len = strlen(argv[index].a_w.w_symbol->s_name);
                subbuf = malloc(len-2);
                sprintf(subbuf, "%.*s", len-2, argv[index].a_w.w_symbol->s_name+1);
                result = atof(subbuf);
                free (subbuf);
                break;
        default:
                error ("Unexpected input");
        }
        return result;
}

void swephemer_tilde_bang(t_swephemer_tilde *x)
{
        double x2[6];
        char serr[256];
        long iflgret = swe_calc(x->x_jul_date, x->x_body, x->x_iflag, x2, serr);
        if (iflgret < 0) error("%s", serr);
        else
        {
                outlet_float (x->x_longitude_out, x2[0]);
                outlet_float (x->x_latitude_out, x2[1]);
                outlet_float (x->x_distance_out, x2[2]);
                outlet_float (x->x_date_out, x->x_jul_date);
        }
}

void swephemer_tilde_float(t_swephemer_tilde *x, t_float f)
{
        x->x_jul_date += f;
        //post ("Received adjustment value %+g. New date is %10.4f", f, x->x_jul_date);
        swephemer_tilde_bang(x);
}

void swephemer_tilde_step(t_swephemer_tilde *x, t_float f)
{
        x->x_step = f;
        post ("Set step to %g", x->x_step);
}

void swephemer_tilde_loop(t_swephemer_tilde *x, t_symbol *s, short argc, t_atom *argv)
{
        UNUSED(s);
        switch (argc)
        {
        case 3:
                x->x_step = (double)argv[2].a_w.w_float;
        case 2:
                x->x_loop_end = atom2double(argv,1);
        case 1:
                x->x_loop_start = atom2double(argv,0);
                x->x_loop_val = x->x_loop_start;
                x->x_loop_end = x->x_loop_start + 365242.2; // 1000 years
                x->x_loop = true;
                break;
        case 0:
                x->x_loop = true;
                break;
        default:
                error("too many arguments to -loop");
        }
        post ("Looping from %10.2f to %10.2f by %g step", x->x_loop_start, x->x_loop_end, x->x_step);
}

void swephemer_tilde_loopoff(t_swephemer_tilde *x)
{
        x->x_loop = false;
        post("looping disabled");
}

void swephemer_tilde_bj(t_swephemer_tilde *x, t_symbol *s, short argc, t_atom *argv)
{
        UNUSED(s);
        UNUSED(argc);
        x->x_jul_date = atom2double(argv,0);
        x->x_loop_val = x->x_jul_date;
        post ("Set date to %10.4f", x->x_jul_date);
        outlet_float (x->x_date_out, x->x_jul_date);
}

void swephemer_tilde_b(t_swephemer_tilde *x, t_symbol *s, short argc, t_atom *argv)
{
        UNUSED(s);
        long y = (argc > 0 && argv[0].a_type == A_FLOAT) ? argv[0].a_w.w_float : 0;
        long m = (argc > 1 && argv[1].a_type == A_FLOAT) ? argv[1].a_w.w_float : 1;
        long d = (argc > 2 && argv[2].a_type == A_FLOAT) ? argv[2].a_w.w_float : 1;
        double h = (argc > 3 && argv[3].a_type == A_FLOAT) ? argv[3].a_w.w_float : 0;
        x->x_jul_date = swe_julday(y,m,d,h,SE_GREG_CAL);
        x->x_loop_val = x->x_jul_date;
        post ("Set %d.%d.%d.%f to Julian date %10.4f", y, m, d, h, x->x_jul_date);
        outlet_float (x->x_date_out, x->x_jul_date);
}

void swephemer_tilde_path(t_swephemer_tilde *x, t_symbol *s)
{
        UNUSED(x);
        swe_set_ephe_path(s->s_name); // tell Pd to look for ephemera files in same dir as this external
        post("Swephemer will look for ephemera files at %s", s->s_name);
}

void swephemer_tilde_body(t_swephemer_tilde *x, t_float f)
{
        x->x_body = (long)f;
        char snam[40];
        swe_get_planet_name(x->x_body, snam);
        post ("Celestial body: %s", snam);
}

void swephemer_tilde_audioflag(t_swephemer_tilde *x, t_float f)
{
        int v = (int)f;
        x->x_dsp_flag = (v == 0) ? false : true;
}

void swephemer_tilde_iflag(t_swephemer_tilde *x, t_symbol *s, short argc, t_atom *argv)
{
        UNUSED(s);
        x->x_iflag = 0;
        unsigned short i;
        for (i = 0; i < argc; i++)
        {
                switch(argv[i].a_type)
                {
                case A_FLOAT:
                        x->x_iflag += (int)argv[i].a_w.w_float;
                        break;
                case A_SYMBOL:
                        if (argv[i].a_w.w_symbol == gensym("DEFAULT"))
                        {
                                x->x_iflag = 0;
                                post ("resetting iflags");
                                break;
                        }
                        if (argv[i].a_w.w_symbol == gensym("JPLEPH"))
                        {
                                x->x_iflag += SEFLG_JPLEPH;
                                post ("adding JPLEPH flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("SWIEPH"))
                        {
                                x->x_iflag += SEFLG_SWIEPH;
                                post ("adding SWIEPH flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("MOSEPH"))
                        {
                                x->x_iflag += SEFLG_MOSEPH;
                                post ("adding MOSEPH flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("HELCTR"))
                        {
                                x->x_iflag += SEFLG_HELCTR;
                                post ("adding HELCTR flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("TRUEPOS"))
                        {
                                x->x_iflag += SEFLG_TRUEPOS;
                                post ("adding TRUEPOS flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("J2000"))
                        {
                                x->x_iflag += SEFLG_J2000;
                                post ("adding J2000 flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("NONUT"))
                        {
                                x->x_iflag += SEFLG_NONUT;
                                post ("adding NONUT flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("SPEED3"))
                        {
                                x->x_iflag += SEFLG_SPEED3;
                                post ("adding SPEED3 flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("SPEED"))
                        {
                                x->x_iflag += SEFLG_SPEED;
                                post ("adding SPEED flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("NOGDEFL"))
                        {
                                x->x_iflag += SEFLG_NOGDEFL;
                                post ("adding NOGDEFL flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("NOABERR"))
                        {
                                x->x_iflag += SEFLG_NOABERR;
                                post ("adding NOABERR flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("ASTROMETRIC"))
                        {
                                x->x_iflag += SEFLG_ASTROMETRIC;
                                post ("adding ASTROMETRIC flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("EQUATORIAL"))
                        {
                                x->x_iflag += SEFLG_EQUATORIAL;
                                post ("adding EQUATORIAL flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("XYZ"))
                        {
                                x->x_iflag += SEFLG_XYZ;
                                post ("adding XYZ flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("RADIANS"))
                        {
                                x->x_iflag += SEFLG_RADIANS;
                                post ("adding RADIANS flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("BARYCTR"))
                        {
                                x->x_iflag += SEFLG_BARYCTR;
                                post ("adding BARYCTR flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("TOPOCTR"))
                        {
                                x->x_iflag += SEFLG_TOPOCTR;
                                post ("adding TOPOCTR flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("SIDEREAL"))
                        {
                                x->x_iflag += SEFLG_SIDEREAL;
                                post ("adding SIDEREAL flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("ICRS"))
                        {
                                x->x_iflag += SEFLG_ICRS;
                                post ("adding ICRS flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("DPSIDEPS_1980"))
                        {
                                x->x_iflag += SEFLG_DPSIDEPS_1980;
                                post ("adding DPSIDEPS_1980 flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("JPLHOR"))
                        {
                                x->x_iflag += SEFLG_JPLHOR;
                                post ("adding JPLHOR flag");
                        }
                        else if (argv[i].a_w.w_symbol == gensym("JPLHOR_APPROX"))
                        {
                                x->x_iflag += SEFLG_JPLHOR_APPROX;
                                post ("adding JPLHOR_APPROX flag");
                        }
                        else error ("unrecognized iflag: %s", argv[i].a_w.w_symbol->s_name);
                        break;
                default:
                        error ("unrecognized type for -iflag");
                }
        }
}

void swephemer_tilde_setup(void) {
        swephemer_tilde_class = class_new(gensym("swephemer~"),
                                          (t_newmethod)swephemer_tilde_new,
                                          (t_method)swephemer_tilde_free,
                                          sizeof(t_swephemer_tilde), 0, A_GIMME, 0);
        class_addbang(swephemer_tilde_class, swephemer_tilde_bang);
        class_addfloat(swephemer_tilde_class, swephemer_tilde_float);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_bj, gensym("-bj"), A_GIMME, 0);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_bj, gensym("-j"), A_GIMME, 0);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_b, gensym("-b"), A_GIMME, 0);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_body, gensym("-p"), A_FLOAT, 0);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_step, gensym("-step"), A_FLOAT, 0);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_path, gensym("-path"), A_SYMBOL, 0);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_iflag, gensym("-iflag"), A_GIMME, 0);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_audioflag, gensym("-audio"), A_FLOAT, 0);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_loop, gensym("-loop"), A_GIMME, 0);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_loopoff, gensym("-loopoff"), 0);
        class_addmethod(swephemer_tilde_class, (t_method)swephemer_tilde_dsp, gensym("dsp"), 0);
}

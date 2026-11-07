/* DUO Aquarium: experimental visualization of interacting three-bit cycles.
 * Standalone demonstration. No files from the DUE library are modified.
 * Copyright (c) 2026. Released under the MIT license (see LICENSE).
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <time.h>
#ifdef WITH_X11
#include <X11/Xlib.h>
#include <X11/keysym.h>
#endif

static const uint8_t DUO[2][4] = {{0,6,3,2},{7,1,4,5}};
static const uint8_t ORBIT_OF[8] = {0,1,0,0,1,1,0,1};
static const uint8_t PHASE_OF[8] = {0,1,3,2,2,3,1,0};

typedef struct {
    unsigned orbit, phase;
    double clock, spin_rate;
    double x,y,vx,vy;
    unsigned long collisions;
} Particle;
typedef struct {
    Particle *p;
    unsigned count;
    double width, height;
    uint64_t rng;
    unsigned long long tick, impacts;
    FILE *events, *samples;
    unsigned sample_every;
    bool periodic;
} World;

static uint64_t rng_next(World *w) {
    uint64_t x = w->rng;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    w->rng = x;
    return x * UINT64_C(2685821657736338717);
}
static double uniform(World *w) { return (rng_next(w) >> 11) * 0x1.0p-53; }
static double random_range(World *w,double lo,double hi){return lo+(hi-lo)*uniform(w);}
static unsigned state_of(const Particle *p){return DUO[p->orbit][p->phase];}
static void assign_state(Particle *p, unsigned state) {
    p->orbit=ORBIT_OF[state]; p->phase=PHASE_OF[state];
    /* A collision installs the XOR pattern; the independent fractional clocks persist. */
}
static void sample(World *w) {
    if (!w->samples) return;
    unsigned orbit_a=0, states[8]={0};
    double sum=0,sum2=0,min=INFINITY,max=-INFINITY;
    for (unsigned i=0;i<w->count;i++) {
        Particle *p=&w->p[i];
        orbit_a += p->orbit==0;
        states[state_of(p)]++;
        sum += p->spin_rate;
        sum2 += p->spin_rate*p->spin_rate;
        if(p->spin_rate<min)min=p->spin_rate;
        if(p->spin_rate>max)max=p->spin_rate;
    }
    double avg=sum/w->count, variance=fmax(0,sum2/w->count-avg*avg);
    fprintf(w->samples,"%llu,%llu,%u,%u,%.9f,%.9f,%.9f,%.9f,%.9f",
        w->tick,w->impacts,orbit_a,w->count-orbit_a,sum,avg,sqrt(variance),min,max);
    for(unsigned s=0;s<8;s++)fprintf(w->samples,",%u",states[s]);
    fputc('\n',w->samples);
}
static void collision(World *w,unsigned i,unsigned j) {
    Particle *a=&w->p[i], *b=&w->p[j];
    const unsigned before_a=state_of(a),before_b=state_of(b);
    const double ea=a->spin_rate,eb=b->spin_rate;
    const unsigned result=before_a^before_b;
    /* Transfer 1/8 of the signed difference: conserve sum, don't fully average. */
    const double transfer=(ea-eb)*0.125;
    a->spin_rate=ea-transfer;
    b->spin_rate=eb+transfer;
    assign_state(a,result);assign_state(b,result);
    w->impacts++;a->collisions++;b->collisions++;
    if(w->events) fprintf(w->events,
      "%llu,%llu,%u,%u,%u,%u,%u,%u,%u,%.9f,%.9f,%.9f,%.9f,%.4f,%.4f\n",
      w->tick,w->impacts,i,j,before_a,before_b,result,a->orbit,a->phase,
      ea,eb,a->spin_rate,b->spin_rate,(a->x+b->x)/2,(a->y+b->y)/2);
}
#define RADIUS 8.0
static void update(World *w) {
    const double dt=1.0/60.0;
    w->tick++;
    for(unsigned i=0;i<w->count;i++){
        Particle *p=&w->p[i];
        p->clock += p->spin_rate*dt;
        unsigned advance=(unsigned)p->clock;
        if(advance){p->phase=(p->phase+advance)%4;p->clock-=advance;}
        p->x+=p->vx*dt;p->y+=p->vy*dt;
        if(w->periodic){
            if(p->x<0)p->x+=w->width;
            if(p->x>=w->width)p->x-=w->width;
            if(p->y<0)p->y+=w->height;
            if(p->y>=w->height)p->y-=w->height;
        }else{
            if(p->x<RADIUS){p->x=RADIUS;p->vx=fabs(p->vx);}
            if(p->x>w->width-RADIUS){p->x=w->width-RADIUS;p->vx=-fabs(p->vx);}
            if(p->y<RADIUS){p->y=RADIUS;p->vy=fabs(p->vy);}
            if(p->y>w->height-RADIUS){p->y=w->height-RADIUS;p->vy=-fabs(p->vy);}
        }
    }
    /* Pairwise equal-mass elastic collisions of motion; XOR changes internal state.
       Speeds of spatial travel and internal spin are distinct quantities. */
    for(unsigned i=0;i<w->count;i++)for(unsigned j=i+1;j<w->count;j++){
        Particle *a=&w->p[i],*b=&w->p[j];
        double dx=b->x-a->x,dy=b->y-a->y;
        if(w->periodic){
            if(dx>w->width/2)dx-=w->width;
            if(dx<-w->width/2)dx+=w->width;
            if(dy>w->height/2)dy-=w->height;
            if(dy<-w->height/2)dy+=w->height;
        }
        double dsq=dx*dx+dy*dy;
        if(dsq>=(2*RADIUS)*(2*RADIUS))continue;
        double dist=sqrt(dsq);
        if(dist<1.e-9){dx=1;dy=0;dist=1;}
        double nx=dx/dist,ny=dy/dist;
        double approach=(b->vx-a->vx)*nx+(b->vy-a->vy)*ny;
        double overlap=2*RADIUS-dist;
        a->x-=nx*overlap*.5;a->y-=ny*overlap*.5;
        b->x+=nx*overlap*.5;b->y+=ny*overlap*.5;
        if(w->periodic){
            if(a->x<0) a->x+=w->width;
            if(a->x>=w->width) a->x-=w->width;
            if(b->x<0) b->x+=w->width;
            if(b->x>=w->width) b->x-=w->width;
            if(a->y<0) a->y+=w->height;
            if(a->y>=w->height) a->y-=w->height;
            if(b->y<0) b->y+=w->height;
            if(b->y>=w->height) b->y-=w->height;
        }
        if(approach>=0)continue; /* Avoid repeated impact on already separating objects. */
        a->vx+=approach*nx;a->vy+=approach*ny;
        b->vx-=approach*nx;b->vy-=approach*ny;
        collision(w,i,j);
    }
    if(w->sample_every && w->tick%w->sample_every==0)sample(w);
}
static void world_init(World *w,unsigned count,uint64_t seed,double width,double height,bool periodic){
    memset(w,0,sizeof(*w));
    w->width=width;w->height=height;w->count=count;
    w->rng=seed?seed:UINT64_C(0x9e3779b97f4a7c15);
    w->periodic=periodic;
    w->sample_every=60;
    w->p=calloc(count,sizeof(Particle));
    if(!w->p){perror("calloc");exit(1);}
    for(unsigned i=0;i<count;i++){
        Particle *p=&w->p[i];
        assign_state(p,(unsigned)(rng_next(w)%8));
        p->clock=uniform(w);
        p->spin_rate=random_range(w,0.6,9.0);
        p->x=random_range(w,RADIUS*2,width-RADIUS*2);
        p->y=random_range(w,RADIUS*2,height-RADIUS*2);
        double ang=random_range(w,0,2*3.14159265358979323846);
        double speed=random_range(w,25,85);
        p->vx=cos(ang)*speed;p->vy=sin(ang)*speed;
    }
}
static void world_log(World *w,const char *prefix){
    char path[1024];
    snprintf(path,sizeof(path),"%s_events.csv",prefix);
    w->events=fopen(path,"w");
    if(!w->events){perror(path);exit(1);}
    fprintf(w->events,"tick,impact,id_a,id_b,state_a_before,state_b_before,xor_result,new_orbit,new_phase,spin_a_before,spin_b_before,spin_a_after,spin_b_after,impact_x,impact_y\n");
    snprintf(path,sizeof(path),"%s_samples.csv",prefix);
    w->samples=fopen(path,"w");
    if(!w->samples){perror(path);exit(1);}
    fprintf(w->samples,"tick,impacts,orbit_A,orbit_B,total_spin,mean_spin,stddev_spin,min_spin,max_spin,state_000,state_001,state_010,state_011,state_100,state_101,state_110,state_111\n");
    sample(w);
}
static void close_log(World *w){if(w->events)fclose(w->events);if(w->samples)fclose(w->samples);}
static void results(World *w){
    unsigned na=0;double sum=0,min=INFINITY,max=-INFINITY;
    for(unsigned i=0;i<w->count;i++){Particle *p=&w->p[i];na+=p->orbit==0;sum+=p->spin_rate;min=fmin(min,p->spin_rate);max=fmax(max,p->spin_rate);}
    printf("Ticks: %llu  Collisions: %llu  Orbit A/B: %u/%u  Total spin: %.9f  Spin range: %.4f..%.4f\n",w->tick,w->impacts,na,w->count-na,sum,min,max);
}
#ifdef WITH_X11
static unsigned long color(Display *d,int screen,const char *name){
    XColor c,exact;
    if(XAllocNamedColor(d,DefaultColormap(d,screen),name,&c,&exact))return c.pixel;
    return BlackPixel(d,screen);
}
static void draw(World *w,Display *d,Window win,GC gc,int screen,bool paused){
    unsigned long bg=color(d,screen,"#111b29"),white=color(d,screen,"#edf2f7");
    unsigned long blue=color(d,screen,"#54a4e1"),gold=color(d,screen,"#f1b75b");
    unsigned long cold=color(d,screen,"#5fcac5"),hot=color(d,screen,"#ec735f");
    XSetForeground(d,gc,bg);XFillRectangle(d,win,gc,0,0,(unsigned)w->width,(unsigned)w->height+66);
    for(unsigned i=0;i<w->count;i++){
        Particle *p=&w->p[i];int x=(int)p->x,y=(int)p->y;
        XSetForeground(d,gc,p->orbit?gold:blue);
        XFillArc(d,win,gc,x-8,y-8,16,16,0,360*64);
        XSetForeground(d,gc,p->spin_rate>4.8?hot:cold);
        XDrawArc(d,win,gc,x-11,y-11,22,22,0,360*64);
        XSetForeground(d,gc,white);
        char bit[8];snprintf(bit,sizeof(bit),"%u",state_of(p));
        XDrawString(d,win,gc,x+12,y-10,bit,(int)strlen(bit));
    }
    unsigned a=0;double sum=0;
    for(unsigned i=0;i<w->count;i++){a+=w->p[i].orbit==0;sum+=w->p[i].spin_rate;}
    char buf[256];snprintf(buf,sizeof(buf),"%s  ticks %llu | impacts %llu | orbit A %u / B %u | total spin %.2f",paused?"PAUSED |":"RUNNING |",w->tick,w->impacts,a,w->count-a,sum);
    XSetForeground(d,gc,white);XDrawString(d,win,gc,12,(int)w->height+23,buf,(int)strlen(buf));
    const char *controls="Space pause   N single-step   Esc/Q quit    Blue=A, Gold=B; ring=spin level";
    XDrawString(d,win,gc,12,(int)w->height+46,controls,(int)strlen(controls));
    XFlush(d);
}
static int run_window(World *w,unsigned long long steps){
    Display *d=XOpenDisplay(NULL);
    if(!d){fprintf(stderr,"Cannot open X11 display. Try --headless.\n");return 1;}
    int screen=DefaultScreen(d);
    Window win=XCreateSimpleWindow(d,RootWindow(d,screen),80,70,(unsigned)w->width,(unsigned)w->height+66,1,WhitePixel(d,screen),BlackPixel(d,screen));
    XStoreName(d,win,"DUO Aquarium — three-bit demonstration");
    XSelectInput(d,win,ExposureMask|KeyPressMask|StructureNotifyMask);
    Atom close=XInternAtom(d,"WM_DELETE_WINDOW",False);
    XSetWMProtocols(d,win,&close,1);
    GC gc=XCreateGC(d,win,0,NULL);XMapWindow(d,win);
    bool paused=false,done=false,redraw=true;
    struct timespec frame={0,16666667};
    while(!done){
        while(XPending(d)){
            XEvent e;XNextEvent(d,&e);
            if(e.type==ClientMessage&&(Atom)e.xclient.data.l[0]==close)done=true;
            if(e.type==Expose)redraw=true;
            if(e.type==KeyPress){
                KeySym key=XLookupKeysym(&e.xkey,0);
                if(key==XK_q||key==XK_Escape)done=true;
                if(key==XK_space)paused=!paused;
                if(key==XK_n&&paused){update(w);redraw=true;}
            }
        }
        if(!paused &&(!steps||w->tick<steps)){update(w);redraw=true;}
        if(steps&&w->tick>=steps)paused=true;
        if(redraw){draw(w,d,win,gc,screen,paused);redraw=false;}
        nanosleep(&frame,NULL);
    }
    XFreeGC(d,gc);XDestroyWindow(d,win);XCloseDisplay(d);
    return 0;
}
#endif
static void usage(const char *name){
    printf("Usage: %s [--headless] [--steps N] [--objects N] [--seed N] [--prefix PATH] [--periodic] [--no-log]\n",name);
    puts("  --headless      Run without graphics (default if compiled without X11)");
    puts("  --steps N       Simulation steps, 60 per simulated second (default: 36000 headless, unlimited GUI)");
    puts("  --objects N     Number of objects, 2..128 (default: 32)");
    puts("  --seed N        Exact reproducible random seed (default: 20261107)");
    puts("  --prefix PATH   Output log prefix (default: aquarium)");
    puts("  --periodic      Exit one side of space, return at opposite side (default: reflecting walls)");
    puts("  --no-log        Do not generate CSV files");
    puts("  --help          Show help");
}
int main(int argc,char **argv){
    unsigned count=32;uint64_t seed=20261107;const char *prefix="aquarium";
    unsigned long long steps=0;
    bool headless=false,periodic=false,logging=true;
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"--help")){usage(argv[0]);return 0;}
        else if(!strcmp(argv[i],"--headless"))headless=true;
        else if(!strcmp(argv[i],"--periodic"))periodic=true;
        else if(!strcmp(argv[i],"--no-log"))logging=false;
        else if(i+1<argc&&!strcmp(argv[i],"--steps"))steps=strtoull(argv[++i],NULL,10);
        else if(i+1<argc&&!strcmp(argv[i],"--objects"))count=(unsigned)strtoul(argv[++i],NULL,10);
        else if(i+1<argc&&!strcmp(argv[i],"--seed"))seed=strtoull(argv[++i],NULL,10);
        else if(i+1<argc&&!strcmp(argv[i],"--prefix"))prefix=argv[++i];
        else {fprintf(stderr,"Unknown or incomplete argument: %s\n",argv[i]);usage(argv[0]);return 2;}
    }
    if(count<2||count>128){fprintf(stderr,"--objects must be 2..128\n");return 2;}
#ifndef WITH_X11
    headless=true;
#endif
    if(headless&&!steps)steps=36000;
    World w;world_init(&w,count,seed,800,520,periodic);
    if(logging)world_log(&w,prefix);
    printf("DUO Aquarium | %u objects | seed=%" PRIu64 " | %s boundaries | %s\n",count,seed,periodic?"periodic":"reflecting",headless?"headless":"X11 window");
    int rc=0;
    if(headless){for(unsigned long long i=0;i<steps;i++)update(&w);}
#ifdef WITH_X11
    else rc=run_window(&w,steps);
#endif
    results(&w);close_log(&w);free(w.p);
    return rc;
}

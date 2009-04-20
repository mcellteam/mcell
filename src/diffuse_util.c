/**************************************************************************\
 ** File: diffuse_util.c                                                 **
 **                                                                      **
 ** Purpose: Sets up tables of diffusion directions and distances        **
 **                                                                      **
 ** Testing status: essentially unchanged from MCell2.  Compiles.        **
\**************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "logging.h"
#include "diffuse_util.h"
#include "mcell_structs.h"

extern struct volume *world;

#ifndef MY_PI
#define MY_PI 3.14159265358979324
#endif


/***************************************************************************
dgammln:
  In: double > 0.0
  Out: ln( gamma( input ) )
  Note: From Numerical Recipes in C, 2nd ed., p.214, adapted for double
***************************************************************************/

double dgammln(double xx)
{
  double tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			  -1.231739516,0.120858003e-2,-0.536382e-5};
  double stp=2.50662827465;
  int j;

  xx=xx-1.0;
  tmp=xx+5.5;
  tmp=(xx+0.5)*log(tmp)-tmp;
  ser=1.0;
  for (j=0;j<6;j++) {
    ser=ser+cof[j]/++xx;
  }
  xx=tmp+log(stp*ser);
  return(xx);
}


/***************************************************************************
dgser:
  In: double > 0.0
      double >= 0.0
  Out: double containing incomplete gamma of first argument, second
       i.e. int(exp(-t)*t^(arg1-1),t=0..arg2) / gamma(arg1)
  Note: adapted from Numerical Recipes in C, 2nd ed., p.218
        This is the series expansion form.
***************************************************************************/

double dgser(double aa, double xx)
{
  double y,ap,sum,del,eps,gln;
  int itmax,n;

  itmax=100;
  eps=3.0e-7;
  gln=dgammln(aa);
  if (xx<=0.0) {
    if (xx<0.0) mcell_log("DGSER: xx < 0");
    y=0.0;
    return(y);
  }
  ap=aa;
  sum=1.0/aa;
  del=sum;
  for (n=0;n<itmax;n++) {
    del=del*xx/++ap;
    sum=sum+del;
    if (fabs(del)<fabs(sum)*eps) {
      y=sum*exp(-xx+aa*log(xx)-gln);
      return(y);
    }
  }
  mcell_log("DGSER: aa too large, itmax too small");
  y=sum*exp(-xx+aa*log(xx)-gln);
  return(y);
}


/***************************************************************************
dgcf:
  In: double > 0.0
      double >= 0.0
  Out: double containing other incomplete gamma of first argument, second
       i.e. 1 - int(exp(-t)*t^(arg1-1),t=0..arg2) / gamma(arg1)
  Note: adapted from Numerical Recipes in C, 2nd ed., p.219
        This is the continued fraction form.
***************************************************************************/
double dgcf(double aa, double xx)
{
  double y,gold,a0,a1,b0,b1,fac,an,ana,anf,itmax,eps,gln;
  double g = 0.0;

  itmax=100;
  eps=3.0e-7;
  gln=dgammln(aa);
  gold=0.0;
  a0=1.0;
  a1=xx;
  b0=0.0;
  b1=1.0;
  fac=1.0;
  for (an=1.0;an<itmax+1;an++) {
    ana=an-aa;
    a0=(a1+a0*ana)*fac;
    b0=(b1+b0*ana)*fac;
    anf=an*fac;
    a1=xx*a0+anf*a1;
    b1=xx*b0+anf*b1;
    if (a1!=0.0) {
      fac=1./a1;
      g=b1*fac;
      if (fabs((g-gold)/g)<eps) {
	y=g*exp(-xx+aa*log(xx)-gln);
	return(y);
      }
      gold=g;
    }
  }
  mcell_log("DGCF: aa too large, itmax too small");
  y=g*exp(-xx+aa*log(xx)-gln);
  return(y);
}


/***************************************************************************
dgammp:
  In: double > 0.0
      double >= 0.0
  Out: double containing incomplete gamma of first argument, second
       i.e. int(exp(-t)*t^(arg1-1),t=0..arg2) / gamma(arg1)
  Note: adapted from Numerical Recipes in C, 2nd ed., p.218
        Calls dgser or dgcf, depending on which is faster/more accurate
***************************************************************************/

double dgammp(double aa, double xx)
{
  if (xx<aa+1.0) return dgser(aa,xx);
  else return 1.0-dgcf(aa,xx);
}


/***************************************************************************
derf:
  In: double
  Out: double containing error function of input
  Note: adapted from Numerical Recipes in C, 2nd ed., p.220
***************************************************************************/

double derf(double xx)
{
  if (xx<0.0) return -dgammp(0.5,xx*xx);
  else return dgammp(0.5,xx*xx);
}


/***************************************************************************
inverf:
  In: double
  Out: inverse error function of input
  Note: calculated using Euler's method and derf function
***************************************************************************/

#define XPI 1.12837916709551257
/* XPI is 2.0/sqrt(pi) */

double inverf(double xx)
{
  double y,ynew;

  y=-1.0;
  ynew=xx;
  while (fabs(ynew-y)>=1.0e-4) {
    y=ynew;
    ynew=y-(derf(y)-xx)/(XPI*exp(-y*y));
  }
  return(ynew);
}

#undef XPI


/***************************************************************************
r_func:
  In: double containing distance (arbitrary units, mean=1.0)
  Out: double containing probability of diffusing that distance
***************************************************************************/

double r_func(double s)
{   
  double f,s_sqr,val;
    
  f=2.2567583341910251478;  /* 4.0/sqrt(pi) */
  s_sqr=s*s;
  val=f*s_sqr*exp(-s_sqr);
    
  return(val);
}   


/***************************************************************************
init_r_step:
  In: number of desired radial subdivisions
  Out: pointer to array of doubles containing those subdivisions
       returns NULL on malloc failure
  Note: This is for 3D diffusion from a point source (molecule movement)
***************************************************************************/

double* init_r_step(int radial_subdivisions)
{   
  double inc,target,accum,r,r_max,delta_r,delta_r2;
  double *r_step = NULL;
  int j;
    
  r_step = CHECKED_MALLOC_ARRAY_NODIE(double,
                                      radial_subdivisions,
                                      "radial step length table");
  if (r_step == NULL)
    return NULL;
      
  inc=1.0/radial_subdivisions;
  accum=0;
  r_max=3.5;
  delta_r=r_max/(1000*radial_subdivisions);
  delta_r2=0.5*delta_r;
  r=0;
  target=0.5*inc;
  j=0;
  while (j<radial_subdivisions) {
   accum=accum+(delta_r*r_func(r+delta_r2));
   r=r+delta_r;
   if (accum>=target) {
     r_step[j]=r;
     target=target+inc;
     j++;
   }
  }
  no_printf("Min r step = %20.17g   Max r step = %20.17g\n",
            r_step[0],r_step[radial_subdivisions-1]);
  
  return r_step;
}

/***************************************************************************
init_r_step_surface:
  In: number of desired radial subdivisions
  Out: pointer to array of doubles containing those subdivisions
       returns NULL on malloc failure
  Note: This is for 3D molecules emitted from a plane
***************************************************************************/

double* init_r_step_surface(int radial_subdivisions)
{
  double *r_step_s = NULL;
  double p,r_max,r_min,step,r,cdf;
  int i,j;
  static const double sqrt_pi = 1.7724538509055160273;
  
  r_step_s = CHECKED_MALLOC_ARRAY_NODIE(double,
                                        radial_subdivisions,
                                        "radial step length table (surface)");
  if (r_step_s == NULL)
    return NULL;
  
  step = 1.0/radial_subdivisions;
  for ( i=0 , p=(1.0-1e-6)*step ; p<1.0 ; p+=step,i++ )
  {
    r_min = 0;
    r_max = 3.0; /* 17 bit high-end CDF cutoff */
    for (j=0;j<20;j++) /* 20 bits of accuracy */
    {
      r = 0.5*(r_min+r_max);
      cdf = 1.0-exp(-r*r)+sqrt_pi*r*erfc(r);
      if (cdf>p) r_max=r;
      else r_min=r;
    }
    r_step_s[i] = r;
  }
  
  return r_step_s;
}


/***************************************************************************
 init_r_step_3d_release:
  In: number of desired radial subdivisions
  Out: pointer to array of doubles containing those subdivisions
       returns NULL on malloc failure
  Note: This is for 3D molecules emitted from a plane
***************************************************************************/

double* init_r_step_3d_release(int radial_subdivisions)
{
  double *r_step_r = NULL;
  double p,r_max,r_min,step,r,cdf;
  int i,j;
  static const double sqrt_pi_over_2 = 0.886226925452758015;
  
  r_step_r = CHECKED_MALLOC_ARRAY_NODIE(double,
                                        radial_subdivisions,
                                        "radial step length table (3d release)");
  if (r_step_r == NULL)
    return NULL;
  
  step = 1.0/radial_subdivisions;
  for (i=0,p=step*0.5 ; i<radial_subdivisions ; p+=step,i++)
  {
    r_min=0;
    r_max=3.5; /* 18 bit high-end CDF cutoff */
    for (j=0;j<20;j++) /* 20 bits accuracy */
    {
      r = 0.5*(r_min+r_max);
      cdf = 1.0-exp(-r*r)+sqrt_pi_over_2*r*erfc(r);
      if (cdf>p) r_max=r;
      else r_min=r;
    }
    r_step_r[i] = r;
  }
    
  return r_step_r;
}



/***************************************************************************
init_d_step:
  In: number of desired angular subdivisions
      pointer to an int that will contain the actual number of subdivisions
  Out: returns a pointer to array of doubles containing the subdivisions
       returns NULL on malloc failure
  Note: no longer calls init_r_step itself, so you must call both
        data contains three doubles (x,y,z) per requested subdivision
***************************************************************************/

#define DEG_2_RAD  0.01745329251994329576
#define RAD_2_DEG 57.29577951308232087684
/* Multiply by this factor (Pi/180) to convert from degrees to radians */

double* init_d_step(int radial_directions, unsigned int *actual_directions)
{   
  double z;
  double d_phi,phi_mid,phi_edge_prev,phi_edge_approx,phi_factor,theta_mid;
  double *phi_edge = NULL;
  double x_bias,y_bias,z_bias;
  double x,y;
  int i,j,k,n_tot,n_edge,n_patches;
  int *n = NULL;
  double *d_step = NULL;
    
  n_edge=(int) sqrt(radial_directions*MY_PI/2.0);
  n_patches=(int) (2*(n_edge*n_edge)/MY_PI);
  no_printf("desired n_patches in octant = %d\n",radial_directions);
  no_printf("approximate n_patches in octant = %d\n",n_patches);
  phi_edge = CHECKED_MALLOC_ARRAY_NODIE(double, n_edge, "directional step table");
  if (phi_edge == NULL)
    return NULL;

  n=CHECKED_MALLOC_ARRAY_NODIE(int, n_edge, "directional step table");
  if (n == NULL)
  {
    free(phi_edge);
    return NULL;
  }

  for (i=0;i<n_edge;i++) {
    phi_edge[i]=0;
    n[i]=0;
  }

  x_bias=0;
  y_bias=0;
  z_bias=0;
      
  phi_edge_prev=0;
  d_phi=MY_PI/(2.0*n_edge);
  n_tot=0;
  for (i=0;i<n_edge;i++) {
    phi_edge_approx=phi_edge_prev+d_phi; 
    phi_mid=phi_edge_prev+(d_phi/2.0); 
    if (phi_mid<60*DEG_2_RAD) {
      n[i] = (int) ((n_patches*(cos(phi_edge_prev)-cos(phi_edge_approx)))+0.5);
    }
    else {
      n[i] = (int) ((n_patches*(cos(phi_edge_prev)-cos(phi_edge_approx)))-0.5);
    }
    n_tot+=n[i];
    no_printf("i = %d  phi mid = %f  n = %d  n_tot = %d\n",i,phi_mid*RAD_2_DEG,n[i],n_tot);
    phi_edge[i]=acos(cos(phi_edge_prev)-((double)n[i]/n_patches));
    phi_edge_prev=phi_edge[i];
  }
  phi_factor=phi_edge[n_edge-1]*2.0/MY_PI;
  for (i=0;i<n_edge;i++) {
    phi_edge[i]=phi_edge[i]/phi_factor;
  }
  radial_directions=n_tot;
  *actual_directions = n_tot;
  no_printf("actual n_patches in octant = %d\n",radial_directions);
  no_printf("phi factor = %f\n",phi_factor);

  d_step=CHECKED_MALLOC_ARRAY(double, 3*n_tot, "directional step table");
  if (d_step == NULL)
  {
    free(n);
    free(phi_edge);
    return NULL;
  }
  k=0;
  phi_edge_prev=0;
  for (i=0;i<n_edge;i++) {
    phi_mid=(phi_edge_prev+phi_edge[i])/2.0;
    for (j=0;j<n[i];j++) {
      theta_mid=(j*(MY_PI/(2.0*n[i])))+(0.5*MY_PI/(2.0*n[i]));
      x=sin(phi_mid)*cos(theta_mid);
      y=sin(phi_mid)*sin(theta_mid);
      z=cos(phi_mid);
      d_step[k++]=x;
      d_step[k++]=y;
      d_step[k++]=z;
      x_bias+=x;
      y_bias+=y;
      z_bias+=z;
    }
    phi_edge_prev=phi_edge[i];
  }
#ifdef DEBUG
  {
  FILE *fp;
  mcell_log("x_bias = %.17g\n",x_bias);
  mcell_log("y_bias = %.17g\n",y_bias);
  mcell_log("z_bias = %.17g\n",z_bias);
  fp=fopen("vector_table.rib","w");
  for (i=0;i<radial_directions;i++) {
    j=3*i;
/*
    fprintf(fp,"[ %8.5e %8.5e %8.5e ]\n",d_step[j],d_step[j+1],d_step[j+2]);
*/
    fprintf(fp,"AttributeBegin\n");
    fprintf(fp,"Translate %.9g %.9g %.9g\n",
      d_step[j],d_step[j+1],d_step[j+2]);
    fprintf(fp,"ObjectInstance 1\n");
    fprintf(fp,"AttributeEnd\n");
  }
  fclose(fp);
  }
#endif

  free(n);
  free(phi_edge);
  return d_step;
}

#undef DEG_2_RAD
#undef RAD_2_DEG

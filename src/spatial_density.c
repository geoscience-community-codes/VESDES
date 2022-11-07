/*############################################################################
# VESDES (Volcanic Event Spatial Density Estimation Suite)
#
# VESDES includes a number of kernel functions and bandwidth 
# estimation techniques used to estimate spatial density of 
# volcanic events. VESDES is written in the C programming
# language and includes calls to the R statistical package.
#
#    Copyright (C) 2021-2022  
#    Laura Connor (ljconnor@gmail.com)
#    Charles Connor (chuck.connor@gmail.com)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################*/

/* This file is spatial_density.c */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <gc.h>
#define CR 13            /* Decimal code of Carriage Return char */
#define LF 10            /* Decimal code of Line Feed char */

#define PI 3.1415926535897932384626433832795029

/* specifies which kernel method to use for calculating the spatial density */
enum KERNEL {NONE, GAUSS, EPAN, CAUCHY} ; 

typedef struct param {
  int spd;
  int kernel;
  char key[7];
  double west;
  double east;
  double south;
  double north;
  double aoi_northing;
  double aoi_easting;
  char aoi_out[12];
  int utm_zone;
  char *datum;
  char *ellps;
  char *proj;
  double grid_spacing;
  char *event_file;
  int samse;
  double smooth_x;
  double smooth_y;
  double covariance;
  char *bandwidth_file;
  int plot;
  char *plot_dir;
  double map_scale;
  double tick_scale;
  char *map_title;
  char *output_file;
  char *aoi_in; 
  char site[7]; 
} Param;

typedef struct Matrix {
  double h1;
  double h2;
  double h3;
  double h4;
} Matrix;

typedef struct Location {
  double east;
  double north;
} Location;

/*
#############################################
# INPUTS: 
# (1) IN File name of vent locations
# (2) IN/OUT array of vent locations
# OUTPUTS: 
# (1) Number of data lines in the file
############################################# */
Location *load_file(Param *p, int *num_vents, FILE *log) {
  
  FILE *VENTS;
  int maxLineLength = 256;
  char line[256];           /* Line from config file */
  double east;            /* (type=string) Parameter Name  in each line*/
  double north;           /* (type=string) Parameter Value in each line*/
  int i;
  Location *vent;          
  
  /* Open input file of vent locations and create array */
  VENTS = fopen(p->event_file, "r");
  if (VENTS == NULL) {
    fprintf(stderr, "\nERROR : Cannot open %s [%s]!\n", p->event_file, strerror(errno));
    return NULL;
  }
  fprintf(log, "Opening events file %s:\n", p->event_file);
  fflush(log);
  
  *num_vents = 0;
  while (fgets(line, maxLineLength, VENTS) != NULL) {
    if (line[0] == '#' || line[0] == '\n' || line[0] == ' ') continue;
    if (strlen(line) > 2) (*num_vents)++;
  }
  rewind(VENTS);
  
  vent = (Location *) GC_MALLOC ( (size_t)*num_vents * sizeof(Location));
  if (vent == NULL) {
    fprintf(stderr, "Unable to allocate memory for vent array! Exiting.\n");
    return NULL;
  }
  
  for (i = 0; i < *num_vents; i++) {
    if (fgets(line, maxLineLength, VENTS) == NULL) break;
		
    /*if first character is comment, new line, space, return to next line*/
    if (line[0] == '#' || line[0] == '\n' || line[0] == ' ') continue;
  
    sscanf (line,"%lf %lf", &east, &north); /*split line into easting and northing*/
    /* Convert vent location to km */
    vent[i].east = east/1000.0;
    vent[i].north = north/1000.0;
    fprintf (log, "%lf %lf\n", vent[i].east, vent[i].north); 
  }  
  
  fclose (VENTS);
  fprintf (stderr, "Loaded %d vents from %s\n", *num_vents, p->event_file);
  return vent;
}

double determinant(Matrix *m) {
  return ((m->h1 * m->h4) - (m->h2 * m->h3)); 
}

int sqrtM(Matrix *m, Matrix *sqrtm, double det) {
  
  double trace, t, s, t_inv;
  
  trace = m->h1 + m->h4;
  s = sqrt(det);
  t = sqrt(trace + 2*s);
  t_inv = 1.0/t;
  
  if (!t) {
    return 1;
  }
  else {
    
    sqrtm->h1 = t_inv * (m->h1 + s);
    sqrtm->h2 = t_inv * m->h2;
    sqrtm->h3 = t_inv * m->h3;
    sqrtm->h4 = t_inv * (m->h4 + s);
  }
  return 0; 
}

int inv_matrix(Matrix *sqrtm, Matrix *sqrtmi) {
  
  double det, inv_det;
  det = determinant(sqrtm);
  inv_det = 1.0 / det;
  
  if (!inv_det) {
    return 1;
  }
  else {
    sqrtmi->h1 = inv_det * sqrtm->h4;
    sqrtmi->h2 = inv_det * sqrtm->h2 * -1.0;
    sqrtmi->h3 = inv_det * sqrtm->h3 * -1.0;
    sqrtmi->h4 = inv_det * sqrtm->h1;
  }
  
  return 0;
}

/*##################################################################
# Function gauss()
# INPUTS:
# (1) X (km) current grid location
# (2) Y (km) current grid location
# (3) reference to array of volcanic vent locations (km)
# (4) number of vent locations
# (5) Matrix: this is the inverse of the square root
#     of the (2 x 2)bandwidth matrix
# (6) Const: this is 2*pi*determinant(H)
#     note: spd = 1 for spatial intensity
# OUTPUTS:
# (1) lambda (i.e. spatial intensity/density at the current grid location)
####################################################################*/
double gauss(double x1, double y1, Location vent[], int num_vents, double Const, Matrix *ma) {
  
  int i;
  double dx, dy, lambda, dist, edist;
  Matrix dxdy, Mdxdy;
  double sum = 0.0;
  
  for (i = 0; i < num_vents; i++) { 
    
  /* For each event 
     Get distance (km) from event to grid point
  */
      dx = x1 - (vent+i)->east;
      dy = y1 - (vent+i)->north;
       
      dxdy.h1 = dx;
      dxdy.h2 = dy;
      
      Mdxdy.h1 = ma->h1 * dxdy.h1 + ma->h2 * dxdy.h2;
      Mdxdy.h2 = ma->h3 * dxdy.h1 + ma->h4 * dxdy.h2;
      
      dist = (Mdxdy.h1 * Mdxdy.h1) + (Mdxdy.h2 * Mdxdy.h2);
      
      dist *= -0.5;
      edist = exp(dist);
      sum += edist;
  }
  lambda = sum/Const;
  return lambda;
}

double gaussL(double x1, double y1, Location vent[], int num_vents, double Const, Matrix *ma) {

  int i;
  double dx, dy, lambda, dist, edist;
  Matrix dxdy, Mdxdy;
  double sum = 0.0;

  for (i = 0; i < num_vents; i++) {
    if (x1 != (vent+i)->east && y1 != (vent+i)->north) {
    /* For each event
       Get distance (km) from current event to all other events
     */
      dx = x1 - (vent+i)->east;
      dy = y1 - (vent+i)->north;

      dxdy.h1 = dx;
      dxdy.h2 = dy;

      Mdxdy.h1 = ma->h1 * dxdy.h1 + ma->h2 * dxdy.h2;
      Mdxdy.h2 = ma->h3 * dxdy.h1 + ma->h4 * dxdy.h2;

      dist = (Mdxdy.h1 * Mdxdy.h1) + (Mdxdy.h2 * Mdxdy.h2);

      dist *= -0.5;
      edist = exp(dist);
      sum += edist;
    }
  }
  lambda = sum/Const;
  return lambda;
}

/*##################################################################
# Function epan()
# INPUTS: 
# (1) X (km) current grid location
# (2) Y (km) current grid location
# (3) reference to array of volcanic vent locations (km)
# (4) number of vent locations
# (5) Matrix: this is the inverse of the square root 
#     of the (2 x 2)bandwidth matrix
# (6) Const: PI * sqrt_detH * (double)spd;
#     note: spd = 1 for spatial intensity
# OUTPUTS:
# (1) lambda (i.e. spatial intensity/density at the current grid location)
####################################################################*/
double epan(double x1, double y1, Location vent[], int num_vents, double Const, Matrix *ma) {
  int i;
  double dx, dy, lambda, dist;
  Matrix dxdy, Mdxdy;
  double sum = 0.0;
    
  /* For each event 
     Get distance (km) from event to grid point
   */
  for (i = 0; i < num_vents; i++) {
    dx = x1 - (vent+i)->east;
    dy = y1 - (vent+i)->north;
       
    dxdy.h1 = dx;
    dxdy.h2 = dy;
    Mdxdy.h1 = ma->h1 * dxdy.h1 + ma->h2 * dxdy.h2;
    Mdxdy.h2 = ma->h3 * dxdy.h1 + ma->h4 * dxdy.h2;
    dist = (Mdxdy.h1 * Mdxdy.h1) + (Mdxdy.h2 * Mdxdy.h2);
    
    if ( dist < 1 ) sum += (1.0 - dist); 
  }
  lambda = Const * sum;
  return lambda;
  
}

double epanL(double x1, double y1, Location vent[], int num_vents, double Const, Matrix *ma) {
  int i;
  double dx, dy, lambda, dist;
  Matrix dxdy, Mdxdy;
  double sum = 0.0;

  /* For each event
     Get distance (km) from event to current event
   */
  for (i = 0; i < num_vents; i++) {
    if (x1 != (vent+i)->east && y1 != (vent+i)->north) {
      dx = x1 - (vent+i)->east;
      dy = y1 - (vent+i)->north;

      dxdy.h1 = dx;
      dxdy.h2 = dy;
      Mdxdy.h1 = ma->h1 * dxdy.h1 + ma->h2 * dxdy.h2;
      Mdxdy.h2 = ma->h3 * dxdy.h1 + ma->h4 * dxdy.h2;
      dist = (Mdxdy.h1 * Mdxdy.h1) + (Mdxdy.h2 * Mdxdy.h2);

      if ( dist < 1 ) sum += (1.0 - dist);

    }
  }
  lambda = Const * sum;
  return lambda;

}

/*##################################################################
# Function cauchy()
# INPUTS: 
# (1) X (km) current grid location
# (2) Y (km) current grid location
# (3) reference to array of volcanic vent locations (km)
# (4) number of vent locations
# (5) Matrix: this is the inverse of the square root 
#     of the (2 x 2)bandwidth matrix
# (6) Const: 2 * PI * (double)spd
#     note: spd = 1 for spatial intensity
# OUTPUTS:
# (1) lambda (i.e. spatial intensity/density at the current grid location)
####################################################################*/

double cauchy(double x1, double y1, Location vent[], int num_vents, double Const, Matrix *ma) {
  
  int i;
  double dx, dy, lambda, dist;
  double sum = 0.0;
  Matrix dxdy, Mdxdy;
  double K_t;
  
 for (i = 0; i < num_vents; i++) {
    dx = x1 - (vent+i)->east;
    dy = y1 - (vent+i)->north;
       
    dxdy.h1 = dx;
    dxdy.h2 = dy;
    Mdxdy.h1 = ma->h1 * dxdy.h1 + ma->h2 * dxdy.h2;
    Mdxdy.h2 = ma->h3 * dxdy.h1 + ma->h4 * dxdy.h2;
    dist = 1.0 + (Mdxdy.h1 * Mdxdy.h1) + (Mdxdy.h2 * Mdxdy.h2);
    
    K_t = 1.0 / pow(dist, 1.5);
    sum += K_t;    
  }
  lambda = sum * Const;
  return lambda;
}

double cauchyL(double x1, double y1, Location vent[], int num_vents, double Const, Matrix *ma) {

  int i;
  double dx, dy, lambda, dist;
  double sum = 0.0;
  Matrix dxdy, Mdxdy;
  double K_t;

 for (i = 0; i < num_vents; i++) {
   if (x1 != (vent+i)->east && y1 != (vent+i)->north) {
     //fprintf(stderr, "%f=%f, %f=%f\n", x1, (vent+i)->east, y1, (vent+i)->north);
     dx = x1 - (vent+i)->east;
     dy = y1 - (vent+i)->north;

     dxdy.h1 = dx;
     dxdy.h2 = dy;
     Mdxdy.h1 = ma->h1 * dxdy.h1 + ma->h2 * dxdy.h2;
     Mdxdy.h2 = ma->h3 * dxdy.h1 + ma->h4 * dxdy.h2;
     dist = 1.0 + (Mdxdy.h1 * Mdxdy.h1) + (Mdxdy.h2 * Mdxdy.h2);

     K_t = 1.0 / pow(dist, 1.5);
     sum += K_t;
   }
 }
  lambda = sum * Const;
  return lambda;

}

int main(int argc, char *argv[]) {
  
  FILE *log;
  FILE *CONF = NULL;
  FILE *KOP = NULL;
  FILE *BW = NULL;
  FILE *OUT = NULL;
  FILE *AOI_SPD = NULL;
  FILE *AOI_IN = NULL;
  int maxLineLength = 256;
  char line[256];             /* Line from config file */
  char var[64];               /* (type=string) Parameter Name  in each line*/
  char value[256];            /* (type=string) Parameter Value in each line*/
  Param P;
  Matrix H, sqrtH, sqrtHi;
  Location *vent = NULL;
  double (*kernel_function)(double, double, Location vent[], int, double, Matrix *);
  int i;
  int spd;
  int num_vents;
  double detH, sqrt_detH, Const = 0, grid2, site_const = 0;
  double grid_total = 0.0, 
         pdf = 0.0, sum_logL = 0.0, log_pdf = 0.0,
         X_easting, Y_northing, XX, YY, 
	 x_smoothing, y_smoothing, covariance;

  if (argc < 1) {
    fprintf(stderr, "USAGE: %s <file.conf>\n\n", argv[0]);
    return 1;
  }

  fprintf(stderr, "Opening and appending run info to to: logfile\n");
  log = fopen("logfile", "w");
  if (log == NULL) {
    fprintf(stderr, "\nERROR : Cannot open logfile [%s]!\n", strerror(errno));
    return 1;
  } 
  
  fprintf(log, "Parameters:\n");
  fflush(log); 

  CONF = fopen(argv[1], "r"); /*open configuration file*/
  if (CONF == NULL) {
    fprintf(stderr, "\nERROR : Cannot open %s [%s]!\n", argv[1], strerror(errno));
    fclose (log);
    return 1;
  }

  /* Initialize */
  kernel_function = &gauss;
  
  P.spd = 0;
  P.kernel = NONE;
  P.west = (double)0;
  P.east = (double)0;
  P.south = (double)0;
  P.north = (double)0;
  P.aoi_northing = (double)0;
  P.aoi_easting = (double)0;
  strcpy( P.aoi_out, "aoi_spd.txt");
  P.grid_spacing = (double)0;
  P.event_file = NULL;
  P.samse = 0;
  P.smooth_x = (double)0;
  P.smooth_y = (double)0;
  P.covariance = (double)0;
  P.bandwidth_file = NULL;
  P.output_file = NULL;
  strcpy(P.site, "");
  

  while (fgets(line, maxLineLength, CONF) != NULL) {
		
    /*if first character is comment, new line, space, return to next line*/
    if (line[0] == '#' || line[0] == '\n' || line[0] == ' ') continue;
    
    /*print incoming parameter*/
    var[0] = '\0';
    value[0] = '\0';
    sscanf (line,"%s =%s", var, value); /*split line into before ' = ' and after*/
    fprintf(stdout, "%s = %s\n",var, value); /*print incoming parameter value*/
    fflush(stdout); 
    
    if (!strncmp(var, "SPD", strlen("SPD"))) { 
      sscanf(value, "%d", &P.spd);
      if (P.spd != 1 && P.spd != 2) { 
        fprintf(stderr, "Parameter SPD[%d] must = one of these options: 1 (calculate spatial density) or 2 (calculate spatial intensity)\n", P.spd);
        fclose(CONF);
        fclose(log);
        return 1;
      }
      fprintf(log, "%s = %d\n", var, P.spd); /*print incoming parameter value*/
    }
    else if (!strncmp(var, "KERNEL", strlen("KERNEL"))) {
      sscanf(value, "%d", &P.kernel);
      if (P.kernel != 1 && P.kernel != 2 && P.kernel != 3 && P.kernel != 4) {
        fprintf(stderr, 
        "Parameter KERNEL[%d] must = one of these options: 1 (use gaussian kernel function) or 2 (use epanechnikov kernel function) or 3 (use cauchy kernel function) \n", P.kernel);
        fclose(CONF);
        fclose(log);
        return 1;
      }
      else {
        
        if (P.kernel == 1) {
          kernel_function = &gauss;
          sprintf (P.key, "gauss");
        }
        else if (P.kernel == 2) {
          kernel_function = &epan;
          sprintf (P.key, "epan");
        }
        else if (P.kernel == 3) {
          kernel_function = &cauchy;
          sprintf (P.key, "cauchy");
        }
        fprintf(log, "%s = %s(%d)\n", var, P.key, P.kernel); /*print incoming parameter value*/
      }
    }
    else if  (!strncmp(var, "WEST", strlen("WEST"))) {
      sscanf(value, "%lf", &P.west);
      /* Convert map coords to km */
      P.west /= 1000.0;
      fprintf(log, "%s = %lf (km)\n",var, P.west); /*print incoming parameter value*/
    }
    else if  (!strncmp(var, "EAST", strlen("EAST"))) {
      sscanf(value, "%lf", &P.east);
      /* Convert map coords to km */
      P.east /= 1000.0;
      fprintf(log, "%s = %lf (km)\n",var, P.east); /*print incoming parameter value*/
    }
    else if  (!strncmp(var, "SOUTH", strlen("SOUTH"))) {
      sscanf(value, "%lf", &P.south);
      /* Convert map coords to km */
      P.south /= 1000.0;
      fprintf(log, "%s = %lf (km)\n",var, P.south); /*print incoming parameter value*/
    }
    else if  (!strncmp(var, "NORTH", strlen("NORTH"))) {
      sscanf(value, "%lf", &P.north);
      /* Convert map coords to km */
      P.north /= 1000.0;
      fprintf(log, "%s = %lf (km)\n",var, P.north); /*print incoming parameter value*/
    }
    else if (!strncmp(var, "SAMSE", strlen("SAMSE"))) { 
      sscanf(value, "%d", &P.samse);
      if (P.samse != 1 && P.samse != 0) { 
        fprintf(stderr, "Parameter SAMSE[%d] must = one of these options: 1 (use SAMSE method) or 0 (input X- and Y- smoothing and covartiance directly)\n", P.samse);
        fclose(CONF);
        fclose(log);
        return 1;
      }
      fprintf(log, "%s = %d\n", var, P.samse); /*print incoming parameter value*/
    }
    else if  (!strncmp(var, "SMOOTH_X", strlen("SMOOTH_X"))) {
      sscanf(value, "%lf", &P.smooth_x);
      /* Convert smoothing to km */
      P.smooth_x /= 1000.0;
      fprintf(log, "%s = %lf (km)\n",var, P.smooth_x); /*print incoming parameter value*/
    }
    else if  (!strncmp(var, "SMOOTH_Y", strlen("SMOOTH_Y"))) {
      sscanf(value, "%lf", &P.smooth_y);
      /* Convert smoothing to km */
      P.smooth_y /= 1000.0;
      fprintf(log, "%s = %lf (km)\n",var, P.smooth_y); /*print incoming parameter value*/
    }
    else if  (!strncmp(var, "AOI_NORTHING", strlen("AOI_NORTHING"))) {
      sscanf(value, "%lf", &P.aoi_northing);
      /* Convert map coords to km */
      P.aoi_northing /= 1000.0;
      fprintf(log, "%s = %lf (km)\n",var, P.aoi_northing); /*print incoming parameter value*/
    }
    else if  (!strncmp(var, "AOI_EASTING", strlen("AOI_EASTING"))) {
      sscanf(value, "%lf", &P.aoi_easting);
      /* Convert map coords to km */
      P.aoi_easting /= 1000.0;
      fprintf(log, "%s = %lf (km)\n",var, P.aoi_easting); /*print incoming parameter value*/
    }
    else if  (!strncmp(var, "COVARIANCE", strlen("COVARIANCE"))) {
      sscanf(value, "%lf", &P.covariance);
      fprintf(log, "%s = %lf\n",var, P.covariance); /*print incoming parameter value*/
    }
    else if  (!strncmp(var, "GRID_SPACING", strlen("GRID_SPACING"))) {
      sscanf(value, "%lf", &P.grid_spacing);
      /* Convert grid spacing to km */
      P.grid_spacing /= 1000.0;
      fprintf(log, "%s = %lf (km)\n",var, P.grid_spacing); /*print incoming parameter value*/
    }
    if (!strncmp(var, "EVENT_FILE", strlen("EVENT_FILE"))) { 
      P.event_file = (char*) GC_MALLOC(sizeof(char) * (strlen(value)+1));
      if (P.event_file == NULL) {
        fprintf(stderr, "\n[INITIALIZE] Out of Memory assigning EVENT_FILE!\n");
        fclose(CONF);
        fclose(log);
        return 1;
      }
      sscanf(value, "%s", P.event_file);
      fprintf(log, "%s = %s\n", var, P.event_file); /*print incoming parameter value*/
    }
    if (!strncmp(var, "BANDWIDTH_FILE", strlen("BANDWIDTH_FILE"))) { 
      P.bandwidth_file = (char*) GC_MALLOC(sizeof(char) * (strlen(value)+1));
      if (P.bandwidth_file == NULL) {
        fprintf(stderr, "\n[INITIALIZE] Out of Memory assigning BANDWIDTH_FILE!\n");
        fclose(CONF);
        fclose(log);
        return 1;
      }
      sscanf(value, "%s", P.bandwidth_file);
      fprintf(log, "%s = %s\n", var, P.bandwidth_file); /*print incoming parameter value*/
    }
    if (!strncmp(var, "AOI_FILE", strlen("AOI_FILE"))) {
     P.aoi_in = (char*) GC_MALLOC(sizeof(char) * (strlen(value)+1));
     if (P.aoi_in == NULL) {
        fprintf(stderr, "\n[INITIALIZE] Out of Memory assigning AOI_IN!\n");
        fclose(CONF);
        fclose(log);
        return 1;
      }
      sscanf(value, "%s", P.aoi_in);
      fprintf(log, "%s = %s\n", var, P.aoi_in); /*print incoming parameter value*/ 
    }
  }
  fflush(stdout);
  fclose(CONF); 
  
  if (P.samse > 0) {
    /* FIND BANDWIDTH using SAMSE bandwidth from R */
    fprintf (log, "\nOptimizing Pilot Bandwidth (SAMSE)\n");
    system ("touch bandwidth.dat");
    KOP = fopen("R-samse", "w");
    if (KOP == NULL) {
      fprintf (stderr, "\nERROR : Cannot create R-script file [%s]!\n", strerror(errno));
      fclose (log);
      return 1;
    }
    fprintf (log, "Creating R script file: R-samse\n");

    fprintf (KOP, "library(ks)\n");
    fprintf (KOP, "vents<-read.table(file=\"%s\") [1:2]\n", P.event_file); /* locations are UTM meters */
    fprintf (KOP, "bd <- Hpi(x=vents,nstage=2,pilot=\"samse\",pre=\"sphere\", binned=FALSE, amise=FALSE, deriv.order=0, verbose=FALSE,optim.fun=\"nlm\")\n"); /* command to run samse */
    fprintf (KOP, "sink(\"%s\")\n", P.bandwidth_file); /* designates write-to file */
    fprintf (KOP, "show(bd)\n"); 	/* should be 2x2 matrix */
    fprintf (KOP, "sink()\n"); /* clears sink */
    fclose (KOP);

    system("R CMD BATCH R-samse");
    BW = fopen (P.bandwidth_file, "r");
    if (BW == NULL) {
      fprintf (stderr, "\nERROR : Cannot open %s file [%s]!\n", P.bandwidth_file, strerror(errno));
      fprintf (stderr, "Check bandwidth parameters in configuration file (spatial_density.conf)\n");
      fclose (log);
      return 1;
    }

    i = 0;
    H.h1 = H.h2 = H.h3 = H.h4 = 0;
    
    fprintf(log, "'R' bandwidth matrix in sq m\n");
    while (fgets (line, maxLineLength, BW) != NULL) {
      if (i > 2) break;
      if (!i) fprintf(stdout, "%s", line);
      if (i == 1) {
        fprintf (log, "%s", line);
        sscanf (line, "[1,] %lf %lf", &H.h1, &H.h2);
      }
      if (i == 2) {
        fprintf (log, "%s", line);
        sscanf (line, "[2,] %lf %lf", &H.h3, &H.h4);
      }
      i++;
    }
    fclose (BW);

    /* Convert bandwidth matrix from 'R' to KM */
    H.h1 /= 1e6;
    H.h2 /= 1e6;
    H.h3 /= 1e6;
    H.h4 /= 1e6;
  }
  else if (P.smooth_x > 0 && P.smooth_y > 0) {
    /* FIND BANDWIDTH using user defined parameters */
    fprintf (log, "\nUsing user defined smoothing Bandwidth\n");
    x_smoothing = P.smooth_x;
    y_smoothing = P.smooth_y;
    covariance= P.covariance;
    
    BW = fopen (P.bandwidth_file, "w");
    if (BW == NULL) {
      fprintf (stderr, "\nERROR : Cannot open %s file [%s]!\n", P.bandwidth_file, strerror(errno));
      fclose (log);
      return 1;
    }
    fprintf (BW, "Using user defined bandwidth\n");
    fclose(BW);

    H.h1 = x_smoothing * x_smoothing;
    H.h2 = H.h3 = covariance * x_smoothing * y_smoothing;
    H.h4 = y_smoothing * y_smoothing;
  }
   
  fprintf(log, "Bandwidth in sq km\n");
  fprintf (log, "%lf %lf\n", H.h1, H.h2);
  fprintf (log, "%lf %lf\n", H.h3, H.h4);

  /* Create output file */
  P.output_file = (char*) GC_MALLOC(sizeof(char) * (strlen(P.event_file) + strlen(".samse.xyz") + 1));
  if (P.output_file == NULL) {
    fprintf(stderr, "\n[INITIALIZE] Out of Memory assigning output filename!\n");
    fclose(log);
    return 1;
  }

  sprintf (P.output_file, "%s.%s.xyz", P.event_file, P.key);
  OUT = fopen (P.output_file, "w");
  if (OUT == NULL) {
    fprintf (stderr, "\nERROR : Cannot open output file %s [%s]!\n", P.output_file, strerror(errno));
    fclose (log);
    return 1;
  }
  fprintf(log, "Opening %s for output\n", P.output_file);
  fflush(log);
  
  /* Load vent locations (meters) */
  vent = load_file(&P, &num_vents, log);
  if (vent == NULL) {
    fprintf (stderr, "[ERROR]No vents read! Exiting.");
    fclose(log);
    return 1;
  }

  /* Calculate spatial density or intensity */
  spd = num_vents;

  if (P.spd == 2) { /* Calculate spatial intensity */
    fprintf (log, "Calculating spatial intensity.\n");
    spd = 1; 
  } else {
    fprintf (log, "Calculating spatial density; grid should sum to 1.\n");
    fprintf (log, "Number of vents = %d\n", spd);
  }
  
  /* Calculate necessary constants for 
     the bivariate anisotropic  kernel fuctions:
   */

  /* determinant of the bandwidth matrix in km*/
  detH = determinant(&H);
  fprintf (log, "Determinant: %lf (km)\n", detH);
  
  /* square root of the bandwidth matrix in km*/
  if (sqrtM(&H, &sqrtH, detH)) {
    fprintf(stderr, "ERROR calculating square root matrix\n");
    fclose(log);
    return 0;
  }
  fprintf (log, "Square Root Matrix (km):\n");
  fprintf (log, "%lf %lf\n", sqrtH.h1, sqrtH.h2);
  fprintf (log, "%lf %lf\n", sqrtH.h3, sqrtH.h4);
  
  /* square root of the determinant in km*/
  sqrt_detH = sqrt(detH);
  fprintf (log, "sqrt(Determinant (km)): %lf\n", sqrt_detH);
  
  /* determinant of sqrtH 
  det_sqrtH = determinant(&sqrtH);
  fprintf (log, "Determinant(sqrtH): %lf\n", det_sqrtH);
  fprintf (stderr, "Determinant(sqrtH): %lf\n", det_sqrtH);
  */
  
  /* inverse of the square root matrix */
  if (inv_matrix(&sqrtH, &sqrtHi)) {
    fprintf(stderr, "ERROR calculating inverse of square root matrix\n");
    fclose(log);
    return 0;
  }
  fprintf (log, "Inverse Square Root Matrix (km):\n");
  fprintf (log, "%lf %lf\n", sqrtHi.h1, sqrtHi.h2);
  fprintf (log, "%lf %lf\n", sqrtHi.h3, sqrtHi.h4);
  

/* 
 This is to calculate spatial density / spatial intensity
 that is derived from the number of vents. Do the
 summation then multiply by this constant. 
*/
  if (P.kernel == 1) {
  /* gaussian constant */
    Const = 2.0 * PI * sqrt_detH * (double)spd;
    site_const = 2.0 * PI * sqrt_detH * (double)num_vents;
  }
  else if (P.kernel == 2) {
  /* Epanechnikov constant */
    Const = 2.0 / (PI * sqrt_detH * (double)spd);
    site_const = 2.0 / (PI * sqrt_detH * (double)num_vents);
  }
  else if (P.kernel == 3) {
  /* Cauchy constant */
    Const = 1.0 / (2.0 * PI * sqrt_detH * (double)spd);
    site_const = 1.0 / (2.0 * PI * sqrt_detH * (double)num_vents);
  }
  else 
    fprintf (stderr, "valid kernel function[%d] not selected\n", P.kernel);
 
  fprintf (log, "Const: %lf\n", Const);
  
  /* Create the spatial intensity grid in km*/
  grid_total = 0.0;
  grid2 = P.grid_spacing * P.grid_spacing;
  
  /* First calculate the spatial density at the sites */
  AOI_SPD = fopen (P.aoi_out, "w");
  if (AOI_SPD == NULL) {
     fprintf (stderr, "\nERROR : Cannot open output file %s [%s]!\n", P.aoi_out, strerror(errno));
  } else {
    fprintf(log, "Opening %s for output\n", P.aoi_out);
    fflush(log);
  }
  /* if site file is specified
   340972 4828222 7,Helvetica-Narrow-Bold,0 0 TR ATRC
   */
  AOI_IN = fopen(P.aoi_in, "r");
  if (AOI_IN == NULL) {
    fprintf(stderr, "Cannot open Sites file %s [%s]!\n", P.aoi_in, strerror(errno));
    /* Othrwise, check for a single site */
    if (P.aoi_easting && P.aoi_northing) {
      pdf = (*kernel_function)(P.aoi_easting, P.aoi_northing, vent, num_vents, site_const, &sqrtHi);
      XX = P.aoi_easting * 1000.0;
      YY = P.aoi_northing * 1000.0;        
      /* Print out locations in meters; pdf in sq km); */
      fprintf (AOI_SPD, "%lf %lf %g\n", XX, YY, pdf);
    }
  } else {
    while (fgets (line, maxLineLength, AOI_IN) != NULL) {
      fprintf (log, "%s", line);
      sscanf (line, "%lf %lf %*s %*d %*s %s", &P.aoi_easting, &P.aoi_northing, P.site);
      P.aoi_easting /= 1000.0;
      P.aoi_northing /= 1000.0;
      pdf = (*kernel_function)(P.aoi_easting, P.aoi_northing, vent, num_vents, site_const, &sqrtHi);
      XX = P.aoi_easting * 1000.0;
      YY = P.aoi_northing * 1000.0;
        
      /* Print out locations in meters; pdf in sq km); */
      fprintf (AOI_SPD, "%lf %lf %g %s\n", XX, YY, pdf, P.site);
    }
    fclose(AOI_IN);
  }
  fclose(AOI_SPD);
  
  /* Create grid of spatial densities */
  X_easting = P.west - P.grid_spacing; 
  do {
    X_easting += P.grid_spacing;
    Y_northing = P.south - P.grid_spacing;
    do {
      Y_northing += P.grid_spacing;
      pdf = (*kernel_function)(X_easting, Y_northing, vent, num_vents, Const, &sqrtHi);
      XX = X_easting * 1000.0;
      YY = Y_northing * 1000.0;
      pdf *= grid2;
      grid_total += pdf;

      /* Print out locations in meters; pdf in sq km); */
      fprintf (OUT, "%lf %lf %g\n", XX, YY, pdf);
    }while (Y_northing < P.north);
   
  }while (X_easting < P.east);

  fprintf (log, "Grid totals %g Finished Calculations.\n", grid_total); 

   /*
  For each vent location, calculate the spatial density
  at this location, due to the remaining vent locations.
  Calculate the log of this value.
  Sum the logs.
  */

  /* Calculate spatial density */
  spd = num_vents - 1;

  /*
   This is to calculate spatial density
   that is derived from the number of vents. Do the
   summation then multiply by this constant.
   */
  
  if (P.kernel == GAUSS) {
    kernel_function = &gaussL;
  /* gaussian constant */
    Const = 2.0 * PI * sqrt_detH * (double)spd;
  }
  else if (P.kernel == EPAN) {
    kernel_function = &epanL;
  /* Epanechnikov constant */
    Const = 2.0 / (PI * sqrt_detH * (double)spd);
  }
  else if (P.kernel == CAUCHY) {
    kernel_function = &cauchyL;
  /* Cauchy constant */
    Const = 1.0 / (2.0 * PI * sqrt_detH * (double)spd);
  }
  else {
    fprintf (stderr, "valid kernel function[%d] not selected\n", P.kernel);
    return 1;
  }


  for ( i = 0; i < num_vents; i++ ) {
      log_pdf = 0.0;
      X_easting = (vent+i)->east; //vent[i].east;
      Y_northing = (vent+i)->north; //vent[i].north;
      pdf = (*kernel_function)(X_easting, Y_northing, vent, num_vents, Const, &sqrtHi);
      //XX = X_easting * 1000.0;
      //YY = Y_northing * 1000.0;
      if (pdf <= 0.0) log_pdf = -100;
      else log_pdf = log10(pdf);
      sum_logL += log_pdf;
      //fprintf(log, "[%d] %0.f %0.f %g %g\n", i, XX, YY, pdf, log_pdf);
  }
  sum_logL = abs(sum_logL);
  fprintf(log, "Log-likelihood = %g",sum_logL);

  fclose(log);

  fprintf (stdout, "Grid totals %g Finished Calculations.\n", grid_total);
  fprintf(stdout, "Log-likelihood (%s) = %g\n", P.key, sum_logL);
  return 1;
}



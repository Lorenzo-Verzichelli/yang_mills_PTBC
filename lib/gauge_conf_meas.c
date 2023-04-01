#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include"../include/macro.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/gparam.h"
#include"../include/function_pointers.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"
#include"../include/geometry.h"


#include<time.h> // DEBUG

// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double plaquettep(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  long r,
                  int i,
                  int j)
   {
   GAUGE_GROUP matrix;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param; // just to avoid warning at compile time
   #endif

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

   equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(&matrix, &(GC->lattice[r][i]));
   times_equal(&matrix, &(GC->lattice[r][j]));

   return retr(&matrix);
   }


// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double complex plaquettep_complex(Gauge_Conf const * const GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  long r,
                                  int i,
                                  int j)
   {
   GAUGE_GROUP matrix;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param; // just to avoid warning at compile time
   #endif

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

   equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(&matrix, &(GC->lattice[r][i]));
   times_equal(&matrix, &(GC->lattice[r][j]));

   return retr(&matrix)+I*imtr(&matrix);
   }


// computation of the plaquette (matrix) in position r and positive directions i,j
void plaquettep_matrix(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r,
                       int i,
                       int j,
                       GAUGE_GROUP *matrix)
   {
   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param; // just to avoid warning at compile time
   #endif

//
//       ^ j
//       |   (3)
//       +---<---+
//       |       |
//   (4) V       ^ (2)
//       |       |
//       +--->---+---> i
//       r   (1)
//

   equal(matrix, &(GC->lattice[r][i]));
   times_equal(matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(matrix, &(GC->lattice[r][j]));
   }


// compute the four-leaf clover in position r, in the plane i,j and save it in M
void clover(Gauge_Conf const * const GC,
            Geometry const * const geo,
            GParam const * const param,
            long r,
            int i,
            int j,
            GAUGE_GROUP *M)
   {
   GAUGE_GROUP aux;
   long k, p;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(i >= STDIM || j >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param;
   #endif

   zero(M);

//
//                   i ^
//                     |
//             (14)    |     (3)
//         +-----<-----++-----<-----+
//         |           ||           |
//         |           ||           |
//   (15)  V      (13) ^V (4)       ^ (2)
//         |           ||           |
//         |   (16)    || r   (1)   |
//    p    +----->-----++----->-----+------>   j
//         +-----<-----++-----<-----+
//         |    (9)    ||   (8)     |
//         |           ||           |
//    (10) V      (12) ^V (5)       ^ (7)
//         |           ||           |
//         |           ||           |
//         +------>----++----->-----+
//              (11)   k      (6)
//
   // avanti-avanti
   equal(&aux, &(GC->lattice[r][i]) );                           // 1
   times_equal(&aux, &(GC->lattice[nnp(geo, r, i)][j]) );        // 2
   times_equal_dag(&aux, &(GC->lattice[nnp(geo, r, j)][i]) );    // 3
   times_equal_dag(&aux, &(GC->lattice[r][j]) );                 // 4
   plus_equal(M, &aux);

   k=nnm(geo, r, j);

   // avanti-indietro
   equal_dag(&aux, &(GC->lattice[k][j]) );                       // 5
   times_equal(&aux, &(GC->lattice[k][i]) );                     // 6
   times_equal(&aux, &(GC->lattice[nnp(geo, k, i)][j]) );        // 7
   times_equal_dag(&aux, &(GC->lattice[r][i]) );                 // 8
   plus_equal(M, &aux);

   p=nnm(geo, r, i);

   // indietro-indietro
   equal_dag(&aux, &(GC->lattice[p][i]) );                       // 9
   times_equal_dag(&aux, &(GC->lattice[nnm(geo, k, i)][j]) );    // 10
   times_equal(&aux, &(GC->lattice[nnm(geo, k, i)][i]) );        // 11
   times_equal(&aux, &(GC->lattice[k][j]) );                     // 12
   plus_equal(M, &aux);

   // indietro-avanti
   equal(&aux, &(GC->lattice[r][j]) );                            // 13
   times_equal_dag(&aux, &(GC->lattice[nnp(geo, p, j)][i]) );     // 14
   times_equal_dag(&aux, &(GC->lattice[p][j]) );                  // 15
   times_equal(&aux, &(GC->lattice[p][i]) );                      // 16
   plus_equal(M, &aux);
   }


// compute the mean plaquettes (spatial, temporal)
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               GParam const * const param,
               double *plaqs,
               double *plaqt)
   {
   long r;
   double ps, pt;

   ps=0.0;
   pt=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt) reduction(+ : ps)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      int i, j;
      i=0;
      for(j=1; j<STDIM; j++)
         {
         pt+=plaquettep(GC, geo, param, r, i, j);
         }
     
      for(i=1; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            ps+=plaquettep(GC, geo, param, r, i, j);
            }
         }
      }

   if(STDIM>2)
     {
     ps*=param->d_inv_vol;
     ps/=((double) (STDIM-1)*(STDIM-2)/2);
     }
   else
     {
     ps=0.0;
     }

   pt*=param->d_inv_vol;
   pt/=((double) STDIM-1);

   *plaqs=ps;
   *plaqt=pt;
   }


// compute the clover discretization of
// sum_{\mu\nu}  Tr(F_{\mu\nu}F_{\mu\nu})/2
void clover_disc_energy(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        double *energy)
  {
  long r;
  double ris;

  ris=0.0;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
  #endif
  for(r=0; r<param->d_volume; r++)
     {
     int i, j;
     GAUGE_GROUP aux1, aux2;

     for(i=0; i<STDIM; i++)
        {
        for(j=i+1; j<STDIM; j++)
           {
           clover(GC, geo, param, r, i, j, &aux1);

           ta(&aux1);
           equal(&aux2, &aux1);
           times_equal(&aux1, &aux2);
           ris+=-NCOLOR*retr(&aux1)/16.0;
           }
        }
     }

  *energy=ris*param->d_inv_vol;
  }


// compute the mean Polyakov loop (the trace of)
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *repoly,
              double *impoly)
   {
   long rsp;
   double rep, imp;

   rep=0.0;
   imp=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int i;
      GAUGE_GROUP matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(i=0; i<param->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

      rep+=retr(&matrix);
      imp+=imtr(&matrix);
      }

   *repoly=rep*param->d_inv_space_vol;
   *impoly=imp*param->d_inv_space_vol;
   }


// compute the mean Polyakov loop in the adjoint representation (the trace of)
void polyakov_adj(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  double *repoly,
                  double *impoly)
   {
   long rsp;
   double rep, imp;
   double complex tr;

   rep=0.0;
   imp=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int i;
      GAUGE_GROUP matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(i=0; i<param->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }
      tr=NCOLOR*retr(&matrix)+NCOLOR*imtr(&matrix)*I;

      #if NCOLOR==1
        rep+=0.0;
      #else
        rep+=(cabs(tr)*cabs(tr)-1)/(NCOLOR*NCOLOR-1);
      #endif

      imp+=0.0;
      }

   *repoly=rep*param->d_inv_space_vol;
   *impoly=imp*param->d_inv_space_vol;
   }



// compute the mean Polyakov loop and its powers (trace of) in the presence of trace deformation
void polyakov_with_tracedef(Gauge_Conf const * const GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            double *repoly,
                            double *impoly)
   {
   long rsp;
   double **rep, **imp;
   int j, err;
   long i;

   for(j=0;j<(int)floor(NCOLOR/2);j++)
      {
      repoly[j]=0.0;
      impoly[j]=0.0;
      }

   err=posix_memalign((void**)&rep, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double*));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   err=posix_memalign((void**)&imp, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double*));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   for(i=0; i<param->d_space_vol; i++)
      {
      err=posix_memalign((void**)&(rep[i]), (size_t)DOUBLE_ALIGN, (size_t) (int)floor(NCOLOR/2) * sizeof(double));
      if(err!=0)
        {
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      err=posix_memalign((void**)&(imp[i]), (size_t)DOUBLE_ALIGN, (size_t) (int)floor(NCOLOR/2) * sizeof(double));
      if(err!=0)
        {
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   for(i=0; i<param->d_space_vol; i++)
      {
      for(j=0; j<(int)floor(NCOLOR/2); j++)
         {
         rep[i][j] = 0.0;
         imp[i][j] = 0.0;
         }
      }

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int k;
      GAUGE_GROUP matrix, matrix2;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(k=0; k<param->d_size[0]; k++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

       rep[rsp][0] = retr(&matrix);
       imp[rsp][0] = imtr(&matrix);

       equal(&matrix2, &matrix);

      for(k=1; k<(int)floor(NCOLOR/2.0); k++)
         {
         times_equal(&matrix2, &matrix);
         rep[rsp][k] = retr(&matrix2);
         imp[rsp][k] = imtr(&matrix2);
         }
      }

    for(j=0; j<(int)floor(NCOLOR/2); j++)
       {
       for(i=0; i<param->d_space_vol; i++)
          {
          repoly[j] += rep[i][j];
          impoly[j] += imp[i][j];
          }
       }

   for(j=0; j<(int)floor(NCOLOR/2.0); j++)
      {
      repoly[j] *= param->d_inv_space_vol;
      impoly[j] *= param->d_inv_space_vol;
      }

   for(i=0; i<param->d_space_vol; i++)
      {
      free(rep[i]);
      free(imp[i]);
      }
   free(rep);
   free(imp);
   }


// compute the local topological charge at point r
// see readme for more details
double loc_topcharge(Gauge_Conf const * const GC,
                     Geometry const * const geo,
                     GParam const * const param,
                     long r)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     (void) GC;
     (void) geo;
     (void) param;
     (void) r;
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double ris;

   #if (STDIM==4 && NCOLOR>1)
     GAUGE_GROUP aux1, aux2, aux3;
     double real1, real2, loc_charge;
     const double chnorm=1.0/(128.0*PI*PI);
     int i, dir[4][3], sign;

     dir[0][0] = 0;
     dir[0][1] = 0;
     dir[0][2] = 0;

     dir[1][0] = 1;
     dir[1][1] = 2;
     dir[1][2] = 3;

     dir[2][0] = 2;
     dir[2][1] = 1;
     dir[2][2] = 1;

     dir[3][0] = 3;
     dir[3][1] = 3;
     dir[3][2] = 2;

     sign=-1;
     loc_charge=0.0;

     for(i=0; i<3; i++)
        {
        clover(GC, geo, param, r, dir[0][i], dir[1][i], &aux1);
        clover(GC, geo, param, r, dir[2][i], dir[3][i], &aux2);

        times_dag2(&aux3, &aux2, &aux1); // aux3=aux2*(aux1^{dag})
        real1=retr(&aux3)*NCOLOR;

        times(&aux3, &aux2, &aux1); // aux3=aux2*aux1
        real2=retr(&aux3)*NCOLOR;

        loc_charge+=((double) sign*(real1-real2));
        sign=-sign;
        }
     ris=(loc_charge*chnorm);
   #endif

   #if (STDIM==2 && NCOLOR==1)
     GAUGE_GROUP u1matrix;
     double angle;

     plaquettep_matrix(GC, geo, param, r, 0, 1, &u1matrix);
     angle=atan2(cimag(u1matrix.comp), creal(u1matrix.comp))/PI2;

     ris=angle;
   #endif

   return ris;
   }


// compute the topological charge
// see readme for more details
double topcharge(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double ris;
   long r;

   ris=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      ris+=loc_topcharge(GC, geo, param, r);
      }

   return ris;
   }

// chi^\prime = (1/8) int d^4x |x|^2 <q(x)q(0)> = < (1/8) int d^4x |x|^2 q(x) q(0) > = < G2 >
// This function computes the quantity (q(0)/8) sum_{x} d(x,0)^2 q(x) = a^2 G2, whose mean over the ensamble is a^2 chi^\prime
// d(x,y) = lattice distance between sites x and y keeping periodic boundary conditions into account (i.e., the shortest distance between x and y)
double topo_chi_prime(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	double ris=0.0, factor=0.125; // factor = 1/(2D) = 1/8
	long r;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:ris)
	#endif
	for(r=0; r<(param->d_volume); r++)
	{
		double d2 = square_distance(r, 0, param); // d(r,0)^2
		ris += d2 * loc_topcharge(GC, geo, param, r);
	}
	ris *=  loc_topcharge(GC, geo, param, 0) * factor; // ris *= q(0) / 8
	
	return ris;
}

void topcharge_timeslices(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param, double *ris, int ncool, FILE *topchar_tcorr_filep)
{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	long r;
	int N_t = param->d_size[0];

	for (int i=0; i<N_t; i++) ris[i]=0.0;
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:ris[:N_t])
	#endif
	for(r=0; r<(param->d_volume); r++)
	{
		int t = geo->d_timeslice[r];
		ris[t] += loc_topcharge(GC, geo, param, r);
	}

	fprintf(topchar_tcorr_filep, "%ld %d ", GC->update_index, ncool);
	for (int i=0; i<param->d_size[0]; i++) fprintf(topchar_tcorr_filep, " %.12g", ris[i]);
	fprintf(topchar_tcorr_filep, "\n");
}

void topcharge_timeslices_cooling(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param, FILE *topchar_tcorr_filep)
{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	
	double *sum_q_timeslices;

	int err=posix_memalign((void**) &(sum_q_timeslices), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0]*sizeof(double));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating the aux vector for topcharge tcorr meas! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}	

	if(param->d_coolsteps>0)  // if using cooling
	{  
		Gauge_Conf helperconf;
		int iter;

		// measure no cooling
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep); 
		// conf that will be cooled
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param); // helperconf is a copy of the configuration
		// measure with cooling
		for(iter=0; iter<(param->d_coolrepeat); iter++)
		{
			cooling(&helperconf, geo, param, param->d_coolsteps);
			topcharge_timeslices(&helperconf, geo, param, sum_q_timeslices, (iter+1)*param->d_coolsteps, topchar_tcorr_filep);
		}
		free_gauge_conf(&helperconf, param);
		fflush(topchar_tcorr_filep);
	}
	else // no cooling
	{
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep);
		fflush(topchar_tcorr_filep);
	}
	free(sum_q_timeslices);
}

void topcharge_timeslices_cooling_beta_pt(Gauge_Conf const * const GC_vec,
                                          Geometry const * const geo,
                                          GParam const * const param_vec, FILE *topchar_tcorr_filep)
{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
for (int rep_index = 0; rep_index < param_vec->d_N_replica_pt; rep_index++) {
   Gauge_Conf const * GC = GC_vec + rep_index;
   GParam const * param = param_vec + rep_index;
	
	double *sum_q_timeslices;

	int err=posix_memalign((void**) &(sum_q_timeslices), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0]*sizeof(double));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating the aux vector for topcharge tcorr meas! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}	

	if(param->d_coolsteps>0)  // if using cooling
	{  
		Gauge_Conf helperconf;
		int iter;

		// measure no cooling
      fprintf(topchar_tcorr_filep, "    beta = %.12g\n", param->d_beta);
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep); 
		// conf that will be cooled
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param); // helperconf is a copy of the configuration
		// measure with cooling
		for(iter=0; iter<(param->d_coolrepeat); iter++)
		{
			cooling(&helperconf, geo, param, param->d_coolsteps);
			topcharge_timeslices(&helperconf, geo, param, sum_q_timeslices, (iter+1)*param->d_coolsteps, topchar_tcorr_filep);
		}
		free_gauge_conf(&helperconf, param);
		fflush(topchar_tcorr_filep);
	}
	else // no cooling
	{
		topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tcorr_filep);
		fflush(topchar_tcorr_filep);
	}
	free(sum_q_timeslices);
}
}

// compute topological observables (Q, chi_prime) after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topo_obs_cooling(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
											 double *chi_prime,
                       double *meanplaq)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   if(param->d_coolsteps>0)  // if using cooling
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
     // helperconf is a copy of the configuration
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        cooling(&helperconf, geo, param, param->d_coolsteps);

        ris=topcharge(&helperconf, geo, param);
        charge[iter]=ris;

				ris=topo_chi_prime(&helperconf, geo, param);
				chi_prime[iter]=ris;

        plaquette(&helperconf, geo, param, &plaqs, &plaqt);
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }

     free_gauge_conf(&helperconf, param); 
     }
   else   // no cooling
     {
     double ris, ris2, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, geo, param);
		 ris2=topo_chi_prime(GC, geo, param);
     plaquette(GC, geo, param, &plaqs, &plaqt);
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        charge[iter]=ris;
				chi_prime[iter]=ris2;
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }
     } 
   }

/*---------------------------------------------*/
// OBSERVABLE NEEDED JUST TO CHECK HOW COOLING DESTROYS TOPOLOGICAL CORRELATIONS
void check_correlation_decay_cooling(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param, double *ratio)
{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	if(param->d_coolsteps>0)  // if using cooling
	{  
		Gauge_Conf helperconf; 
		double Q, satd;
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		for(int iter=0; iter<(param->d_coolrepeat); iter++)
		{
			cooling(&helperconf, geo, param, param->d_coolsteps);
			Q = fabs(topcharge(&helperconf, geo, param));
			satd = sum_abs_topcharge_dens(&helperconf, geo, param);
			ratio[iter] = (satd-Q)/satd;
		}
		free_gauge_conf(&helperconf, param); 
	}
}

double sum_abs_topcharge_dens(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param)
{
	double sum=0.0;
	for (long r=0; r<(param->d_volume); r++)
	{
		sum += fabs(loc_topcharge(GC, geo, param, r));
	}
	return sum;
}

/*---------------------------------------------*/

// compute the topological charge after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topcharge_cooling(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
                       double *meanplaq)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   if(param->d_coolsteps>0)  // if using cooling
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
     // helperconf is a copy of the configuration
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        cooling(&helperconf, geo, param, param->d_coolsteps);

        ris=topcharge(&helperconf, geo, param);
        charge[iter]=ris;

        plaquette(&helperconf, geo, param, &plaqs, &plaqt);
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }

     free_gauge_conf(&helperconf, param); 
     }
   else   // no cooling
     {
     double ris, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, geo, param);
     plaquette(GC, geo, param, &plaqs, &plaqt);
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        charge[iter]=ris;
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }
     } 
   }


// compute the correlator of the local topological charge
// after "ncool" cooling steps up to spatial distance "dist"
void loc_topcharge_corr(Gauge_Conf const * const GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    int ncool,
                    int dist,
                    double *ris)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double *topch;
   long r;
   int err, i;

   err=posix_memalign((void**) &(topch), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating memory! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   // compute the local topological charge
   if(ncool>0)
     {
     Gauge_Conf helperconf;

     // helperconf is a copy of GC
     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);

     // cool helperconf
     cooling(&helperconf, geo, param, ncool);

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_volume; r++)
        {
        topch[r]=loc_topcharge(&helperconf, geo, param, r);
        }

     // free helperconf
     free_gauge_conf(&helperconf, param);
     }
   else
     {
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_volume; r++)
        {
        topch[r]=loc_topcharge(GC, geo, param, r);
        }
     }

   // compute correlators
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(i)
   #endif
   for(i=0; i<dist; i++)
      {
      double aux;
      long r1, r2;
      int j, dir;

      ris[i]=0.0;

      for(r1=0; r1<param->d_volume; r1++)
         {
         aux=0.0;

         for(dir=1; dir<STDIM; dir++)
            {
            r2=r1;
            for(j=0; j<i; j++)
               {
               r2=nnp(geo, r2, dir);
               }
            aux+=topch[r2];
            }
         aux/=(double)(STDIM-1);

         ris[i]+=aux*topch[r1];
         }
      ris[i]*=param->d_inv_vol;
      }

   // free memory
   free(topch);
   }

//stupid, NOT PARALLELIZED implementation
void perform_measures_beta_pt_replica(Gauge_Conf *GC_vec,
                                      Geometry const * const geo,
                                      GParam const * const param_vec,
                                      FILE* datafilep, FILE* chiprimefilep)
{
   for (int rep_index = 0; rep_index < param_vec->d_N_replica_pt; rep_index++) {
     GParam const * param = param_vec + rep_index;
     Gauge_Conf* GC = GC_vec + rep_index;
#if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
     int i, err;
     double plaqs, plaqt, polyre, polyim, *charge, *chi_prime, *meanplaq, charge_nocooling, chi_prime_nocooling;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);
		 
     charge_nocooling=topcharge(GC, geo, param);

		 // refresh topological charge of periodic replica (only for multicanonic)
		 GC->stored_topo_charge = charge_nocooling;

     fprintf(datafilep, "%.12g %ld %.12g %.12g %.12g %.12g %.12g ", param->d_beta, GC->update_index, plaqs, plaqt, polyre, polyim, charge_nocooling);
		 if (param->d_chi_prime_meas == 1 ) {
         chi_prime_nocooling=topo_chi_prime(GC, geo, param);
         fprintf(chiprimefilep, "%.12g %ld 0 %.12lg\n",param->d_beta, GC->update_index, chi_prime_nocooling);
       }

     err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

		 if (param->d_chi_prime_meas == 1)
		 {
     	err=posix_memalign((void**)&chi_prime, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     	if(err!=0)
       	{
       	fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       	exit(EXIT_FAILURE);
       	}
		 }

     err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     if (param->d_chi_prime_meas == 1 ) topo_obs_cooling(GC, geo, param, charge, chi_prime, meanplaq);
		 else topcharge_cooling(GC, geo, param, charge, meanplaq);
     for(i=0; i<param->d_coolrepeat; i++)
        {
        fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
				if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld %d %.12lg\n", GC->update_index, (i+1)*param->d_coolsteps, chi_prime[i]);
        }
     fprintf(datafilep, "\n");

     fflush(datafilep);
		 if (param->d_chi_prime_meas == 1 ) fflush(chiprimefilep);

     free(charge);
		 if (param->d_chi_prime_meas == 1 ) free(chi_prime);
     else 
			{
				(void) chiprimefilep;
			}
		 free(meanplaq);

#else

     double plaqs, plaqt, polyre, polyim;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);

     fprintf(datafilep, "%.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim);
     fprintf(datafilep, "\n");
     fflush(datafilep);
#endif
   }
}                                      

void perform_measures_localobs(Gauge_Conf *GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep, FILE *chiprimefilep, FILE *topchar_tcorr_filep)
   {
   #if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
     int i, err;
     double plaqs, plaqt, polyre, polyim, *charge, *chi_prime, *meanplaq, charge_nocooling, chi_prime_nocooling;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);
		 
     charge_nocooling=topcharge(GC, geo, param);

		 // refresh topological charge of periodic replica (only for multicanonic)
		 GC->stored_topo_charge = charge_nocooling;

     fprintf(datafilep, "%ld %.12g %.12g %.12g %.12g %.12g ", GC->update_index, plaqs, plaqt, polyre, polyim, charge_nocooling);
		 if (param->d_chi_prime_meas == 1 ) {
         chi_prime_nocooling=topo_chi_prime(GC, geo, param);
         fprintf(chiprimefilep, "%ld 0 %.12lg\n", GC->update_index, chi_prime_nocooling);
       }

     err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

		 if (param->d_chi_prime_meas == 1)
		 {
     	err=posix_memalign((void**)&chi_prime, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     	if(err!=0)
       	{
       	fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       	exit(EXIT_FAILURE);
       	}
		 }

     err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

		 if (param->d_topcharge_tcorr_meas == 1 ) topcharge_timeslices_cooling(GC, geo, param, topchar_tcorr_filep);
		 else {(void) topchar_tcorr_filep;}
     if (param->d_chi_prime_meas == 1 ) topo_obs_cooling(GC, geo, param, charge, chi_prime, meanplaq);
		 else topcharge_cooling(GC, geo, param, charge, meanplaq);
     for(i=0; i<param->d_coolrepeat; i++)
        {
        fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
				if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld %d %.12lg\n", GC->update_index, (i+1)*param->d_coolsteps, chi_prime[i]);
        }
     fprintf(datafilep, "\n");

     fflush(datafilep);
		 if (param->d_chi_prime_meas == 1 ) fflush(chiprimefilep);

     free(charge);
		 if (param->d_chi_prime_meas == 1 ) free(chi_prime);
     else 
			{
				(void) chiprimefilep;
			}
		 free(meanplaq);

   #else

     double plaqs, plaqt, polyre, polyim;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);

     fprintf(datafilep, "%.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim);
     fprintf(datafilep, "\n");
     fflush(datafilep);
   #endif
   }

// perform local observables in the case of trace deformation, it computes all the order parameters
void perform_measures_localobs_with_tracedef(Gauge_Conf const * const GC,
                                             Geometry const * const geo,
                                             GParam const * const param,
                                             FILE *datafilep)
   {
   #if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )

     int i, err;
     double plaqs, plaqt, polyre[NCOLOR/2+1], polyim[NCOLOR/2+1]; // +1 just to avoid warning if NCOLOR=1
     double *charge, *meanplaq, charge_nocooling;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov_with_tracedef(GC, geo, param, polyre, polyim);
     charge_nocooling=topcharge(GC, geo, param);

     fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);

     for(i=0; i<(int)floor(NCOLOR/2); i++)
        {
        fprintf(datafilep, "%.12g %.12g ", polyre[i], polyim[i]);
        }
     fprintf(datafilep, "%.12g ", charge_nocooling);

     err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     topcharge_cooling(GC, geo, param, charge, meanplaq);
     for(i=0; i<param->d_coolrepeat; i++)
        {
        fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
        }
     fprintf(datafilep, "\n");

     fflush(datafilep);

     free(charge);
     free(meanplaq);

   #else

     double plaqs, plaqt, polyre, polyim;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);

     fprintf(datafilep, "%.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim);
     fprintf(datafilep, "\n");
     fflush(datafilep);

   #endif
   }


// to optimize the number of hits to be used in multilevel
void optimize_multihit_polycorr(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
  {
  const int max_hit=50;
  const int dir=1;

  int i, mh, t_tmp, err;
  long r, r1, r2;
  double complex poly_corr;
  double poly_corr_abs, poly_corr_fluct, diff_sec;
  double complex *poly_array;
  time_t time1, time2;
  GAUGE_GROUP matrix, tmp;

  err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  #ifdef THETA_MODE
   compute_clovers(GC, geo, param, 0);
  #endif

  fprintf(datafilep, "Multihit optimization: \n");
  fprintf(datafilep, "the smaller the value the better the multihit\n");

  for(mh=1; mh<max_hit; mh++)
     {
     time(&time1);

     // polyakov loop computation
     for(r=0; r<param->d_space_vol; r++)
        {
        one(&matrix);
        for(i=0; i<param->d_size[0]; i++)
           {
           multihit(GC,
                    geo,
                    param,
                    sisp_and_t_to_si(geo, r, i),
                    0,
                    mh,
                    &tmp);
           times_equal(&matrix, &tmp);
           }
        poly_array[r]=retr(&matrix)+I*imtr(&matrix);
        }

     // average correlator computation
     poly_corr=0.0+I*0.0;
     poly_corr_abs=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

        poly_corr += poly_array[r]*conj(poly_array[r2]);
        poly_corr_abs += cabs(poly_array[r]*conj(poly_array[r2]));
        }
     poly_corr*=param->d_inv_space_vol;
     poly_corr_abs*=param->d_inv_space_vol;

     // fluctuation of the average correlator computation
     poly_corr_fluct=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1
        poly_corr_fluct+=cabs( poly_array[r]*conj(poly_array[r2]) - poly_corr );
        }
     poly_corr_fluct*=param->d_inv_space_vol;


     time(&time2);
     diff_sec = difftime(time2, time1);

     fprintf(datafilep, "%d  %.12g  %.12g (time:%g)\n", mh, poly_corr_abs*sqrt(mh), poly_corr_fluct*sqrt(mh), diff_sec);

     fflush(datafilep);
     }

  free(poly_array);
  }


// to optimize the multilevel
void optimize_multilevel_polycorr(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep)
   {
   int i, err;
   long r;
   double complex poly_corr;
   double poly_corr_abs, poly_corr_fluct;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel_polycorr(GC,
                       geo,
                       param,
                       param->d_size[0]);
   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   // averages
   poly_corr=0.0+I*0.0;
   poly_corr_abs=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=cabs(poly_array[r]);
      }
   poly_corr*=param->d_inv_space_vol;
   poly_corr_abs*=param->d_inv_space_vol;

   // fluctuations
   poly_corr_fluct=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_corr_fluct += cabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=param->d_inv_space_vol;

   // normalizations
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*= sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*= sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g ", poly_corr_abs);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "\n");
   fflush(datafilep);

   free(poly_array);
   }


// perform the computation of the polyakov loop correlator with the multilevel algorithm
void perform_measures_polycorr(Gauge_Conf *GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep)
   {
   #ifndef OPT_MULTIHIT
   #ifndef OPT_MULTILEVEL
     double ris;
     long r;
     int i;

     multilevel_polycorr(GC,
                geo,
                param,
                param->d_size[0]);

     for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   #endif

   #ifdef OPT_MULTIHIT
     optimize_multihit_polycorr(GC, geo, param, datafilep);
   #endif

   #ifdef OPT_MULTILEVEL
     optimize_multilevel_polycorr(GC, geo, param, datafilep);
   #endif
   }


// to optimize the number of hits to be used in multilevel for the adjoint representation
void optimize_multihit_polycorradj(Gauge_Conf *GC,
                                   Geometry const * const geo,
                                   GParam const * const param,
                                   FILE *datafilep)
  {
  const int max_hit=50;
  const int dir=1;

  int i, mh, t_tmp, err;
  long r, r1, r2;
  double poly_corr, poly_corr_abs, poly_corr_fluct, diff_sec;
  double complex tr;
  double *poly_array;
  time_t time1, time2;
  GAUGE_GROUP matrix, tmp;

  err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  #ifdef THETA_MODE
   compute_clovers(GC, geo, param, 0);
  #endif

  fprintf(datafilep, "Multihit optimization: \n");
  fprintf(datafilep, "the smaller the value the better the multihit\n");

  for(mh=1; mh<max_hit; mh++)
     {
     time(&time1);

     // polyakov loop in the adjoint representation computation
     for(r=0; r<param->d_space_vol; r++)
        {
        one(&matrix);
        for(i=0; i<param->d_size[0]; i++)
           {
           multihit(GC,
                    geo,
                    param,
                    sisp_and_t_to_si(geo, r, i),
                    0,
                    mh,
                    &tmp);
           times_equal(&matrix, &tmp);
           }

        //trace of the matix in the fundamental representation
        tr=NCOLOR*(retr(&matrix)+I*imtr(&matrix));

        //trace of the matrix in adjoint representation
        poly_array[r]=cabs(tr*conj(tr))-1.0;
        #if NCOLOR != 1
          poly_array[r]/=(NCOLOR*NCOLOR-1);
        #endif
        }

     // average correlator computation
     poly_corr=0.0;
     poly_corr_abs=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

        poly_corr+= poly_array[r]*poly_array[r2];
        poly_corr_abs+=fabs(poly_array[r]*poly_array[r2]);
        }
     poly_corr*=param->d_inv_space_vol;
     poly_corr_abs*=param->d_inv_space_vol;

     // fluctuation of the average correlator computation
     poly_corr_fluct=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

        poly_corr_fluct+=fabs(poly_array[r]*poly_array[r2]-poly_corr);
        }
     poly_corr_fluct*=param->d_inv_space_vol;

     time(&time2);
     diff_sec = difftime(time2, time1);

     fprintf(datafilep, "%d  %.12g  %.12g (time:%g)\n", mh, poly_corr_abs*sqrt(mh), poly_corr_fluct*sqrt(mh), diff_sec);

     fflush(datafilep);
     }

  free(poly_array);
  }


// to optimize the multilevel (adjoint representation)
void optimize_multilevel_polycorradj(Gauge_Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   int i, err;
   long r;
   double poly_corr, poly_corr_abs, poly_corr_fluct;
   double *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel_polycorradj(GC,
                          geo,
                          param,
                          param->d_size[0]);

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProdAdj(&(GC->ml_polycorradj[0][0][r]), &(GC->ml_polycorradj[0][i][r]) );
         }
      }

   // averages
   poly_corr=0.0;
   poly_corr_abs=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProdAdj(&(GC->ml_polycorradj[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=fabs(poly_array[r]);
      }
   poly_corr*=param->d_inv_space_vol;
   poly_corr_abs*=param->d_inv_space_vol;

   // fluctuations
   poly_corr_fluct=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_corr_fluct+=fabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=param->d_inv_space_vol;

   // normalizations
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*=sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*=sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g %.12g ", poly_corr_abs, poly_corr_fluct);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "\n");
   fflush(datafilep);

   free(poly_array);
   }


// perform the computation of the polyakov loop correlator in the adjoint representation with the multilevel algorithm
void perform_measures_polycorradj(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep)
   {
   #ifndef OPT_MULTIHIT
   #ifndef OPT_MULTILEVEL
     double ris;
     long r;
     int i;

     multilevel_polycorradj(GC,
                            geo,
                            param,
                            param->d_size[0]);

     for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           times_equal_TensProdAdj(&(GC->ml_polycorradj[0][0][r]), &(GC->ml_polycorradj[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProdAdj(&(GC->ml_polycorradj[0][0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   #endif

   #ifdef OPT_MULTIHIT
     optimize_multihit_polycorradj(GC, geo, param, datafilep);
   #endif

   #ifdef OPT_MULTILEVEL
     optimize_multilevel_polycorradj(GC, geo, param, datafilep);
   #endif
   }


// to optimize the multilevel
void optimize_multilevel_polycorr_long(Gauge_Conf *GC,
                                       GParam const * const param,
                                       FILE *datafilep)
   {
   int i, err;
   long r;
   double poly_corr_abs, poly_corr_fluct;
   double complex poly_corr;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   // average
   poly_corr=0.0+I*0.0;
   poly_corr_abs=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=cabs(poly_array[r]);
      }
   poly_corr*=param->d_inv_space_vol;
   poly_corr_abs*=param->d_inv_space_vol;

   // fluctuation
   poly_corr_fluct=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_corr_fluct+=cabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=param->d_inv_space_vol;

   // normalization
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*=sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*=sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_ml_level0_repeat);
   poly_corr_fluct*=sqrt(param->d_ml_level0_repeat);

   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g %.12g ", poly_corr_abs, poly_corr_fluct);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "(%d) ", param->d_ml_level0_repeat);
   fprintf(datafilep, "\n");

   fflush(datafilep);

   free(poly_array);
   }


// print the value of the polyakov loop correlator that has been computed by multilevel
void perform_measures_polycorr_long(Gauge_Conf *GC,
                                    GParam const * const param,
                                    FILE *datafilep)
   {
   #ifdef OPT_MULTILEVEL
      optimize_multilevel_polycorr_long(GC, param, datafilep);
   #else
     double ris;
     long r;
     int i;

     for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   }


// perform the computation of the string width with the
// disconnected correlator using the multilevel algorithm
void perform_measures_tube_disc(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   multilevel_tube_disc(GC,
                        geo,
                        param,
                        param->d_size[0]);

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }


// perform the computation of the string width with the
// connected correlator using the multilevel algorithm
void perform_measures_tube_conn(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   multilevel_tube_conn(GC,
                        geo,
                        param,
                        param->d_size[0]);

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }


// print the value of the the string width with the
// connected correlator that has been computed by multilevel
void perform_measures_tube_conn_long(Gauge_Conf *GC,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }

#endif

void measure_poly_profile(Gauge_Conf *GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          FILE* datafilep)
{
   double* poly_re, * poly_im;
   int x_size = param->d_size[1];
   int err = posix_memalign((void**) &poly_re, DOUBLE_ALIGN, (size_t) x_size * sizeof(double));
   if (err != 0) {
      fprintf(stderr, "Unable to allocate memory");
      exit(EXIT_FAILURE);
   }
   err = posix_memalign((void**) &poly_im, DOUBLE_ALIGN, (size_t) x_size * sizeof(double));
   if (err != 0) {
      fprintf(stderr, "Unable to allocate memory");
      exit(EXIT_FAILURE);
   }

   int x;
   for (x = 0; x < x_size; x++) {
      poly_re[x] = 0;
      poly_im[x] = 0;
   }

   long rsp;
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+:poly_re[:x_size]) reduction(+:poly_im[:x_size])
   #endif
   for (rsp = 0; rsp < param->d_space_vol; rsp++) {
      long r;
      GAUGE_GROUP matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(int i=0; i<param->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

      int cart[STDIM];
      si_to_cart(cart, r, param);
      poly_re[cart[1]] += retr(&matrix);
      poly_im[cart[1]] += imtr(&matrix);
   }

   fprintf(datafilep, "%ld ", GC->update_index);
   for(x = 0; x < param->d_size[1]; x++) {
      fprintf(datafilep, "%.15f %.15f ", poly_re[x], poly_im[x]);
   }
   fprintf(datafilep, "\n");
   fflush(datafilep);

   free(poly_im); free(poly_re);
}                          

void perform_measure_spectrum(Gauge_Conf* GC,
                              Geometry const * const geo,
                              GParam const * const param,
                              FILE* poly_profile_filep,
                              FILE* plaq_profile_filep)
{
   double* poly_re, * poly_im;
   int t_size = param->d_size[0];
   int err = posix_memalign((void**) &poly_re, DOUBLE_ALIGN, (size_t) t_size * sizeof(double));
   if (err != 0) {
      fprintf(stderr, "Unable to allocate memory (%s, %d)", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   err = posix_memalign((void**) &poly_im, DOUBLE_ALIGN, (size_t) t_size * sizeof(double));
   if (err != 0) {
      fprintf(stderr, "Unable to allocate memory (%s, %d)", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   double* plaq;
   err = posix_memalign((void**) &plaq, DOUBLE_ALIGN, (size_t) t_size * sizeof(double));
   if (err != 0) {
      fprintf(stderr, "Unable to allocate memory (%s, %d)", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

   // number of blocking steps, must be >1
   #define NUMBLOCK 2
   // smearing parameter
   const double alpha=0.5;
   // smearing steps
   const int smearing_steps = 3;

   Gauge_Conf smeared_GC;
   init_gauge_conf_from_gauge_conf(&smeared_GC, GC, param);
   spatial_smearing(&smeared_GC, geo, param, alpha, smearing_steps);

   poly_profile_mean_winding(&smeared_GC, geo, param, poly_re, poly_im);
   plaq_profile(&smeared_GC, geo, param, plaq);

   print_poly_profile(poly_re, poly_im, t_size, GC->update_index, smearing_steps, 0, poly_profile_filep);
   print_plaq_profile(plaq, t_size, GC->update_index, smearing_steps, 0, plaq_profile_filep);

   GParam blocked_param[NUMBLOCK];
   Geometry blocked_geo[NUMBLOCK];
   Gauge_Conf blocked_GC[NUMBLOCK];

   init_spatial_blocked_conf(blocked_GC, blocked_geo, blocked_param,
                             GC, geo, param, alpha);

   poly_profile_mean_winding(blocked_GC, blocked_geo, blocked_param, poly_re, poly_im);
   plaq_profile(blocked_GC, blocked_geo, blocked_param, plaq);

   print_poly_profile(poly_re, poly_im, t_size, GC->update_index, smearing_steps, 1, poly_profile_filep);
   print_plaq_profile(plaq, t_size, GC->update_index, smearing_steps, 1, plaq_profile_filep);

   for (int i = 1; i < NUMBLOCK; i++) {
      init_spatial_blocked_conf(blocked_GC + i, blocked_geo + i, blocked_param + i, 
                                blocked_GC + i-1, blocked_geo + i-1, blocked_param + i-1,
                                alpha);

      poly_profile_mean_winding(blocked_GC + i, blocked_geo + i, blocked_param + i, poly_re, poly_im);
      plaq_profile(blocked_GC + i, blocked_geo + i, blocked_param + i, plaq);

      print_poly_profile(poly_re, poly_im, t_size, GC->update_index, smearing_steps, i+1, poly_profile_filep);
      print_plaq_profile(plaq, t_size, GC->update_index, smearing_steps, i+1, plaq_profile_filep);
   }
   
   free_gauge_conf(&smeared_GC, param);
   
   for (int i = 0; i < NUMBLOCK; i++) {
      free_gauge_conf(blocked_GC + i, blocked_param + i);
      free_geometry(blocked_geo + i, blocked_param + i);
   }
}

void print_poly_profile(double* poly_re, double* poly_im, int t_size,
                        long update, int smearing, int blocking, FILE* datafilep) {
   fprintf(datafilep, "%ld %d %d ", update, smearing, blocking);
   for (int t = 0; t < t_size; t++)
      fprintf(datafilep, "%.15f %.15f ", poly_re[t], poly_im[t]);
   fprintf(datafilep, "\n");
   fflush(datafilep);
}

void print_plaq_profile(double* plaq, int t_size,
                        long update, int smearing, int blocking, FILE* datafilep) {
   fprintf(datafilep, "%ld %d %d ", update, smearing, blocking);
   for (int t = 0; t < t_size; t++)
      fprintf(datafilep, "%.15f ", plaq[t]);
   fprintf(datafilep, "\n");
   fflush(datafilep);
}

void poly_profile_mean_winding(Gauge_Conf* GC,
                               Geometry const * const geo, 
                               GParam const * const param,
                               double* result_re,
                               double* result_im)
{
   int t_size = param->d_size[0];
   for(int t = 0; t < t_size; t++) {
      result_re[t] = 0;
      result_im[t] = 0;
   }

   double* poly_re, * poly_im;

   int err = posix_memalign(&poly_re, (size_t) DOUBLE_ALIGN, (size_t) t_size * sizeof(double));
   if (err != 0) {
      fprintf(stderr, "problems allocating memoery (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   err = posix_memalign(&poly_im, (size_t) DOUBLE_ALIGN, (size_t) t_size * sizeof(double));
   if (err != 0) {
      fprintf(stderr, "problems allocating memoery (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

   long r_hs;
   for (int dir = 1; dir < STDIM; dir++) { //looping over winding dimension
      long hypersurf = t_size; // computing the hypersurface volume 
      for (int j_dir = 1; j_dir < STDIM; j_dir++) {
         if (j_dir != dir) hypersurf *= param->d_size[j_dir];
      }

      for (int t = 0; t < t_size; t++) {
         poly_re[t] = 0.;
         poly_im[t] = 0.;
      }
      
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r_hs) reduction(+:poly_re[:t_size]) reduction(+:poly_im[:t_size])
      #endif
      for (r_hs = 0; r_hs < hypersurf; r_hs++) {
         int cart[STDIM];
         long aux_hs = r_hs;
         for (int j_dir = 0; j_dir < STDIM; j_dir++) { //compute cartesian
            if (j_dir == dir) cart[j_dir] = 0;
            else {
               cart[j_dir] = (int)(aux_hs % param->d_size[j_dir]);
               aux_hs /= param->d_size[j_dir];
            }
         }
         long r;
         GAUGE_GROUP matrix;
         r = cart_to_si(cart, param);

         one(&matrix);
         for (int i = 0; i < param->d_size[dir]; i++) { //winding
            times_equal(&matrix, &(GC->lattice[r][dir]));
            r=nnp(geo, r, dir);
         }

         poly_re[cart[0]] += retr(&matrix);
         poly_im[cart[0]] += imtr(&matrix);
      }
      for (int t = 0; t < t_size; t++) { //Different ways to do this. Ask Bonanno which is better
         result_re[t] += poly_re[t] / hypersurf;
         result_im[t] += poly_im[t] / hypersurf;
      }
   }
   for (int t = 0; t < t_size; t++) { //normalize O(t) = sum_i O_i(t) / (dim - 1)  
      result_re[t] /= STDIM - 1;
      result_im[t] /= STDIM - 1;
   }
}

void plaq_profile(Gauge_Conf* GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  double* result)
{
   int t_size = param->d_size[0];
   for (int t = 0; t < t_size; t++)
      result[t] = 0.;
   
   long r;
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:plaq[:t_size])
   #endif
   for (r = 0; r < param->d_volume; r++) {
      long rsp;
      int t;
      si_to_sisp_and_t(&rsp, &t, geo, r);
      for (int i = 0; i < STDIM; i++) {
         for (int j = i+1; j < STDIM; j++) {
            result[t] += plaquettep(GC, geo, param, r, i, j);
         }
      }
   }

   for (int t = 0; t < t_size; t++)
      result[t] /= 0.5 * (double)(param->d_space_vol * (STDIM - 1) * (STDIM - 2));
}                  

// perform smearing on spatial links
void spatial_smearing(Gauge_Conf* GC,
                      Geometry const * const geo,
                      GParam const * const param,
                      double alpha,
                      int smearing_steps)
  {
  int i, step;
  long r;
  Gauge_Conf staple_GC;
  
  init_gauge_conf_from_gauge_conf(&staple_GC, GC, param);

  for(step=0; step<smearing_steps; step++)
     {
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r = 0; r < param->d_volume; r++)
        {
        GAUGE_GROUP matrix;
        for(i = 1; i < STDIM; i++)
           {
           calcstaples_wilson_no_time(GC, geo, r, i, &matrix);
           equal(&(staple_GC.lattice[r][i]), &matrix);
           }
        }
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r = 0; r < param->d_volume; r++)
        {
        for(i = 1; i < STDIM; i++)
           {
           times_equal_real(&(staple_GC.lattice[r][i]), alpha);
           plus_equal_dag(&GC->lattice[r][i], &(staple_GC.lattice[r][i]));
           unitarize(&GC->lattice[r][i]); 
           }
        }
     }

  free_gauge_conf(&staple_GC, param);
  }

// compute blocked link for a given site
void spatialblocking_singlelink(Gauge_Conf const * const GC,
                                Geometry const * const geo,
                                long r,
                                int i,
                                double blockcoeff,
                                GAUGE_GROUP *M)
  {
  int j, l;
  long k, k1;

  GAUGE_GROUP link1, link2, link3, link4, staptmp, stap;

  #if DEBUG
  if(r >= geo->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i == 0)
    {
    fprintf(stderr, "time direction selected: i=%d (%s, %d)\n", i, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  zero(&stap); // M=0

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

     if(j!=0)
       {
//
//       i ^
//         |    (1)
//         +----->-----+
//         |           |
//         |           V (2)
//         |           |
//         |           |
//       k +-----------+
//         |           |
//         |           |
//         |           V (3)
//         |           |
//         +-----<-----+-->   j
//       r     (4)
//

       k=nnp(geo, r, i);
       equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
       equal(&link2, &(GC->lattice[nnp(geo, k, j)][i]));  // link2 = (2)
       equal(&link3, &(GC->lattice[nnp(geo, r, j)][i]));  // link3 = (3)
       equal(&link4, &(GC->lattice[r][j]));               // link3 = (4)

       times_dag2(&staptmp, &link1, &link2);   // staptmp=link1*link2^{dag}
       times_equal_dag(&staptmp, &link3);      // staptmp*=link3^{dag}
       times_equal_dag(&staptmp, &link4);      // staptmp*=link4^{dag}

       plus_equal(&stap, &staptmp);

//
//       i ^
//         |   (1)
//         +----<------+
//         |           |
//     (2) V           |
//         |           |
//      k1 +-----------+
//         |           |
//     (3) V           |
//         |           |
//         +------>----+--->j
//        k     (4)    r
//

       k=nnm(geo, r, j);
       k1=nnp(geo, k, i);

       equal(&link1, &(GC->lattice[nnp(geo, k1, i)][j]));  // link1 = (1)
       equal(&link2, &(GC->lattice[k1][i]));               // link2 = (2)
       equal(&link3, &(GC->lattice[k][i]));                // link3 = (3)
       equal(&link4, &(GC->lattice[k][j]));                // link4 = (4)

       times_dag12(&staptmp, &link1, &link2); // stap=link1^{dag}*link2^{dag}
       times_equal_dag(&staptmp, &link3);     // stap*=link3^{dag}
       times_equal(&staptmp, &link4);         // stap*=link4

       plus_equal(&stap, &staptmp);
       }
     }

   equal(M, &(GC->lattice[r][i]));
   times_equal(M, &(GC->lattice[nnp(geo, r, i)][i]));

   times_equal_real(&stap, blockcoeff);
   plus_equal_dag(M, &stap);
   unitarize(M);
   }

// compute staples with no time direction
void calcstaples_wilson_no_time(Gauge_Conf const * const GC,
                                Geometry const * const geo,
                                long r,
                                int i,
                                GAUGE_GROUP *M)
  {
  int j, l;
  long k;
  GAUGE_GROUP link1, link2, link3, link12, stap;

  #if DEBUG
  if(r >= geo->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, geo->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i == 0)
    {
    fprintf(stderr, "time direction selected: i=%d (%s, %d)\n", i, __FILE__, __LINE__); 
    exit(EXIT_FAILURE);
    }
  #endif

  zero(M); // M=0

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

     if(j!=0)
       { 

//
//       i ^
//         |   (1)
//         +----->-----+
//         |           |
//                     |
//         |           V (2)
//                     |
//         |           |
//         +-----<-----+-->   j
//       r     (3)
//

       equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
       equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
       equal(&link3, &(GC->lattice[r][j]));               // link3 = (3)

       times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
       times_dag2(&stap, &link12, &link3);   // stap=link12*stap^{dag}

       plus_equal(M, &stap);

//
//       i ^
//         |   (1)
//         |----<------+
//         |           |
//         |
//     (2) V           |
//         |
//         |           |
//         +------>----+--->j
//        k     (3)    r
//

       k=nnm(geo, r, j);

       equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
       equal(&link2, &(GC->lattice[k][i]));               // link2 = (2)
       equal(&link3, &(GC->lattice[k][j]));               // link3 = (3)

       times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
       times(&stap, &link12, &link3);        // stap=link12*link3

       plus_equal(M, &stap);
       }
     }
   }


void init_spatial_blocked_conf(Gauge_Conf* blocked_GC,
                               Geometry* blocked_geo,
                               GParam* blocked_param,
                               Gauge_Conf const * const source_GC,
                               Geometry const * const sourge_geo,
                               GParam const * const source_param,
                               double alpha)
{
   *blocked_param = *source_param;
   for (int dir = 0; dir < STDIM; dir++) {
      if (source_param->d_size[dir] % 2 != 0){
         fprintf(stderr, "Problem with spatial size not even: %d ! (%s, %d)\n", source_param->d_size[dir], __FILE__, __LINE__);
         exit(EXIT_FAILURE);
      }
      blocked_param->d_size[dir] = source_param->d_size[dir] / 2;
   }
   init_derived_constants(blocked_param);

   init_geometry(blocked_geo, blocked_param);
   
   // allocate the lattice
   int err=posix_memalign((void**)&(blocked_GC->lattice), (size_t) DOUBLE_ALIGN, (size_t) blocked_param->d_volume * sizeof(GAUGE_GROUP *));
   if(err!=0)
      {
      fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   for(long r = 0; r < (blocked_param->d_volume); r++)
      {
      err=posix_memalign((void**)&(blocked_GC->lattice[r]), (size_t) DOUBLE_ALIGN, (size_t) STDIM * sizeof(GAUGE_GROUP));
      if(err!=0)
         {
         fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
      }
   
   long rb; // blocked conf position
   /* Commented waiting for Bonanno approoval
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rb)
   #endif */
   for (rb = 0; rb < blocked_param->d_volume; rb++) {
      GAUGE_GROUP matrix;
      for (int dir = 0; dir < STDIM; dir++) {
         int cart[STDIM]; // converting to source conf position
         si_to_cart(cart, rb, blocked_param);
         for (int i = 1; i < STDIM; i++) cart[i] *= 2;
         long r = cart_to_si(cart, source_param); //source conf position
         spatialblocking_singlelink(source_GC, sourge_geo, r, dir, alpha, &matrix);
         equal(&(blocked_GC->lattice[rb][dir]), &matrix);
      }
   }

   blocked_GC->update_index=source_GC->update_index;
}                               
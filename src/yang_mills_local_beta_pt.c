#ifndef YM_LOCAL_BETA_PT_C
#define YM_LOCAL_BETA_PT_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
  #include<omp.h>
#endif

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"

void real_main(char* in_file) {
    Gauge_Conf *GC;
    Geometry geo;
    GParam param;
    Acc_Utils acc_counters;

    char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
    int count;
    FILE *datafilep, *chiprimefilep, *swaptrackfilep, *topchar_tcorr_filep;
    time_t time1, time2;

    // to disable nested parallelism
    #ifdef OPENMP_MODE
      // omp_set_nested(0); // deprecated
			omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
    #endif

    // read input file
    readinput(in_file, &param);

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
	init_data_file(&datafilep, &chiprimefilep, &topchar_tcorr_filep, &param);

	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);

    //init gauge confs replica
    init_gauge_conf_beta_pt(&GC, &param);

    // init acceptances array
	init_swap_acc_arrays(&acc_counters, &param);

	if (param.d_beta_pt == NULL) {
		fprintf(stderr, "Input file didn't contain beta pt info (%s, %d)\n", __FILE__, __LINE__);
		print_template_input();
		exit(EXIT_FAILURE);
	}

    // MONTECARLO BEGIN
    time(&time1);

    for(count = 0; count < param.d_sample; count++) {
        update_beta_pt_replica(GC, &geo, &param);
		if(GC[0].update_index % param.d_measevery == 0 && GC[0].update_index >= param.d_thermal)
		{
			perform_measures_beta_pt_replica(GC, &geo, &param, datafilep, chiprimefilep, topchar_tcorr_filep);
		}

		if (count % param.d_beta_pt_swap_every == 0) {
			beta_pt_swap(GC, &geo, &param, &acc_counters);
			print_conf_labels(swaptrackfilep, GC, &param);
		}
    }

    time(&time2);
    //MONTECARLO END

    // close data file
    fclose(datafilep);
		if (param.d_chi_prime_meas==1) fclose(chiprimefilep);
		if (param.d_topcharge_tcorr_meas==1) fclose(topchar_tcorr_filep);

    // close swap tracking file
	if (param.d_N_replica_pt > 1) fclose(swaptrackfilep);

    // save configurations
    if(param.d_saveconf_back_every!=0)
      {
      write_replica_on_file_beta_pt(GC, &param);
      }

    // print simulation details
    print_parameters_local_pt(&param, time1, time2);
		
	// print acceptances of parallel tempering
	print_acceptances(&acc_counters, &param);

    // free gauge configurations
    free_replica_beta_pt(GC, &param);

    // free geometry
    free_geometry(&geo, &param);

    // free acceptances array
	end_swap_acc_arrays(&acc_counters, &param);
}

void print_template_input(void) {
    //TO DO
}

int main (int argc, char **argv)
{
	char in_file[STD_STRING_LENGTH];
	if(argc != 2)
	{
		printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
		printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
		printf("Usage: %s input_file\n\n", argv[0]);

		printf("Compilation details:\n");
		printf("\tN_c (number of colors): %d\n", NCOLOR);
		printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
		printf("\tNum_levels (number of levels): %d\n", NLEVELS);
		printf("\n");
		printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
		printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

		#ifdef DEBUG
			printf("\n\tDEBUG mode\n");
		#endif

		#ifdef OPENMP_MODE
			printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
		#endif

		#ifdef THETA_MODE
			printf("\n\tusing imaginary theta\n");
		#endif

		printf("\n");

		#ifdef __INTEL_COMPILER
			printf("\tcompiled with icc\n");
		#elif defined(__clang__)
			printf("\tcompiled with clang\n");
		#elif defined( __GNUC__ )
			printf("\tcompiled with gcc version: %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
		#endif

		print_template_input();

		return EXIT_SUCCESS;
	}
	else
	{
		if(strlen(argv[1]) >= STD_STRING_LENGTH)
		{
			fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
		}
		else
		{
			strcpy(in_file, argv[1]);
		}
	}
	real_main(in_file);
	return EXIT_SUCCESS;
}

#endif
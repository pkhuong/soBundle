#define RUNME /*
gcc -std=gnu99 -O2 -ggdb $CFLAGS $0 -W -Wall -L. -I. -lbundle -o `basename $0 .c`
exit $?

Instances:
MO1: Small SPLP (100x100), upper bound ~1305
MS1: Medium SPLP (1000x1000), upper bound ~5285
MT1: Large SPLP (2000x2000), upper bound ~10070
*/
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "bundle_solver.h"

static const double eps = 5e-7;

#define assert_(EXPR)                                           \
        ({                                                      \
                __typeof__(EXPR) __ret = (EXPR);                \
                if (!__ret) {                                   \
                        perror("Assertion failure: " #EXPR);    \
                        assert(0);                              \
                }                                               \
                __ret;                                          \
        })

typedef struct {
        double opening_cost;
        double arc_cost[];
} facility_t;

typedef struct {
        size_t nfacilities;
        size_t ncities;
        facility_t ** facilities;
} instance_t;

void init_problem (instance_t * instance, size_t nfacilities, size_t ncities)
{
        instance->nfacilities = nfacilities;
        instance->ncities     = ncities;
        facility_t ** facilities
                = assert_(malloc(sizeof(facility_t*)*nfacilities));

        size_t facility_size = 
                sizeof(facility_t) + ncities * sizeof(double);

        for (size_t i = 0; i < nfacilities; i++)
                bzero(facilities[i] = assert_(malloc(facility_size)),
                      facility_size);

        instance->facilities = facilities;
}

instance_t read_problem (FILE * in)
{
        size_t nfacilities, ncities;
        assert_(2 == fscanf(in, "%lu %lu", &nfacilities, &ncities));
        instance_t instance;
        init_problem(&instance, nfacilities, ncities);

        for (size_t i = 0; i < nfacilities; i++) {
                double cap, cost;
                assert_(2 == fscanf(in, "%lf %lf", &cap, &cost));
                instance.facilities[i]->opening_cost = cost;
        }

        for (size_t i = 0; i < ncities; i++) {
                double demand;
                assert_(1 == fscanf(in, "%lf", &demand));
                for (size_t j = 0; j < nfacilities; j++) {
                        double connection_cost;
                        assert_(1 == fscanf(in, "%lf", &connection_cost));
                        instance.facilities[j]->arc_cost[i]
                                = connection_cost;
                }
        }

        return instance;
}

double eval_multipliers (instance_t * instance,
                         const double * multipliers,
                         double * OUT_subgradient)
{
        size_t ncities = instance->ncities;
        double * flow = OUT_subgradient;
        double value = 0;
        
        for (size_t i = 0; i < ncities; i++) {
                flow[i] = -1.0;
                value -= multipliers[i];
        }

        for (size_t i = 0; i < instance->nfacilities; i++) {
                facility_t * facility = instance->facilities[i];

                double cur_value = facility->opening_cost;
                for (size_t j = 0; j < ncities; j++) {
                        double delta = facility->arc_cost[j]
                                + multipliers[j];
                        if (delta < 0) cur_value += delta;
                }

                if (cur_value < 0) {
                        value += cur_value;
                        for (size_t j = 0; j < ncities; j++) {
                                double delta = facility->arc_cost[j]
                                        + multipliers[j];
                                if (delta < 0)
                                        flow[j] += 1.0;
                        }
                }
        }

        printf(" city 0: %f/%f ", multipliers[0], flow[0]);

        return value;
}

double upper_bound_multipliers (instance_t * instance, double * multipliers)
{
        size_t ncities = instance->ncities;
        size_t nfacilities = instance->nfacilities;
        double value = 0;
        int * facility_state
                = assert_(calloc(nfacilities, sizeof(int)));

        for (size_t i = 0; i < nfacilities; i++) {
                facility_t * facility = instance->facilities[i];

                double cur_value = facility->opening_cost;
                for (size_t j = 0; j < ncities; j++) {
                        double delta = facility->arc_cost[j]
                                + multipliers[j];
                        if (delta < 0) cur_value += delta;
                }
                if (cur_value < 0) facility_state[i] = 1;
        }

        for (size_t i = 0; i < ncities; i++) {
                double best_cost     = 1.0/0.0;
                size_t best_facility = -1UL;

                double best_closed_cost     = 1.0/0.0;
                size_t best_closed_facility = -1UL;

                for (size_t j = 0; j < nfacilities; j++) {
                        double cost = instance->facilities[j]->arc_cost[i];

                        if (cost < best_cost) {
                                best_closed_cost     = cost;
                                best_closed_facility = j;
                                if (!facility_state[j]) continue;
                                best_cost     = cost;
                                best_facility = j;
                        }
                }

                assert(best_closed_facility != -1UL);

                if (best_facility != -1UL) {
                        facility_state[best_facility] = 2;
                        value += best_cost;
                } else {
                        facility_state[best_closed_facility] = 2;
                        value += best_closed_cost;
                }
        }

        for (size_t i = 0; i < nfacilities; i++)
                if (2 == facility_state[i]) {
                        value += instance->facilities[i]->opening_cost;
                }

        free(facility_state);

        return value;
}

void solve_splp (bundle_solver_t * solver, instance_t * splp)
{
        size_t nmultipliers = splp->ncities*splp->nfacilities;
        struct bundle_sparse_vector svec;
        double * vec = assert_(calloc(nmultipliers, sizeof(double)));

        bundle_set_lambda(solver, vec);

        bundle_start_solving(solver);
        while (!bundle_solve_done_p(solver)) {
                bundle_get_current_lambda(solver, &svec);
                double * OUT_subg = bundle_get_gi_dest(solver, NULL);
                
                if (svec.indices) {
                        for (size_t i = 0; i < nmultipliers; i++)
                                vec[i] = -svec.values[i];
                } else {
                        memset(vec, 0, nmultipliers*sizeof(double));
                        for (size_t i = 0; i < svec.count; i++)
                                vec[svec.indices[i]] = -svec.values[i];
                }
                
                double value = eval_multipliers(splp, vec, OUT_subg);
                bundle_register_values(solver, -value,
                                       OUT_subg, NULL, eps);

                printf(" value %f\n", value);
        }

        bundle_write_solution(solver, vec);
        for (size_t i = 0; i < nmultipliers; i++)
                vec[i] *= -1;

        printf("\nbest value: %f (%f)\n", 
               -bundle_get_best_Fi(solver),
               upper_bound_multipliers(splp, vec));

        free(vec);
}

int main (int argc, char** argv)
{
        bundle_solver_t * solver;
        instance_t instance;
        {
                char * filepath = "MO1";
                if (argc > 1)
                        filepath = argv[1];
                
                FILE * in = assert_(fopen(filepath, "r"));
                instance = read_problem(in);
                assert_(!fclose(in));

                printf("constraints: %lu\n", instance.ncities);

                struct bundle_config config;
                bundle_config_init(&config);
                solver = bundle_solver_create_from_config(&config,
                                                          instance.ncities*instance.nfacilities);
        }

        {
                struct rusage usage;
                assert_(!getrusage(RUSAGE_SELF, &usage));
                double start = usage.ru_utime.tv_sec 
                        + usage.ru_utime.tv_usec * 1e-6;
                solve_splp(solver, &instance);
                assert_(!getrusage(RUSAGE_SELF, &usage));
                double end = usage.ru_utime.tv_sec 
                        + usage.ru_utime.tv_usec * 1e-6;
                printf("Time: %f\n", end-start);
        }

        bundle_solver_destroy(solver);

        return 0;
}

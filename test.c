#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "bundle_solver.h"

typedef struct fun_state {
        struct bundle_callbacks callbacks;
        double * Gi;
        size_t nv;
        double center;
} fun_state_t;

void SetGi(struct bundle_callbacks * callbacks, double * subg, unsigned name)
{
        (void)name;
        ((fun_state_t*)callbacks)->Gi = subg;
}

double Fi(struct bundle_callbacks * callbacks, 
          const double * Lam_val, const unsigned * Lam_idx,
          unsigned Lam_count)
{
        (void)Lam_count;

        double * Gi     = ((fun_state_t*)callbacks)->Gi;
        size_t   nv     = ((fun_state_t*)callbacks)->nv;
        double   center = ((fun_state_t*)callbacks)->center;

        double z = 0;

        if (Lam_idx){
                for (size_t i = 0; i < nv; i++)
                        Gi[i] = 0.0;
                size_t ptr = 0;
                for (size_t i = 0; i < nv; i++) {
                        double x;
                        if (i == Lam_idx[ptr])
                                x = Lam_val[ptr++];
                        else 
                                x = 0;
                        Gi[i] = (x - center);
                        z += Gi[i]*Gi[i];
                }
        } else {
                printf("x: ");
                for (size_t i = 0; i < nv; i++) {
                        Gi[i] = (Lam_val[i] - center);
                        z += Gi[i]*Gi[i];
                        printf("%f ", Lam_val[i]);
                }
                printf("\n");
        }

        return z/2+1;
}

int GetGi (struct bundle_callbacks * callbacks, double * Eps, const unsigned ** Indices)
{
        (void)callbacks;

        *Eps = 0;
        *Indices = 0;

        return 0;
}

int main ()
{
        const size_t nv = 10;
        const double center = 2.0;
        struct bundle_config config;
        bundle_config_init(&config);
        config.iteration_count = 10;
        bundle_solver_t * bundle
                = bundle_solver_create_from_config(&config, nv);
        double * lambda = calloc(nv, sizeof(double));

        srandom(42);

        for (size_t i = 0; i < nv; i++)
                lambda[i] = 100.0 - random()*(200.0/RAND_MAX);
        bundle_set_lambda(bundle, lambda);

        fun_state_t state;
        state.callbacks.SetGi = SetGi;
        state.callbacks.Fi = Fi;
        state.callbacks.GetGi = GetGi;
        state.callbacks.EveryIteration = 0;
        state.nv = nv;
        state.center = center;

        assert(kOK == bundle_solve_with_callbacks(bundle, (struct bundle_callbacks*)&state));

        return 0;

        for (bundle_start_solving(bundle); !bundle_solve_done_p(bundle);) {
                struct bundle_sparse_vector lambda;
                bundle_get_current_lambda(bundle, &lambda);
                double * Gi = bundle_get_gi_dest(bundle, NULL);
                double z = 0;
                
                printf("iteration\n");

                if (lambda.indices) {
                        for (size_t i = 0; i < nv; i++)
                                Gi[i] = 0.0;
                        size_t ptr = 0;
                        for (size_t i = 0; i < nv; i++) {
                                double x;
                                if (i == lambda.indices[ptr])
                                        x = lambda.values[ptr++];
                                else 
                                        x = 0;
                                Gi[i] = (x - center);
                                z += Gi[i]*Gi[i];
                        }
                } else {
                        printf("x: ");
                        for (size_t i = 0; i < nv; i++) {
                                Gi[i] = (lambda.values[i] - center);
                                z += Gi[i]*Gi[i];
                                printf("%f ", lambda.values[i]);
                        }
                        printf("\n");
                }
         
                printf("z : %f\n", z/2+1);
       
                bundle_register_values(bundle, z/2+1,
                                       Gi, NULL, 0);
        }

        bundle_solver_destroy(bundle);
        return 0;
}

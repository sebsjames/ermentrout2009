module;
/*
 * 2D Ermentrout system deriving from RD_Base.
 */
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>

export module rd_erm;

import sm.hexgrid;
import mplot.readcurves;
import sm.hdfdata;
import sm.hexgrid.hdf;

import rd_base;

/*!
 * Reaction diffusion system; Ermentrout 2009.
 */
export namespace sm
{
    template <class F>
    class rd_erm : public sm::rd_base<F>
    {
    public:
        //! Set N>1 for maintaining multiple expression gradients
        alignas(alignof(unsigned int)) unsigned int N = 1;
        //! The c_i(x,t) variables from the Ermentrout paper (chemoattractant concentration)
        alignas(alignof(std::vector<std::vector<F> >)) std::vector<std::vector<F> > c;
        //! The n_i(x,t) variables from the Ermentrout paper (density of tc axons)
        alignas(alignof(std::vector<std::vector<F> >)) std::vector<std::vector<F> > n;
        //! Holds the Laplacian
        alignas(alignof(std::vector<std::vector<F> >)) std::vector<std::vector<F> > lapl;
        //! Holds the Poisson terms (final non-linear term in Ermentrout equation 1)
        alignas(alignof(std::vector<std::vector<F> >)) std::vector<std::vector<F> > poiss;
        //! Sum of c
        alignas(alignof(std::vector<F>)) std::vector<F> sum_c;
        //! Sum of n
        alignas(alignof(std::vector<F>)) std::vector<F> sum_n;

        // Parameters of the Ermentrout model - default values.
        //! Diffusion constant for n
        alignas(F) F Dn = 0.3;
        //! Diffusion constant for c
        alignas(F) F Dc = Dn * 0.3;
        //! saturation term in function for production of c
        alignas(F) F beta = 5.0;
        //! production of new axon branches
        alignas(F) F a = 1.0;
        //! pruning constant
        alignas(F) F b = 1.0;
        //! decay of chemoattractant constant
        alignas(F) F mu = 1.0;
        //! degree of attraction of chemoattractant
        alignas(F) F chi = Dn;

        //! Frame number, used when saving PNG movie frames.
        unsigned int frameN = 0;

        //! Simple constructor; no arguments.
        RD_Erm() : morph::RD_Base<F>() {}

        //! Perform memory allocations, vector resizes and so on.
        virtual void allocate() { morph::RD_Base<F>::allocate(); }

        //! Initialise HexGrid, variables. Carry out any one-time computations of the model.
        virtual void init()
        {
            // Resize and zero-initialise the various containers
            this->resize_vector_vector (this->c, this->N);
            this->resize_vector_vector (this->n, this->N);
            this->resize_vector_vector (this->lapl, this->N);
            this->resize_vector_vector (this->poiss, this->N);
            this->sum_n.resize (this->N, F{0});
            this->sum_c.resize (this->N, F{0});

            // Initialise with noise
            F noiseamp = F{0.1};
            for (unsigned int i  = 0; i < this->N; ++i) {
                this->noiseify_vector_variable (this->n[i], 1., noiseamp);
                this->noiseify_vector_variable (this->c[i], beta*0.5, noiseamp);
            }
        }

        //! Compute one step of the model
        virtual void step()
        {
            this->stepCount++;

            for (unsigned int i=0; i<this->N; ++i) {

                sum_n[i] = F{0};
                sum_c[i] = F{0};

                this->compute_poiss (n[i],c[i],i);  // compute the non-linear Poission term in Eq1
                this->compute_lapl (n[i], i);       // populate lapl[i] with laplacian of n

                // integrate n
                for (unsigned int h=0; h<this->nhex; ++h) {
                    n[i][h] += (a - b*n[i][h] + Dn*lapl[i][h] - chi*poiss[i][h])*this->dt;
                    sum_n[i] += n[i][h];
                }

                this->compute_lapl (c[i], i);       // populate lapl[i] with laplacian of c

                // integrate c
                F n2;
                for (unsigned int h=0; h<this->nhex; ++h) {
                    n2 = n[i][h]*n[i][h];
                    c[i][h] += (beta*n2/(F{1}+n2) - mu*c[i][h] +Dc*lapl[i][h])*this->dt;
                    sum_c[i] += c[i][h];
                }
            }
        }

        //! Computes the Laplacian Stable with dt = 0.0001;
        void compute_lapl (std::vector<F>& fa, unsigned int i)
        {
#pragma omp parallel for schedule(static) shared(lapl, fa, i)
            for (unsigned int hi=0; hi<this->nhex; ++hi) {

                // 1. The D Del^2 term

                // Compute the sum around the neighbours
                F thesum = F{-6} * fa[hi];
                if (this->hg->has_ne (hi)) {
                    thesum += fa[this->hg->ne (hi)];
                } else {
                    thesum += fa[hi]; // A ghost neighbour-east with same value as Hex_0
                }
                if (this->hg->has_nne (hi)) {
                    thesum += fa[this->hg->nne (hi)];
                } else {
                    thesum += fa[hi];
                }
                if (this->hg->has_nnw (hi)) {
                    thesum += fa[this->hg->nnw (hi)];
                } else {
                    thesum += fa[hi];
                }
                if (this->hg->has_nw (hi)) {
                    thesum += fa[this->hg->nw (hi)];
                } else {
                    thesum += fa[hi];
                }
                if (this->hg->has_nsw (hi)) {
                    thesum += fa[this->hg->nsw (hi)];
                } else {
                    thesum += fa[hi];
                }
                if (this->hg->has_nse (hi)) {
                    thesum += fa[this->hg->nse (hi)];
                } else {
                    thesum += fa[hi];
                }

                this->lapl[i][hi] = this->twoover3dd * thesum;
            }
        }

        //! Computes the Poisson term. Stable with dt = 0.0001;
        void compute_poiss (std::vector<F>& fa1, std::vector<F>& fa2, unsigned int i)
        {
            // Compute non-linear term
#pragma omp parallel for schedule(static) shared(poiss, fa1, fa2)
            for (unsigned int hi=0; hi<this->nhex; ++hi) {

                std::vector<F> dum1(6,fa1[hi]);
                std::vector<F> dum2(6,fa2[hi]);

                if (this->hg->has_ne (hi)) {
                    dum1[0] = fa1[this->hg->ne (hi)];
                    dum2[0] = fa2[this->hg->ne (hi)];
                }
                if (this->hg->has_nne (hi)) {
                    dum1[1] = fa1[this->hg->nne (hi)];
                    dum2[1] = fa2[this->hg->nne (hi)];
                }
                if (this->hg->has_nnw (hi)) {
                    dum1[2] = fa1[this->hg->nnw (hi)];
                    dum2[2] = fa2[this->hg->nnw (hi)];
                }
                if (this->hg->has_nw (hi)) {
                    dum1[3] = fa1[this->hg->nw (hi)];
                    dum2[3] = fa2[this->hg->nw (hi)];
                }
                if (this->hg->has_nsw (hi)) {
                    dum1[4] = fa1[this->hg->nsw (hi)];
                    dum2[4] = fa2[this->hg->nsw (hi)];
                }
                if (this->hg->has_nse (hi)) {
                    dum1[5] = fa1[this->hg->nse (hi)];
                    dum2[5] = fa2[this->hg->nse (hi)];
                }

                // John Brooke's final thesis solution (based on 'finite volume method'
                // of Lee et al. https://doi.org/10.1080/00207160.2013.864392
                F val = (dum1[0]+fa1[hi]) * (dum2[0]-fa2[hi])
                + (dum1[1]+fa1[hi]) * (dum2[1]-fa2[hi])
                + (dum1[2]+fa1[hi]) * (dum2[2]-fa2[hi])
                + (dum1[3]+fa1[hi]) * (dum2[3]-fa2[hi])
                + (dum1[4]+fa1[hi]) * (dum2[4]-fa2[hi])
                + (dum1[5]+fa1[hi]) * (dum2[5]-fa2[hi]);
                this->poiss[i][hi] = val * this->oneover3dd;// / (3 * this->d * this->d);
            }
        }

        //! Save state to HDF5
        void saveState()
        {
            std::string fname = this->logpath + "/2Derm.h5";
            morph::HdfData data (fname);

            // Save some variables
            for (unsigned int i = 0; i<this->N; ++i) {
                std::stringstream vss;
                vss << "c_" << i;
                std::string vname = vss.str();
                data.add_contained_vals (vname.c_str(), this->c[i]);
                vname[0] = 'n';
                data.add_contained_vals (vname.c_str(), this->n[i]);
            }

            // Parameters
            data.add_val ("/Dn", this->Dn);
            data.add_val ("/Dc", this->Dc);
            data.add_val ("/beta", this->beta);
            data.add_val ("/a", this->a);
            data.add_val ("/b", this->b);
            data.add_val ("/mu", this->mu);
            data.add_val ("/chi", this->chi);

            // HexGrid information
            this->saveHexPositions (data);
        }

    }; // RD_Erm
} // namespace

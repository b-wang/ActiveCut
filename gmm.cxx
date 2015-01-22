#include <common.h>
#include <gmm.h>
#include <utility.h>
// one thing I worry is the initial user input. Within the box, all data are
// init'd to forground, while some are not. This may result in inaccurate
// foreground gmm parameter estimaties.

// gmm parameter estimation. (M step).
int gmm_mstep(const vnl_matrix <double> & data, 
	      const vnl_matrix<double> & gmm_labels,
	      const vnl_vector<unsigned> & alpha,
	      unsigned whatground,
	      GMMType & gmm,
	      bool cov_shared,
	      unsigned cov_type) 

{
     unsigned label = 0;
     unsigned n_samples = data.rows();
     unsigned n_channels = data.cols();
     vnl_vector<double> sample_centered;

     // clear to zero before computation. 
     gmm.n_pts = 0;
     for (unsigned comp_id = 0; comp_id < gmm.n_comp; comp_id ++) {
	  gmm.comp[comp_id].label = 0;
	  gmm.comp[comp_id].numPts = 0;
	  gmm.comp[comp_id].mu.fill( 0 );
	  gmm.comp[comp_id].cov.fill( 0 );
     }

     // estimate mu. First compute sum. 
     for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	  // assert (gmm_labels.get_row(sample_id).sum() == 1);
	  if (alpha[sample_id] == whatground) {
	       gmm.n_pts ++;
	       for (unsigned k = 0; k < gmm.n_comp; k ++) {
		    gmm.comp[k].numPts += gmm_labels[sample_id][k];
		    gmm.comp[k].mu +=  gmm_labels[sample_id][k] * data.get_row(sample_id);
	       }
	  }
     } 

     // comptue mean from the sum.
     for (unsigned k = 0; k < gmm.n_comp; k ++) {
     	  gmm.comp[k].mu /= gmm.comp[k].numPts;
     }

     // compute covariance matrix.
     if (cov_shared) {
	  vnl_matrix<double> big_cov(n_channels, n_channels, 0);
	  for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	       if (alpha[sample_id] == whatground) {
		    for (unsigned k = 0; k < gmm.n_comp; k ++) {
			 sample_centered = data.get_row(sample_id) - gmm.comp[k].mu;
			 big_cov += gmm_labels[sample_id][k] * outer_product(sample_centered, sample_centered);
		    } // k
	       } // alpha
	  } // sample_id
	  big_cov /= gmm.n_pts;

	  for (unsigned k = 0; k < gmm.n_comp; k ++) {
	       gmm.comp[k].cov = big_cov;
	  }
     }

     else {	  // each comp has own cov mat.
	  for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	       if (alpha[sample_id] == whatground) {
		    for (unsigned k = 0; k < gmm.n_comp; k ++) {
			 sample_centered = data.get_row(sample_id) - gmm.comp[k].mu;
			 gmm.comp[k].cov += gmm_labels[sample_id][k] * outer_product(sample_centered, sample_centered);
		    }
	       }
	  }

	  // normalize
	  for (unsigned k = 0; k < gmm.n_comp; k ++) {
	       gmm.comp[k].cov /= gmm.comp[k].numPts;
	  }
     }

     // regularize the covariance matrix. 
     if (cov_type == 0 ) {
	  // full covariance matrix. No regularization.
     }
     if (cov_type >= 1) {
	  // diagonal matrix. Set non-diagonal to zero.
	  for (unsigned k = 0; k < gmm.n_comp; k ++) {
	       for (unsigned i = 0; i < n_channels; i ++) {
		    for (unsigned j = 0; j < n_channels; j++) {
			 if (i != j) {
			      gmm.comp[k].cov(i,j) = 0;
			 }
		    }
	       }
	  } // k
     }
     if (cov_type >= 3) {
	  // identity matrix.
	  for (unsigned k = 0; k < gmm.n_comp; k ++) {
	       gmm.comp[k].cov.set_identity();
	  }
     }

     // Compute inverse of cov.
     for (unsigned k = 0; k < gmm.n_comp; k ++) {
	  // check singularity.
	  vnl_diag_matrix<double> diag_filler(n_channels, EPS);
	  // gmm.comp[k].cov.print(std::cout);
	  while (vnl_rank(gmm.comp[k].cov) < n_channels) {
	       gmm.comp[k].cov += diag_filler;
	       printf("gmm_mstep(): inv_cov[%i](FG): add on diag: %f\n", k, diag_filler[0]);
	       diag_filler *= 2;
	  }
	  gmm.comp[k].inv_cov = vnl_matrix_inverse<double>(gmm.comp[k].cov);	  
	  gmm.comp[k].det_cov = vnl_determinant(gmm.comp[k].cov);
     }

     // update pi. (do not need this code since we have a function for updating
     // the random variable pi given the hidden variables z.

     // for (unsigned k = 0; k < gmm.n_comp; k ++) {
     // 	  gmm.pi[k] = gmm.comp[k].numPts;
     // }
     // gmm.pi /= gmm.pi.sum();

     return 0;
}

// estimate the label posterior, given the gmm parameters. (E step).
double gmm_estep(const vnl_matrix <double> & data, 
		 vnl_sparse_matrix <double> & con_map,
		 vnl_matrix<double> & gmm_labels,
		 const vnl_vector<unsigned> & alpha,
		 unsigned whatground,
		 const GMMType & gmm)
{
     unsigned n_samples = gmm_labels.rows();
     unsigned n_channels = data.cols();
     unsigned n_comp = gmm.n_comp;
     double m = 0;
     // vnl_vector <double> zero_vec(gmm_labels.cols(), 0);

     vnl_vector<double> sample_c; // centered.
     vnl_vector<double> exp_term (gmm.n_comp, 0);
     vnl_matrix<double> gmm_oldlabels(gmm_labels);
     unsigned nbr_id = 0;
     vnl_sparse_matrix<double>::row con_map_row;

#pragma omp parallel for schedule(dynamic, 1000)  default(none) private(con_map_row, sample_c, exp_term, nbr_id, m) shared(alpha, whatground, gmm, con_map, n_channels, gmm_labels, data, n_samples, n_comp)
     for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {     
	  exp_term.set_size(gmm.n_comp);
	  exp_term.fill(0);
     	  if (alpha[sample_id] == whatground) {
	       con_map_row = con_map.get_row(sample_id);		    
	       for (unsigned k = 0; k < gmm.n_comp; k ++) {
		    sample_c = data.get_row(sample_id) - gmm.comp[k].mu;
		    exp_term[k] = (- 0.5 * n_channels) * log(2 * PI)
			 - 0.5 * log(gmm.comp[k].det_cov) 
			 - 0.5 * dot_product(gmm.comp[k].inv_cov * sample_c, sample_c)
			 + log(gmm.pi(sample_id, k));

		    vnl_vector<double> z(gmm_labels.cols(), 0);
		    z[k] = 1;
		    for (vnl_sparse_matrix<double>::row::const_iterator col_iter = con_map_row.begin(); col_iter != con_map_row.end();  ++ col_iter) {
			 // go the neighbor sample_id
			 nbr_id = (*col_iter).first;
			 if (alpha[nbr_id] == whatground) { // in same ground.
			      exp_term[k] += gmm.beta * dot_product(z, gmm_labels.get_row(nbr_id));
			      // the above line does not take into account the
			      // patches. Int assume each patch is just a single voxel. To
			      // account for that, need to follow the build_nliks func.
			 }
		    } // con_iter
	       } // k
	       m = exp_term.max_value();
	       exp_term = exp_term - m;

	       for (unsigned k = 0; k < gmm.n_comp; k ++) {
		    gmm_labels(sample_id, k) = exp(exp_term[k]);
	       }	       
	       // normalize so it sum to 1. (unfilled values are just zero so it
	       // doesn't matter) Since each row may be for FG or BG. the col
	       // num may be larger than necessary. make sure unused elements
	       // does not contribute to the sum.
	       gmm_labels.scale_row(sample_id, 1/ gmm_labels.get_row(sample_id).extract(n_comp, 0).sum() );
	  } // alpha
     } // sample_id

     // sum of absolute values / 2 == the change of labels. 
     return (gmm_labels - gmm_oldlabels).array_one_norm() * 0.5;
}

int update_pi(const vnl_vector<unsigned> & alpha,
	      const vnl_matrix<double> & atlas,
	      GMMType & gmm,
	      double lambda)
{
     for (unsigned n = 0; n < atlas.rows(); n ++) {
	  // Important: For points not in this gmm, we still need a pi value
	  // (will be usedin the graphcuts and gmm_eval_ll func). So, no check
	  // for whatground.
	  for (unsigned k = 0; k < gmm.n_comp; k ++) {
	       gmm.pi(n, k) = pow(atlas(n, k) + EPS, lambda) * double(gmm.comp[k].numPts)/double(gmm.n_pts);
	  }

	  // normalize so pi sums to one. 
	  gmm.pi.scale_row(n, 1/gmm.pi.get_row(n).sum() );
     } // n
     return 0;
}

// query the log-likelihood of one data point. 
double gmm_eval_ll(const vnl_vector<double> & data,
		   const GMMType & gmm,
		   unsigned sample_id)
{
     vnl_vector<double> sample_c; // centered.
     vnl_vector<double> exp_term (gmm.n_comp, 0);
     unsigned n_channels = gmm.comp[0].mu.size();
     double sum_term = 0, m = 0;

     for (unsigned k = 0; k < gmm.n_comp; k ++) {
	  sample_c = data - gmm.comp[k].mu;
	  exp_term[k] = (- 0.5 * n_channels) * log(2 * PI)
	       - 0.5 * log(gmm.comp[k].det_cov) 
	       - 0.5 * dot_product(gmm.comp[k].inv_cov * sample_c, sample_c);
     }
     m = exp_term.max_value();
     exp_term = exp_term - m; // to avoid underflow.
     sum_term = 0;
     for (unsigned k = 0; k < gmm.n_comp; k ++) {
	  sum_term += gmm.pi(sample_id, k) * exp(exp_term[k]);
     }	       
     
     return log(sum_term) + m;
}

int kmeans(const vnl_matrix <double> & data,
	   const vnl_vector<unsigned> & alpha,
	   vnl_matrix<double> & gmm_labels, // initial gmm labels.
	   ParType & par,
	   unsigned whatground)
{
     unsigned n_comp = 0;
     if (whatground == ALPHA_FG) n_comp = par.gmm_fg.n_comp;
     else if (whatground == ALPHA_BG) n_comp = par.gmm_bg.n_comp;
     else {
	  printf("whatground must be either FG or BG!.\n");
	  exit (1);
     }

     // cluster centers. 
     std::vector<vnl_matrix<double> > cc(par.kmeansruns);
     vnl_vector<double> mean_ssr(par.kmeansruns, 0);
     double mean_ssr_old = 0;

     // init by kmeans++
#pragma omp parallel for default(none) private(mean_ssr_old) shared(mean_ssr, data, alpha, par, whatground, gmm_labels, n_comp, cc)

     for (unsigned r = 0; r < par.kmeansruns; r ++) {
	  printf("kmeans run: %i begin\n", r);
	  mean_ssr_old = 0;
	  printf("n_comp = %i, par.n_channels = %i\n", n_comp, par.n_channels);
	  (cc[r]).set_size(n_comp, par.n_channels);
	  (cc[r]).fill(0);
	  kmeans_init(data, alpha, cc[r], par, whatground, par.seed + r);

	  do {
	       // update labels. 
	       mean_ssr_old = mean_ssr[r];
	       mean_ssr[r] = kmeans_updatelabels(data, cc[r], alpha, gmm_labels, par, whatground);

	       // update cluster centers.
	       kmeans_updatecc(data, alpha, cc[r], gmm_labels, par, whatground);
	       // mean_ssr = compute_mean_ssr(data, alpha, cc, gmm_labels, par, whatground);

	       if (par.verbose >= 2) 
		    printf("   kmeans run %i,  mean_ssr = %E\n", r, mean_ssr[r]);
	  }
	  while(fabs((mean_ssr[r] - mean_ssr_old) / mean_ssr_old) > 1e-5);
	  
	  // if (best_mean_ssr > mean_ssr) {
	  //      best_mean_ssr = mean_ssr;
	  //      best_cc = cc;
	  // }
	  // if (par.verbose >= 1) 
	  //      printf("  kmeans run %i done. mean_ssr = %E. best_mean_ssr = %E\n", r, mean_ssr[r], best_mean_ssr);
	   printf("  kmeans run %i done. mean_ssr = %E.\n", r, mean_ssr[r]);
     } // kmeans run r.

     // give best cc, re-estimate lables. 
     mean_ssr[0] = kmeans_updatelabels(data, cc[mean_ssr.arg_max()], alpha, gmm_labels, par, whatground);
     printf("kmeans(): best_mean_ssr = %E\n", mean_ssr[0]);
}

int kmeans_init(const vnl_matrix <double> & data,
		const vnl_vector<unsigned> & alpha,
		vnl_matrix<double> & cc,
		ParType & par,
		unsigned whatground,
		unsigned seed)
{
     boost::random::mt19937 rng(42u);         // produces randomness out of thin air
     rng.seed(static_cast<unsigned int>(seed));
     boost::random::uniform_real_distribution<> uni_dist(0, 1);
     unsigned n_comp = cc.rows();

     // define a pdf. it include both fg and bg samples. But unused cell will be
     // zero and have no effect on the results.
     vnl_vector<double> pdf(par.n_samples, 0);
     
     // find the first bogus center. 
     unsigned sample_id = 0;
     uni_dist(rng);

     do {
     	  // find one in FG or BG. 
     	  sample_id =  (int)(floor) (uni_dist(rng) * par.n_samples);
     }
     while (alpha(sample_id ) != whatground);
     cc.set_row(0, data.get_row(sample_id) );

     unsigned clsIdx = 0;
     double rand_num = 0;
     double cdf = 0;
     for (clsIdx = 0; clsIdx < n_comp; clsIdx ++) {
     	  compute_dist(pdf, data, alpha, cc, clsIdx, whatground);
     	  rand_num = uni_dist(rng);
     	  cdf = 0;

     	  sample_id = 0;
     	  while (rand_num > cdf) {
     	       sample_id ++;
     	       cdf += pdf[sample_id];
     	  }

     	  // found the point. 
	  if (par.verbose >= 2) {
	       printf("kmeans_init(): cc[%i] init to sample %i.\n", clsIdx, sample_id);
	  }
     	  cc.set_row(clsIdx, data.get_row(sample_id));
     }

     return 0;
}

int compute_dist(vnl_vector<double> & pdf,
		 const vnl_matrix<double> data, 
		 const vnl_vector<unsigned> & alpha,
		 const vnl_matrix<double> & cc,
		 unsigned cur_comp_id, // current cluster id.
		 unsigned whatground)
{
     double min_dist = 0, this_dist = 0;
     unsigned prev_comp_id = 0;
     unsigned n_samples = data.rows();
     pdf.fill( 0 );

     for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	  if (alpha[sample_id] == whatground) {
	       min_dist = 1e10;
	       if (cur_comp_id == 0) {
		    min_dist = (data.get_row(sample_id) - cc.get_row(0)).squared_magnitude();
	       }
	       else {
		    for (prev_comp_id = 0; prev_comp_id < cur_comp_id; prev_comp_id ++) {
			 this_dist = (data.get_row(sample_id) - cc.get_row(prev_comp_id)).squared_magnitude();
			 if (this_dist < min_dist)
			      min_dist = this_dist;
		    }
	       }

	       // now I have the min_dist for current point. save it. 
	       pdf[sample_id] = min_dist;
	  }  // in whatground.
     }

     // normalize the pdf. the zero entry outside of whatground does not have
     // effect of summation.
     pdf = pdf / pdf.sum();
     return 0;
}

double kmeans_updatelabels(const vnl_matrix <double> & data,
			   const vnl_matrix<double> & cc,
			   const vnl_vector<unsigned> & alpha,
			   vnl_matrix<double> & gmm_labels, // initial gmm labels.
			   ParType & par,
			   unsigned whatground)
{
     double mean_ssr = 0;
     double nearest_dist = 0, this_dist = 0;
     unsigned nearest_label = 0;
     unsigned n_points = 0, n_comp = 0;

     if (whatground == ALPHA_FG) n_comp = par.gmm_fg.n_comp;
     else if (whatground == ALPHA_BG) n_comp = par.gmm_bg.n_comp;

     vnl_vector <double> zero_vec(gmm_labels.cols(), 0);

     for (unsigned sample_id = 0; sample_id < par.n_samples; sample_id ++) {
	  if (alpha[sample_id] == whatground) {
	       n_points ++;
	       nearest_dist = 1e10;
	       for (unsigned clsIdx = 0; clsIdx < n_comp; clsIdx ++) {
		    // printf("squared_mag: %f    ", squared_mag(data[sample_id], cc[clsIdx], par.n_channels));
		    // printf("orig: %f\n", (data.get_row(sample_id) - cc.get_row(clsIdx)).squared_magnitude());
		    // this_dist  = (data.get_row(sample_id) - cc.get_row(clsIdx)).squared_magnitude();
		    this_dist = squared_mag(data[sample_id], cc[clsIdx], par.n_channels);
		    if (this_dist < nearest_dist ) {
			 nearest_dist = this_dist;
			 nearest_label = clsIdx;
		    }
	       } // clsIdx
	       
	       // found the nearest cluster center. 
	       gmm_labels.set_row(sample_id, zero_vec);
	       gmm_labels(sample_id, nearest_label) = 1;
	       mean_ssr += nearest_dist;

	       // assert (gmm_labels.get_row(sample_id).sum() == 1);
	  } // whatground
     } // for
     
     mean_ssr = mean_ssr / n_points;
     return mean_ssr;
}

double squared_mag(const double * A, const double * B, unsigned len)
{
     double r = 0;
     for (unsigned n = 0; n < len; n ++) {
	  r += pow((A[n] - B[n]), 2);
     }
     
     return r;
}
// update cluster center. 
int kmeans_updatecc(const vnl_matrix <double> & data,
		    const vnl_vector<unsigned> & alpha,
		    vnl_matrix<double> & cc,
		    const vnl_matrix<double> & gmm_labels, // initial gmm labels.
		    ParType & par,
		    unsigned whatground)
{
     unsigned clsIdx = 0;
     unsigned n_comp = cc.rows();
     vnl_vector<double> n_points(n_comp, 0);
     vnl_matrix<double> cc_sum(n_comp, par.n_channels, 0);

     for (unsigned sample_id = 0; sample_id < par.n_samples; sample_id ++) {
	  if (alpha[sample_id] == whatground) {     
	       clsIdx = gmm_labels.get_row(sample_id).arg_max();
	       cc_sum.set_row(clsIdx, cc_sum.get_row(clsIdx) + data.get_row(sample_id) );
	       n_points[clsIdx] ++;
	  } // whatground
     } // for

     // compute mean from the sum.
     for (clsIdx = 0; clsIdx < n_comp; clsIdx ++) {
	  cc.set_row(clsIdx, cc_sum.get_row(clsIdx) / n_points[clsIdx]);
     }
     return 0;
}

double compute_mean_ssr(const vnl_matrix <double> & data,
			const vnl_vector<unsigned> & alpha,
			vnl_matrix<double> & cc,
			const vnl_matrix<double> & gmm_labels, // initial gmm labels.
			ParType & par,
			unsigned whatground)
{
     double mean_ssr = 0;
     unsigned n_points = 0;
     unsigned clsIdx = 0;
     for (unsigned sample_id = 0; sample_id < par.n_samples; sample_id ++) {
	  if (alpha[sample_id] == whatground) {
	       n_points ++;
	       clsIdx = gmm_labels.get_row(sample_id).arg_max();
	       // mean_ssr += (data.get_row(sample_id) - cc.get_row(clsIdx)).squared_magnitude();
	       mean_ssr += squared_mag(data[sample_id], cc[clsIdx], par.n_channels);
	  }
     }
     
     mean_ssr = mean_ssr / n_points;
     return mean_ssr;
}

double eval_ll(const vnl_matrix <double> & data, 
	       const vnl_vector<unsigned> & alpha,
	       const GMMType & gmm,
	       unsigned whatground)
{
     unsigned n_samples = data.rows();
     unsigned n_comp = gmm.n_comp;
     unsigned n_channels = data.cols();

     vnl_vector<double> exp_term (gmm.n_comp, 0);
     vnl_vector<double> sample_c (n_channels);
     
     double m = 0; // max value of the exp term.
     double LL = 0, sum_term = 0;
     for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	  if (alpha[sample_id] == whatground) {     
	       for (unsigned k = 0; k < n_comp; k ++) {
		    sample_c = data.get_row(sample_id) - gmm.comp[k].mu;
		    exp_term[k] = (- 0.5 * n_channels) * log(2 * PI)
			 - 0.5 * log(gmm.comp[k].det_cov) 
			 - 0.5 * dot_product(gmm.comp[k].inv_cov * sample_c, sample_c);
	       }
	       m = exp_term.max_value();
	       exp_term = exp_term - m;
	       sum_term = 0;
	       for (unsigned k = 0; k < n_comp; k ++) {
		    sum_term += gmm.pi(sample_id, k) * exp(exp_term[k]);
	       }	       

	       // single data point LL added to total.
	       LL = LL + log(sum_term) + m;
	  } // whatground
     } // for

     return LL;
}

double eval_ell(const vnl_matrix <double> & data, 
	       const vnl_matrix<double> & gmm_labels, // initial gmm labels.
	       const vnl_vector<unsigned> & alpha,
	       const GMMType & gmm,
	       unsigned whatground)
{
     unsigned n_samples = data.rows();
     unsigned n_comp = gmm.n_comp;
     unsigned n_channels = data.cols();

     double exp_term = 0;
     vnl_vector<double> sample_c (n_channels);
     
     double m = 0; // max value of the exp term.
     double LL = 0, sum_term = 0;
     for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	  if (alpha[sample_id] == whatground) {     
	       for (unsigned k = 0; k < n_comp; k ++) { // 
		    sample_c = data.get_row(sample_id) - gmm.comp[k].mu;
		    exp_term = (- 0.5 * n_channels) * log(2 * PI)
			 - 0.5 * log(gmm.comp[k].det_cov) 
			 - 0.5 * dot_product(gmm.comp[k].inv_cov * sample_c, sample_c);
		    
		    LL = LL + gmm_labels(sample_id, k) * (log(gmm.pi(sample_id, k)) + exp_term);
	       }
	  } // whatground
     } // for

     return LL;
}

double low_bound(const vnl_matrix<double> & data, 
		 const vnl_matrix<double> & gmm_labels,
		 const vnl_vector<unsigned> & alpha,
		 const GMMType & gmm,
		 unsigned whatground)
{
     unsigned n_samples = data.rows();
     unsigned n_comp = gmm.n_comp;
     unsigned n_channels = data.cols();

     // lower bound  = Q(theta, theta_old) + entropy
     double q_func = eval_ell(data, gmm_labels, alpha, gmm, whatground);
     
     // compute entropy.
     double exp_term = 0;
     double m = 0; // max value of the exp term.
     double entropy = 0;
     for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	  if (alpha[sample_id] == whatground) {     
	       for (unsigned k = 0; k < n_comp; k ++) {
		    if (gmm_labels(sample_id, k) > 0) {
			 entropy += gmm_labels(sample_id, k) * log (gmm_labels(sample_id, k));
		    }
		    else {
			 // do nothing, since the 0 * log(0) has limit 0.
		    }
	       }
	  } // whatground
     } // for

     return q_func - entropy;
}


double kl_dvg(const vnl_matrix <double> & data, 
	      vnl_sparse_matrix <double> & con_map,
	      const vnl_matrix<double> & gmm_labels,
	      const vnl_vector<unsigned> & alpha,
	      const GMMType & gmm,
	      unsigned whatground)
{
     unsigned n_samples = alpha.size();
     unsigned n_comp = gmm.n_comp;
     unsigned n_channels = gmm.comp[0].mu.size();

     // compute new p(z|x)
     vnl_matrix<double> new_post(n_samples, n_comp, 0);
     gmm_estep(data, con_map, new_post, alpha, whatground, gmm);

     // compute entropy.
     double exp_term = 0;
     double m = 0; // max value of the exp term.
     double KL = 0;
     for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	  if (alpha[sample_id] == whatground) {     
	       for (unsigned k = 0; k < n_comp; k ++) {
		    if (gmm_labels(sample_id,k) != 0 && new_post(sample_id,k) != 0)
			 KL -= gmm_labels(sample_id, k) * (log(new_post(sample_id, k)) - log(gmm_labels(sample_id, k)) );
		    else if (gmm_labels(sample_id,k) != 0 && new_post(sample_id,k) == 0)
		    {
			 // printf("kl_dvg(): p(z|x) = 0 but q(z) != 0. No idea!\n");
		    }
		    else if (gmm_labels(sample_id,k) == 0 && new_post(sample_id,k) != 0) {
			 // these cases, would be zero.
		    }
		    else if (gmm_labels(sample_id,k) == 0 && new_post(sample_id,k) == 0) { // both are zero.
			 // printf("kl_dvg(): sample %i, comp %i, both p(z|x) and q(z) = 0. \n", sample_id, k);
		    }
	       }
	  } // whatground
     }
     return KL;
}


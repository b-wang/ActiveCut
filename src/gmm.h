int gmm_mstep(const vnl_matrix <double> & data, 
	      const vnl_matrix<double> & gmm_labels, // initial gmm labels.
	      const vnl_vector<unsigned> & alpha,
	      unsigned whatground,
	      GMMType & gmm,
     	      bool cov_shared,
     	      unsigned cov_type);

// estimate the label posterior, given the gmm parameters. (E step).
double gmm_estep(const vnl_matrix <double> & data, 
		 vnl_sparse_matrix <double> & con_map,
		 vnl_matrix<double> & gmm_labels,
		 const vnl_vector<unsigned> & alpha,
		 unsigned whatground,
		 const GMMType & gmm);

// query the log-likelihood of one data point. 
double gmm_eval_ll(const vnl_vector<double> & data,
		   const GMMType & gmm,
		   unsigned sample_id);

int kmeans_init(const vnl_matrix <double> & data,
		const vnl_vector<unsigned> & alpha,
		vnl_matrix<double> & cc,
		ParType & par,
		unsigned whatground,
		unsigned seed);

int kmeans(const vnl_matrix <double> & data,
	   const vnl_vector<unsigned> & alpha,
	   vnl_matrix<double> & gmm_labels, // initial gmm labels.
	   ParType & par,
	   unsigned whatground);

int compute_dist(vnl_vector<double> & pdf,
		 const vnl_matrix<double> data, 
		 const vnl_vector<unsigned> & alpha,
		 const vnl_matrix<double> & cc,
		 unsigned cur_comp_id,
		 unsigned whatground);

double kmeans_updatelabels(const vnl_matrix <double> & data,
			   const vnl_matrix<double> & cc,
			   const vnl_vector<unsigned> & alpha,
			   vnl_matrix<double> & gmm_labels, // initial gmm labels.
			   ParType & par,
			   unsigned whatground);

int kmeans_updatecc(const vnl_matrix <double> & data,
		    const vnl_vector<unsigned> & alpha,
		    vnl_matrix<double> & cc,
		    const vnl_matrix<double> & gmm_labels, // initial gmm labels.
		    ParType & par,
		    unsigned whatground);

double compute_mean_ssr(const vnl_matrix <double> & data,
			const vnl_vector<unsigned> & alpha,
			vnl_matrix<double> & cc,
			const vnl_matrix<double> & gmm_labels, // initial gmm labels.
			ParType & par,
			unsigned whatground);


// evaluate the log log-likelihood of all data points in the given
// Fore/background.  Might be outdated since I aded MRF.
double eval_ll(const vnl_matrix <double> & data, 
	       const vnl_vector<unsigned> & alpha,
	       const GMMType & gmm,
	       unsigned whatground);

// evaluate the expectation of log log-likelihood of all data points in the
// given Fore/background. This is for debugging of EM algorithm. Might be outdated since I aded MRF. 
double eval_ell(const vnl_matrix <double> & data, 
	       const vnl_matrix<double> & gmm_labels, // initial gmm labels.
	       const vnl_vector<unsigned> & alpha,
	       const GMMType & gmm,
		unsigned whatground);


double low_bound(const vnl_matrix<double> & data, 
		 const vnl_matrix<double> & gmm_labels,
		 const vnl_vector<unsigned> & alpha,
		 const GMMType & gmm,
		 unsigned whatground);

double kl_dvg(const vnl_matrix <double> & data, 
	      vnl_sparse_matrix <double> & con_map,
	      const vnl_matrix<double> & gmm_labels,
	      const vnl_vector<unsigned> & alpha,
	      const GMMType & gmm,
	      unsigned whatground);

double mrf_term(const vnl_sparse_matrix<double>::row & con_map,
		const vnl_matrix<double> & gmm_labels,
		const vnl_vector<unsigned> & alpha,
		unsigned whatground,
		unsigned sample_id);

// read the atlas prior, and assign to pi. Need to run for foreground gmm and
// background seperately. Althrough we fix the atlas, the pi depends on both the
// atlas and the current proportion of data points in certain components. So we
// need run this func inside EM loop. (Think about atlas is very low confidence
// with small lambda, pi need to be close to the class proportion, just like the
// regular gmm without the atlas prior).
int update_pi(const vnl_vector<unsigned> & alpha,
	      const vnl_matrix<double> & atlas,
	      GMMType & gmm,
	      double lambda);

double squared_mag(const double * A, const double * B, unsigned len);

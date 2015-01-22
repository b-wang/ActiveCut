double compute_beta(const vnl_matrix <double> & data, 
		    vnl_sparse_matrix <double> & nlinks_map);
int print_par(ParType par);
int print_gmm(const GMMType & gmm);
int print_vnl_matrix(const vnl_matrix<double> & mat, unsigned r, unsigned c);
int print_vnl_vec(const vnl_vector<double> & vec, unsigned N);

int save_gmm_labelmap(vnl_matrix<double> & gmm_labels, // initial gmm labels.
		      ParType & par,
		      const vnl_vector<unsigned> & alpha,
		      unsigned whatground,
		      ImageType3DU::Pointer lindexPtr,
		      std::string filename);

// save gmm posterior into 4D file.
int save_gmm_posterior(vnl_matrix<double> & gmm_labels, // initial gmm labels.
		      ParType & par,
		      const vnl_vector<unsigned> & alpha,
		      unsigned whatground,
		      ImageType3DU::Pointer lindexPtr,
		       std::string filename);

int save_alpha(const ParType & par,
	       const vnl_vector<unsigned> & alpha,
	       ImageType3DU::Pointer lindexPtr,
	       std::string filename);

// for debug purpose. Save NxP data matrix back into image. The voxels in same
// patch of the output image will have same intensity, which is from the mean of
// the intensities of the patch.
int save_patches(ImageType3DI::Pointer parcelPtr,
	     const vnl_matrix <double> & data,
	     std::string outfile,
	     const ParType & par);

// save the log-likelihood map to 3d file. May be used to detect candidate
// lesion.
int save_llmap(const ParType & par,
	       const vnl_vector<double> & llmap,
	       ImageType3DU::Pointer lindexPtr,
	       std::string filename);

// for debug purpose.
int save_priors(const ParType & par,
	       const vnl_matrix<double> & priors,
	       ImageType3DU::Pointer lindexPtr,
		std::string filename);

int save_image(ImageType3DF::Pointer ptr, std::string filename);

// print the final GMM parameters to file for registration part to use:
int print_gmm_to_file(const GMMType & gmm, bool isFG);

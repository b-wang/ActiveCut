// build a NxN array that contain the NLinks between (i,j) patches.

int build_adjmat(ImageType3DU::Pointer lindexPtr,
		       const vnl_matrix <double> & data, 
		       vnl_sparse_matrix <double> & con_map,
		       ParType & par);


int build_nlinks(const vnl_matrix <double> & data, 
		 const vnl_sparse_matrix <double> & con_map,
		 vnl_sparse_matrix <double> & nlinks_map,
		 ParType & par);

// do the graphcuts and save it in alpha.
unsigned graphcuts(vnl_vector<unsigned> & alpha,
		   const vnl_matrix <double> & data, 
		   const vnl_vector <unsigned> & hard_constraints,
		   const vnl_sparse_matrix <double> & nlinks_map,
		   const ParType & par);

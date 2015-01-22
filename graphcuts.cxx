#include <common.h>
#include <graph.h>
#include <gmm.h>
#include <utility.h>

int build_adjmat(ImageType3DU::Pointer lindexPtr,
		 const vnl_matrix <double> & data, 
		 vnl_sparse_matrix <double> & con_map,
		 ParType & par)
{
     ImageType3DU::IndexType lindexIdx, nbr_idx;

     // unsigned n_samples = data.rows();
     con_map.set_size(par.n_samples, par.n_samples);
     con_map.clear();

     // define neighborhood iterator
     typedef itk::ConstantBoundaryCondition<ImageType3DU>  BoundaryConditionType;
     typedef itk::NeighborhoodIterator< ImageType3DU, BoundaryConditionType > NeighborhoodIteratorType;
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType lindexIt(radius, lindexPtr, lindexPtr->GetLargestPossibleRegion());
     unsigned int nei_set_array[] = {4, 10, 12, 14, 16, 22, // 6 neighborhood
				     1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25, // 18 neighborhood
				     0, 2, 6, 8, 18, 20, 24, 26}; // 26 neighborhood

     if (par.n_nbrs != 6 && par.n_nbrs != 18 && par.n_nbrs != 26) {
	  printf("graphcuts(): number of neighbors must be 6, 18, or 26. Other values may give inacruate results!\n");
	  exit(1);
     }

     int cur_vox_id = 0, nbr_vox_id = 0;
     // if current voxel is boundary between patch i and j, add one on
     // con_map(i,j). 
     unsigned offset = 0;
     for (lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt) {
	  if ( lindexIt.GetCenterPixel() > 0) {
	       cur_vox_id = lindexIt.GetCenterPixel() - 1;
	       
	       for (unsigned neiIdx = 0; neiIdx < par.n_nbrs; neiIdx ++) {
		    offset = nei_set_array[neiIdx];
		    nbr_vox_id =  lindexIt.GetPixel(offset) - 1;

		    // make sure pair (i,j) only count once, since they are
		    // ordered. Only upper triangular part of the matrix is
		    // filled.
		    if ( (nbr_vox_id >= 0) && (nbr_vox_id > cur_vox_id) ) {
			 con_map.put(cur_vox_id, nbr_vox_id, con_map(cur_vox_id, nbr_vox_id) + 1);
		    }
	       } // neiIdx
	  } // in mask
     }
     return 0;
}

int build_nlinks(const vnl_matrix <double> & data, 
		 const vnl_sparse_matrix <double> & con_map,
		 vnl_sparse_matrix <double> & nlinks_map,
		 ParType & par)
{
     // compute beta. 
     double beta_sum = 0, n_edges = 0;
     unsigned row_id, col_id = 0;
     con_map.reset(); 
     while(con_map.next()) {
	  row_id = con_map.getrow();
	  col_id = con_map.getcolumn();
	  beta_sum += (data.get_row(row_id) - data.get_row(col_id)).squared_magnitude();
	  n_edges ++;
     }
     par.beta = 0.5 * (1 / (beta_sum/n_edges) );

     if (par.verbose >= 1) 
	  printf("build_nlinks(): beta = %f\n", par.beta);

     // convert the adjacency information int N-Links weights. 
     con_map.reset(); 
     double new_value = 0;
     nlinks_map.set_size(con_map.rows(), con_map.cols());
     while(con_map.next()) {
	  // get the patch sample id of a pair of neighbors.
	  row_id = con_map.getrow();
	  col_id = con_map.getcolumn();
	  
	  // the affinity value. we need it be roughtly in the range [0,
	  // 1]. First cube rooting the size, and get the edge length of each
	  // patch (roughtly, assuming a cube), then averaging two cube edge
	  // lengths. The rasio of the neighboring voxles (con_map) and the
	  // averaged edge length will be in the range [0, 1], roughly. 

	  // This may not well if two patches are long strips and close to each
	  // other. In such case, the cube root is small and the aff value will
	  // be big, often in range of [1, 20]. So, we should not set the gamma
	  // value too big. Gamma need to be smaller than 2, otherwise the
	  // neighboring links (nlinks_map just too strong such that the small
	  // lesion (foreground) will be absorbed into huge backgorund in order
	  // to acheieve smoothness). 
	  double aff = con_map.value();
	  new_value = par.gamma * aff * exp(- par.beta * par.beta0 * (data.get_row(row_id) - data.get_row(col_id)).squared_magnitude() );
	  nlinks_map.put(row_id, col_id, new_value); 
     };
     return 0;
}

unsigned graphcuts(vnl_vector<unsigned> & alpha,
		   const vnl_matrix <double> & data, 
		   const vnl_vector <unsigned> & hard_constraints,
		   const vnl_sparse_matrix <double> & nlinks_map,
		   const ParType & par)
{
     unsigned sample_idx = 0;
     Graph<double, double, double> G(1000, 1000);
     unsigned n_samples = data.rows();
     std::vector<unsigned> nodes(n_samples);
     for (sample_idx = 0; sample_idx < n_samples; sample_idx ++) {
	  nodes[sample_idx] = G.add_node();
     }
     
     // assign tlink weights. 
     double bg_link = 0, fg_link = 0;
     for (sample_idx = 0; sample_idx < n_samples; sample_idx ++) {
	  if (hard_constraints[sample_idx] == HC_FG) {
	       bg_link = 0;
	       fg_link = 6 * par.gamma + 1;
	  }
	  // either the initial background outside of bounding box (HC_BG), or the new
	  // constraint learnt from self-training/active-learning (HC_BG_NEW).
	  else if ((!par.baseline) && (hard_constraints[sample_idx] == HC_BG || hard_constraints[sample_idx] == HC_BG_NEW) ) {
	       bg_link = 6 * par.gamma + 1;
	       fg_link = 0;
	  }
	  else {
	       // if baseline = true, allow BG-->FG.
	       // must be unknown region.
	       fg_link = - gmm_eval_ll(data.get_row(sample_idx), par.gmm_bg, sample_idx);
	       bg_link = - gmm_eval_ll(data.get_row(sample_idx), par.gmm_fg, sample_idx);
	  }

	  G.add_tweights(nodes[sample_idx], fg_link, bg_link);
     }

     // assign nlink weights. 
     int cur_node = 0, nbr_node = 0, sample_id1 = 0, sample_id2 = 0;
     nlinks_map.reset(); 
     while(nlinks_map.next()) {
	  // get the patch sample id of a pair of neighbors.
	  sample_id1 = nlinks_map.getrow();
	  sample_id2 = nlinks_map.getcolumn();
	  cur_node = nodes[sample_id1];
	  nbr_node = nodes[sample_id2];
	  G.add_edge(cur_node, nbr_node, nlinks_map(sample_id1, sample_id2), nlinks_map(sample_id1, sample_id2));
     };

     // compute min-cut, max-flow
     G.maxflow();

     // update alpha.
     unsigned n_changes = 0;
     int sink_name = Graph<double, double, double>::SINK;
     int source_name = Graph<double, double, double>::SOURCE;
     for (sample_idx = 0; sample_idx < n_samples; sample_idx ++) {
	  unsigned old_value = alpha[sample_idx];

	  if (par.baseline) {
	       // All are unknow regions in baseline mode, so they can freely
	       // move into FG or BG. (actually there is no FG in
	       // hard_constraint map (trimap) for now, so only BG can go into
	       // FG).
	       if (G.what_segment(nodes[sample_idx]) == Graph<double, double, double>::SOURCE)
		    alpha[sample_idx] = ALPHA_FG;
	       else
		    alpha[sample_idx] = ALPHA_BG;
	       
	  }

	  else {
	       if (hard_constraints[sample_idx] == HC_FG) {
		    // make sure maxfow do not change those constraints. 
		    assert (G.what_segment(nodes[sample_idx]) == source_name);
		    alpha[sample_idx] = ALPHA_FG;
	       }
	       else if (hard_constraints[sample_idx] == HC_BG | hard_constraints[sample_idx] == HC_BG_NEW) {
		    assert(G.what_segment(nodes[sample_idx]) == sink_name);
		    alpha[sample_idx] = ALPHA_BG;
	       }
	       else {
		    // unknow regions.
		    if (G.what_segment(nodes[sample_idx]) == Graph<double, double, double>::SOURCE)
			 alpha[sample_idx] = ALPHA_FG;
		    else
			 alpha[sample_idx] = ALPHA_BG;
	       } // unknow regions.

	  } // else baseline. 

	  if (old_value != alpha[sample_idx])
	       n_changes ++;
     }

     return n_changes;
}


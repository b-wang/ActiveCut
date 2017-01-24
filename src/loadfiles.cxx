#include <common.h>
#include <gmm.h>
#include <graphcuts.h>
#include <utility.h>

int load_data(ImageType4DF::Pointer dataPtr,
	      ImageType3DI::Pointer maskPtr,
	      vnl_matrix<double> & D,
	      ImageType3DU::Pointer lindexPtr)
{
     unsigned spixel_id = 0;
     ImageType4DF::SizeType dataSize = dataPtr->GetLargestPossibleRegion().GetSize();
     IteratorType4DF dataIt(dataPtr, dataPtr->GetLargestPossibleRegion());
     ImageType4DF::IndexType dataIdx;

     ImageType3DI::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DI maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType3DI::IndexType maskIdx;

     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());

     unsigned n_channels = dataSize[3];
     unsigned n_samples = 0;

     // calculate number of voxels in brain mask. this code is stupid but it works. 
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       n_samples ++;
	  }
     }

     printf("load_data(): Total number of voxels in mask: %i.\n", n_samples);
     
     D.set_size(n_samples, n_channels);
     D.fill( 0 );
     
     unsigned n = 0; // voxel linear index. 
     for (maskIt.GoToBegin(), lindexIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ lindexIt) {
	  if (maskIt.Get() > 0) {
	       maskIdx = maskIt.GetIndex();
	       dataIdx[0] = maskIdx[0];
	       dataIdx[1] = maskIdx[1];
	       dataIdx[2] = maskIdx[2];
	       for (unsigned channel_id = 0; channel_id < n_channels; channel_id ++) {
		    dataIdx[3] = channel_id;
		    D[n][channel_id] += dataPtr->GetPixel(dataIdx);
	       }
	       n++;
	       lindexIt.Set(n);
	  } // in mask
     }
     return 0;
}

int load_constraints(ImageType3DU::Pointer lindexPtr,
		     ImageType3DC::Pointer initPtr,
		     vnl_vector<unsigned> & hard_constraints,
		     const ParType & par)
{
     unsigned vox_id = 0;
     ImageType3DU::SizeType lindexSize = lindexPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());
     ImageType3DU::IndexType lindexIdx;

     IteratorType3DC initIt(initPtr, initPtr->GetLargestPossibleRegion());
     ImageType3DC::IndexType initIdx;

     // the code is simplified by removing supervoxel support.
     for (lindexIt.GoToBegin(), initIt.GoToBegin(); !initIt.IsAtEnd(); ++ lindexIt, ++ initIt) {
	  if (lindexIt.Get() > 0) { // in mask
	       vox_id = lindexIt.Get() - 1;
	       if (initIt.Get() == HC_FG) { // fg
		    hard_constraints[vox_id] = HC_FG;
	       }
	       // either init'd background, or user given background.
	       else if (initIt.Get() == HC_BG) { 
		    hard_constraints[vox_id] = HC_BG;
	       }
	       else if  (initIt.Get() == HC_BG_NEW) {
		    hard_constraints[vox_id] = HC_BG_NEW;
	       }
	       else { // unknown region will be init'd as foreground!
		    hard_constraints[vox_id] = HC_UNKNOWN;
	       }
	  } // in mask
	  
     }
     return 0;
}

int load_priors(ImageType3DU::Pointer lindexPtr,
		std::string prior_file,
		vnl_matrix<double> & priors,
		const ParType & par)
{
     unsigned vox_id = 0;
     ImageType3DU::SizeType lindexSize = lindexPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());
     ImageType3DU::IndexType lindexIdx;

     unsigned n_comp = 0;

     // printf("fg: %i, bg: %i\n", par.gmm_fg.n_comp, par.gmm_bg.n_comp);
     if (par.gmm_fg.n_comp > par.gmm_bg.n_comp)
	  n_comp = par.gmm_fg.n_comp;
     else
	  n_comp = par.gmm_bg.n_comp;
     
     priors.set_size(par.n_samples, n_comp);
     priors.fill(0);

     // read 4D prior image. (last dim is num of Gaussian components)
     ReaderType4DF::Pointer priorReader = ReaderType4DF::New();
     priorReader->SetFileName(prior_file);
     priorReader->Update();
     ImageType4DF::Pointer priorPtr = priorReader->GetOutput();
     ImageType4DF::SizeType priorSize = priorPtr->GetLargestPossibleRegion().GetSize();
     IteratorType4DF priorIt(priorPtr, priorPtr->GetLargestPossibleRegion());
     ImageType4DF::IndexType priorIdx;

     for (lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt) {
	  if (lindexIt.Get() > 0) { // in mask
	       vox_id = lindexIt.Get() - 1;
	       lindexIdx = lindexIt.GetIndex();
	       priorIdx[0] = lindexIdx[0];
	       priorIdx[1] = lindexIdx[1];
	       priorIdx[2] = lindexIdx[2];
	       for (priorIdx[3] = 0; priorIdx[3] < priorSize[3]; priorIdx[3] ++) {
		    priors(vox_id, priorIdx[3]) = priorPtr->GetPixel(priorIdx);
	       }
	  } // in mask.
     }

     // normalize the atlas
     vnl_vector <double> zero_vec(priors.cols(), 0);
     for (unsigned n = 0; n < par.n_samples; n ++) {
	  if (priors.get_row(n).sum() == 0) {
	       priors.set_row(n, zero_vec);
	  }
	  else {
	       priors.scale_row(n, 1/priors.get_row(n).sum() );
	  }
     }
     return 0;
}

int align_priors(vnl_matrix<double> & gmm_labels,
		 const vnl_matrix<double> & priors,
		 const vnl_vector<unsigned> & alpha,
		 unsigned whatground,
		 const ParType & par)
{
     unsigned n_comp = 0;
     if (whatground == ALPHA_FG) n_comp = par.gmm_fg.n_comp;
     else n_comp = par.gmm_bg.n_comp;
     std::vector<unsigned> permu_map(n_comp, 0);

     for (unsigned k = 0; k < n_comp; k ++) {
	  permu_map[k] = k;
     }
     
     sort(permu_map.begin(), permu_map.end() );
     std::vector<unsigned> permu_map_best(permu_map);
     unsigned k_new = 0, k_old = 0;
     unsigned n_matches = 0;
     unsigned n_matches_best = 0;
     
     do {
	  n_matches = 0;
	  for (unsigned sample_id = 0; sample_id < gmm_labels.rows(); sample_id ++) {
     	       if (alpha[sample_id] == whatground) {
		    k_old = gmm_labels.get_row(sample_id).arg_max();
		    if (permu_map[k_old] == priors.get_row(sample_id).arg_max()) {
			 n_matches ++;
		    }
	       }
	  } // sample id
	  if (n_matches > n_matches_best) {
	       n_matches_best  = n_matches;
	       permu_map_best = permu_map;
	       if (par.verbose >= 2) {
		    for (unsigned k = 0; k < n_comp; k ++) {
			 printf("align_priors(): %i->%i, ", k, permu_map[k]);
		    }
		    printf("n_matches_best: %i\n", n_matches_best);
	       } // verbose
	  }
     } while (next_permutation(permu_map.begin(), permu_map.end()) );

     if (par.verbose >= 0) {
	  printf("align_priors(): Best permutation: ");
	  for (unsigned k = 0; k < n_comp; k ++) {
	       printf("%i->%i, ", k, permu_map_best[k]);
	  }
	  printf("\n");
     } // verbose

     // update gmm_labels. 

     for (unsigned sample_id = 0; sample_id < gmm_labels.rows(); sample_id ++) {     
	  if (alpha[sample_id] == whatground) {
	       k_old = gmm_labels.get_row(sample_id).arg_max();
	       k_new = permu_map_best[k_old];
	       gmm_labels(sample_id, k_old) = 0;
	       gmm_labels(sample_id, k_new) = 1;
	  }
     }
     return 0;
}


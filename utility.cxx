// for saving output images. e.g. convert alpha (seg for now) vector to 3D images. and other utility fun.

#include <common.h>
#include <utility.h>

template <class T>
int save_asimage(itk::Image<T, 3> ptr, 
		 vnl_vector<T> array,
		 std::string fname)
{
     
}


int save_image(ImageType3DF::Pointer ptr, std::string filename)
{

     WriterType3DF::Pointer writer = WriterType3DF::New();
	  
     writer->SetInput(ptr);
     writer->SetFileName(filename);
     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
     } 

     std::cout << "save_gmm_labelmap(): File " << filename << " saved.\n";

     return 0;
}


double compute_beta(const vnl_matrix <double> & data, 
		    vnl_sparse_matrix <double> & nlinks_map)
{
     unsigned n_samples = nlinks_map.rows();
     double beta_sum = 0, n_edges = 0;
     for (unsigned sample_idx = 0; sample_idx < n_samples; sample_idx ++) {
	  for (unsigned nbr_idx = sample_idx; nbr_idx < n_samples; nbr_idx ++) {
	       if (nlinks_map(sample_idx, nbr_idx) > 0) {
		    beta_sum += (data.get_row(sample_idx) - data.get_row(nbr_idx)).squared_magnitude();
		    n_edges ++;
	       }
	  }
     }
     
     return 0.5 * (1 / (beta_sum/n_edges) );
}

int print_par(ParType par)
{
     print_gmm(par.gmm_fg);
     printf("\n");
     print_gmm(par.gmm_bg);
     printf("\n");     
     if (par.verbose >= 1) {
	  printf("n_channels: %i, n_samples: %i, gamma: %.2f, beta: %E, seed: %i, kmeansruns: %i, verbose: %i\n", par.n_channels, par.n_samples, par.gamma, par.beta, par.seed, par.kmeansruns, par.verbose);
     }
	 return 1;
}

int print_gmm(const GMMType & gmm)
{
     unsigned n_channels = gmm.comp[0].mu.size();
     printf("%s n_comp: %i, n_pts: %i\n", gmm.name.c_str(), gmm.n_comp, gmm.n_pts);
     // print pi (pi is now different per voxel) so do not print it.
     // printf("pi = \n");
     // for (unsigned k = 0; k < gmm.n_comp; k ++) 
     // 	  printf("%10.2f ", gmm.pi[k]); 
     // printf("\n");

     printf(" N_k (numPts) = \n");
     for (unsigned k = 0; k < gmm.n_comp; k ++) 
     	  printf("%10.2f ", gmm.comp[k].numPts); 
     printf("\n");

     for (unsigned k = 0; k < gmm.n_comp; k ++) {
	  // print mu
	  printf("mu [%i] = \n", k);
	  for (unsigned i = 0; i < n_channels; i ++) 
	       printf("%10.2f ", gmm.comp[k].mu[i]);
	  printf("\n");


	  // print cov and inv_cov
	  printf("cov[%i] and inv_cov[%i] = \n", k, k);
	  for (unsigned i = 0; i < n_channels; i ++) {
	       for (unsigned j = 0; j < n_channels; j ++) {
		    printf("%10.2f ", gmm.comp[k].cov(i,j));
	       }
	       printf("\t\t");
	       for (unsigned j = 0; j < n_channels; j ++) {
		    printf("%10.2E ", gmm.comp[k].inv_cov(i,j));
	       }
	       printf("\n");
	  }
     }
	 return 1;
}

int print_vnl_matrix(const vnl_matrix<double> & mat, unsigned r, unsigned c)
{

     for (unsigned ri = 0; ri < r; ri ++) {
	  for (unsigned ci = 0; ci < c; ci++) {
	       printf("%f ", mat(ri, ci));
	  }
	  printf("\n");
     }
     // mat.print(vcl_cout);
	 return 1;
}

int print_vnl_vec(const vnl_vector<double> & vec, unsigned N)
{
     for (int i = 0; i < N; i++) {
	  printf("%f ", vec[i]);
     }
     printf("\n");
	 return 1;
}

int save_gmm_labelmap(vnl_matrix<double> & gmm_labels, // initial gmm labels.
		      ParType & par,
		      const vnl_vector<unsigned> & alpha,
		      unsigned whatground,
		      ImageType3DU::Pointer lindexPtr,
		      std::string filename)
{
     // define some variables used for initializing output images. 
     ImageType3DC::IndexType start;
     start[0] = 0;
     start[1] = 0;
     start[2] = 0;

     ImageType3DU::SizeType lindexSize = lindexPtr->GetLargestPossibleRegion().GetSize();
     ImageType3DC::RegionType volRegion;
     volRegion.SetSize(lindexSize);
     volRegion.SetIndex(start);

     // create gmm label image. 
     ImageType3DS::Pointer labelPtr = ImageType3DS::New();
     labelPtr->SetRegions(volRegion);
     labelPtr->Allocate();
     labelPtr->FillBuffer( 0 ); // init to zero.

     labelPtr->SetOrigin( lindexPtr->GetOrigin() );
     labelPtr->SetSpacing(lindexPtr->GetSpacing() );
     labelPtr->SetDirection(lindexPtr->GetDirection() );

     IteratorType3DS labelIt(labelPtr, labelPtr->GetLargestPossibleRegion());
     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());

     unsigned label = 0, sample_id = 0;
     for (labelIt.GoToBegin(), lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt, ++labelIt) {
	  if (lindexIt.Get() > 0) {
	       sample_id = lindexIt.Get() - 1;
	       if (alpha[sample_id] == whatground) {
		    label = gmm_labels.get_row(sample_id).arg_max();
		    labelIt.Set(label + 1); // convert to 1-based.
	       }
	  }
     }

     WriterType3DS::Pointer writer = WriterType3DS::New();
	  
     writer->SetInput(labelPtr);
     writer->SetFileName(filename);
     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
     } 

     std::cout << "save_gmm_labelmap(): File " << filename << " saved.\n";

     return 0;
}

int save_gmm_posterior(vnl_matrix<double> & gmm_labels, // initial gmm labels.
		      ParType & par,
		      const vnl_vector<unsigned> & alpha,
		      unsigned whatground,
		      ImageType3DU::Pointer lindexPtr,
		      std::string filename)
{
     // define some variables used for initializing output images. 
     ImageType4DF::IndexType start;
     start.Fill(0);
     ImageType3DU::SizeType lindexSize = lindexPtr->GetLargestPossibleRegion().GetSize();
     ImageType3DU::IndexType lindexIdx;

     ImageType4DF::SizeType labelSize;
     labelSize[0] = lindexSize[0];
     labelSize[1] = lindexSize[1];
     labelSize[2] = lindexSize[2];
     if (whatground == ALPHA_FG) labelSize[3] = par.gmm_fg.n_comp;
     else labelSize[3] = par.gmm_bg.n_comp;

     ImageType4DF::RegionType volRegion;
     volRegion.SetSize(labelSize);
     volRegion.SetIndex(start);

     // create gmm label image. 
     ImageType4DF::Pointer labelPtr = ImageType4DF::New();
     labelPtr->SetRegions(volRegion);
     labelPtr->Allocate();
     labelPtr->FillBuffer( 0 ); // init to zero.
     ImageType4DF::IndexType labelIdx;

     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());

     unsigned sample_id = 0;
     for (lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt) {
	  if (lindexIt.Get() > 0) {
	       sample_id = lindexIt.Get() - 1;
	       if (alpha[sample_id] == whatground) {
		    lindexIdx = lindexIt.GetIndex();
		    labelIdx[0] = lindexIdx[0];
		    labelIdx[1] = lindexIdx[1];
		    labelIdx[2] = lindexIdx[2];
		    for (unsigned k = 0; k < labelSize[3]; k ++) {
			 labelIdx[3] = k;
			 labelPtr->SetPixel(labelIdx, gmm_labels(sample_id, k) );
		    }
	       }
	  }
     }

     WriterType4DF::Pointer writer = WriterType4DF::New();
	  
     writer->SetInput(labelPtr);
     writer->SetFileName(filename);
     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
     } 

     std::cout << "save_gmm_posterior(): File " << filename << " saved.\n";

     return 0;
}


int save_alpha(const ParType & par,
	       const vnl_vector<unsigned> & alpha,
	       ImageType3DU::Pointer lindexPtr,
	       std::string filename)
{
     ImageType3DU::RegionType lindexRegion = lindexPtr->GetLargestPossibleRegion();

     // create label label image. 
     ImageType3DS::Pointer labelPtr = ImageType3DS::New();
     labelPtr->SetRegions(lindexRegion);
     labelPtr->Allocate();
     labelPtr->FillBuffer( 0 ); // init to zero.

     labelPtr->SetOrigin( lindexPtr->GetOrigin() );
     labelPtr->SetSpacing(lindexPtr->GetSpacing() );
     labelPtr->SetDirection(lindexPtr->GetDirection() );

     IteratorType3DS labelIt(labelPtr, labelPtr->GetLargestPossibleRegion());
     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());

     unsigned label = 0, sample_id = 0;
     for (labelIt.GoToBegin(), lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt, ++labelIt) {
	  if (lindexIt.Get() > 0) {
	       sample_id = lindexIt.Get() - 1;
	       if (alpha[sample_id] == ALPHA_FG)
		    labelIt.Set(1);
	       else 
		    labelIt.Set(0);
	  }
     }

     WriterType3DS::Pointer writer = WriterType3DS::New();
	  
     writer->SetInput(labelPtr);
     writer->SetFileName(filename);
     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
     } 

     std::cout << "save_alpha(): File " << filename << " saved.\n";

     return 0;
}

int save_patches(ImageType3DU::Pointer lindexPtr,
	     const vnl_matrix <double> & data,
	     std::string outfile,
	     const ParType & par)
{
     ImageType3DU::SizeType lindexSize = lindexPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());
     ImageType3DU::IndexType lindexIdx;

     unsigned spixel_id = 0;
     ImageType4DF::IndexType dataIdx;
     ImageType4DF::SizeType dataSize;
     dataSize[0] = lindexSize[0];
     dataSize[1] = lindexSize[1];
     dataSize[2] = lindexSize[2];
     dataSize[3] = data.cols();

     ImageType4DF::IndexType start;
     start.Fill(0);
     
     ImageType4DF::RegionType volRegion;
     volRegion.SetSize(dataSize);
     volRegion.SetIndex(start);

     // create gmm label image. 
     ImageType4DF::Pointer dataPtr = ImageType4DF::New();
     dataPtr->SetRegions(volRegion);
     dataPtr->Allocate();
     dataPtr->FillBuffer( 0 ); // init to zero.

     for (lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt) {
	  if (lindexIt.Get() > 0) {
	       lindexIdx = lindexIt.GetIndex();
	       spixel_id = lindexIt.Get() - 1; // convert to zero-based indexing.
	       dataIdx[0] = lindexIdx[0];
	       dataIdx[1] = lindexIdx[1];
	       dataIdx[2] = lindexIdx[2];
	       for (dataIdx[3] = 0; dataIdx[3] < dataSize[3]; dataIdx[3] ++) {
		    dataPtr->SetPixel(dataIdx, data(spixel_id, dataIdx[3]));
	       }
	  } // in mask
     }

     WriterType4DF::Pointer writer = WriterType4DF::New();
	  
     writer->SetInput(dataPtr);
     writer->SetFileName(outfile);
     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
     } 

     std::cout << "save_patches(): File " << outfile << " saved.\n";

     return 0;
	  
}

int save_llmap(const ParType & par,
	       const vnl_vector<double> & llmap,
	       ImageType3DU::Pointer lindexPtr,
	       std::string filename)
{
     ImageType3DU::RegionType lindexRegion = lindexPtr->GetLargestPossibleRegion();

     // create label label image. 
     ImageType3DF::Pointer llmapPtr = ImageType3DF::New();
     llmapPtr->SetRegions(lindexRegion);
     llmapPtr->Allocate();
     llmapPtr->FillBuffer( 0 ); // init to zero.

     llmapPtr->SetOrigin( lindexPtr->GetOrigin() );
     llmapPtr->SetSpacing(lindexPtr->GetSpacing() );
     llmapPtr->SetDirection(lindexPtr->GetDirection() );

     IteratorType3DF llmapIt(llmapPtr, llmapPtr->GetLargestPossibleRegion());
     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());

     unsigned sample_id = 0;
     // normalize to [0, 1].
     double min = llmap.min_value(), max = llmap.max_value();
     for (llmapIt.GoToBegin(), lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt, ++llmapIt) {
	  if (lindexIt.Get() > 0) {
	       sample_id = lindexIt.Get() - 1;
	       llmapIt.Set( (llmap[sample_id] - min) / (max - min) );
	  }
     }

     WriterType3DF::Pointer writer = WriterType3DF::New();
	  
     writer->SetInput(llmapPtr);
     writer->SetFileName(filename);
     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
     } 

     std::cout << "save_alpha(): File " << filename << " saved.\n";

     return 0;
}

int save_priors(const ParType & par,
	       const vnl_matrix<double> & priors,
	       ImageType3DU::Pointer lindexPtr,
	       std::string filename)
{
     ImageType3DU::RegionType lindexRegion = lindexPtr->GetLargestPossibleRegion();

     // create prior image. 
     ImageType3DU::Pointer priorPtr = ImageType3DU::New();
     priorPtr->SetRegions(lindexRegion);
     priorPtr->Allocate();
     priorPtr->FillBuffer( 0 ); // init to zero.

     priorPtr->SetOrigin( lindexPtr->GetOrigin() );
     priorPtr->SetSpacing(lindexPtr->GetSpacing() );
     priorPtr->SetDirection(lindexPtr->GetDirection() );

     IteratorType3DU priorIt(priorPtr, priorPtr->GetLargestPossibleRegion());
     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());

     unsigned sample_id = 0;
     for (priorIt.GoToBegin(), lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt, ++priorIt) {
	  if (lindexIt.Get() > 0) {
	       sample_id = lindexIt.Get() - 1;
	       priorIt.Set(priors.get_row(sample_id).arg_max() + 1);
	  }
     }

     WriterType3DU::Pointer writer = WriterType3DU::New();
	  
     writer->SetInput(priorPtr);
     writer->SetFileName(filename);
     try 
     { 
	  writer->Update(); 
     } 
     catch( itk::ExceptionObject & err ) 
     { 
	  std::cerr << "ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
     } 

     std::cout << "save_alpha(): File " << filename << " saved.\n";

     return 0;
}


int print_gmm_to_file(const GMMType & gmm, bool isFG)
{
     FILE * pFile;
     char fileName[] = "GMM_parameters_";
     if(isFG)
     {
	  strcat(fileName,"FG.txt");
     }
     else
     {
	  strcat(fileName,"BG.txt");
     }

     std::cout << "saved GMM parameters file name: " << fileName << "\n";

     pFile = fopen (fileName ,"w");

     unsigned n_channels = gmm.comp[0].mu.size();
     fprintf(pFile, "%i\n", gmm.n_comp);
     fprintf(pFile, "\n");

     for (unsigned k = 0; k < gmm.n_comp; k ++) {
	  // print mu
	  for (unsigned i = 0; i < n_channels; i ++) 
	       fprintf(pFile, "%10.2f ", gmm.comp[k].mu[i]);
	  fprintf(pFile, "\n");
	  fprintf(pFile, "\n");

	  // print cov
	  for (unsigned i = 0; i < n_channels; i ++) {
	       for (unsigned j = 0; j < n_channels; j ++) {
		    fprintf(pFile, "%10.2f ", gmm.comp[k].cov(i,j));
	       }
	       fprintf(pFile, "\n");
	  }
	  fprintf(pFile, "\n");
     }

     fclose(pFile);
	 return 1;
}

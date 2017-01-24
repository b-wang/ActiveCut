#include <common.h>
#include <gmm.h>
#include <query.h>
#include <loadfiles.h>
// include libraries and define some struts for CCP
#include <vector>       // std::vector

// the comparator for sorting
bool myCompFunction (int i,int j) { return (i>j); }

// the pair data structure for returning indices in sorting
bool myComparatorforCcp ( const mypair& l, const mypair& r)
   { return l.first > r.first; }
bool myComparatorforCcpFI ( const mypairFI& l, const mypairFI& r)
   { return l.first > r.first; }
 

// for mathematical morphology operations - erosion
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkConnectedComponentImageFilter.h"
// for duplicate image
#include "itkImageDuplicator.h"

using namespace std;

int logistic(vnl_vector<double> & score_map,
	     const vnl_matrix<double> & data,
	     vnl_sparse_matrix <double> & con_map, // can not set const due to a vnl bug.
	     const vnl_vector<unsigned> & alpha,
	     const ParType & par)
{
     unsigned nbr_id = 0;
     vnl_sparse_matrix<double>::row con_map_row;
     vnl_vector<double> exp_term (2, 0);
     vnl_vector<double> score_map_old(score_map);
     double changed_score = 1e7;
     unsigned iter = 0;
     while(changed_score > 0.001 && iter < par.logimax) {
	  iter ++;
	  score_map_old = score_map; // save previous scores.
	  for (unsigned n = 0; n < par.n_samples; n ++) {
	       // compute difference of log-likelihood.
	       exp_term[0] = gmm_eval_ll(data.get_row(n), par.gmm_fg, n);
	       exp_term[1] = gmm_eval_ll(data.get_row(n), par.gmm_bg, n);

	       // compute diff of prior.  we don't bother looping the 2
	       // iterations FG and BG. just repleat it twice.
	       con_map_row = con_map.get_row(n);		    
	       for (vnl_sparse_matrix<double>::row::const_iterator col_iter = con_map_row.begin(); col_iter != con_map_row.end();  ++ col_iter) {
		    nbr_id = (*col_iter).first;
		    // if score = 1 (FG)
		    exp_term[0] += par.eta * (1 * score_map[nbr_id] + 0 *(1-score_map[nbr_id]));
		    // the above line does not take into account the
		    // patches. Int assume each patch is just a single voxel. To
		    // account for that, need to follow the build_nliks func.
	       }

	       for (vnl_sparse_matrix<double>::row::const_iterator col_iter = con_map_row.begin(); col_iter != con_map_row.end();  ++ col_iter) {
		    nbr_id = (*col_iter).first;
		    // if score == 0 (BG)
		    exp_term[1] += par.eta * (0 * score_map[nbr_id] + 1 * (1-score_map[nbr_id]));
	       }

	       score_map[n] = exp_term[0] - exp_term[1];
	       score_map[n] = 1 / (1 + exp(- score_map[n]));
	  } // n i.e. sample id.
	  changed_score = (score_map_old - score_map).two_norm() / score_map.two_norm();
	  if (par.verbose >= 0) {
	       printf("logistic(): iteration %i, changed_score: %f\n", iter, changed_score);
	  }
     }
     return 0;
}


int logistic_init(vnl_vector<double> & score_map,
		  const vnl_matrix<double> & data,
		  vnl_sparse_matrix <double> & con_map, // can not set const due to a vnl bug.
		  const vnl_vector<unsigned> & alpha,
		  const ParType & par,
		  double smalleta,
		  unsigned n_iter)
{
     unsigned nbr_id = 0;
     vnl_sparse_matrix<double>::row con_map_row;
     vnl_vector<double> score_map_old(score_map);
     vnl_vector<double> exp_term (2, 0);
     unsigned iter = 0;
     double changed_score = 1e7;

     for (iter = 0; iter < n_iter; iter ++) {
	  score_map_old = score_map; // save previous scores.
	  for (unsigned n = 0; n < par.n_samples; n ++) {
	       // compute difference of log-likelihood.
	       exp_term[0] = gmm_eval_ll(data.get_row(n), par.gmm_fg, n);
	       exp_term[1] = gmm_eval_ll(data.get_row(n), par.gmm_bg, n);

	       // compute diff of prior.  we don't bother looping the 2
	       // iterations FG and BG. just repleat it twice.
	       con_map_row = con_map.get_row(n);		    
	       for (vnl_sparse_matrix<double>::row::const_iterator col_iter = con_map_row.begin(); col_iter != con_map_row.end();  ++ col_iter) {
		    nbr_id = (*col_iter).first;
		    // if score = 1 (FG)
		    exp_term[0] += smalleta * (1 * score_map[nbr_id] + 0 *(1-score_map[nbr_id]));
		    // the above line does not take into account the
		    // patches. Int assume each patch is just a single voxel. To
		    // account for that, need to follow the build_nliks func.
	       }

	       for (vnl_sparse_matrix<double>::row::const_iterator col_iter = con_map_row.begin(); col_iter != con_map_row.end();  ++ col_iter) {
		    nbr_id = (*col_iter).first;
		    // if score == 0 (BG)
		    exp_term[1] += smalleta * (0 * score_map[nbr_id] + 1 * (1-score_map[nbr_id]));
	       }

	       score_map[n] = exp_term[0] - exp_term[1];
	       score_map[n] = 1 / (1 + exp(- score_map[n]));
	  } // n i.e. sample id.

	  changed_score = (score_map_old - score_map).two_norm() / score_map.two_norm();

	  if (par.verbose >= 1) {
	       printf("logistic_init(): small eta = %f, teration %i, changed score = %f\n", smalleta, iter, changed_score);
	  }
     }
     return 0;
}

void preprocessingBinaryPredicProb(ImageType3DI::Pointer binaryPredictProb)
{
     std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
     std::cout << "Do preprocessing for binary image." << std::endl;

     const ImageType3DI::SizeType sizeOfImage = binaryPredictProb->GetLargestPossibleRegion().GetSize();
     
     int width = sizeOfImage[0];
     int height = sizeOfImage[1];
     int slice = sizeOfImage[2];
     
     ImageType3DI::RegionType outRegion = binaryPredictProb->GetLargestPossibleRegion();

     int i,j,k;

     ImageType3DUC::Pointer binaryPredictProbUC = ImageType3DUC::New();
     binaryPredictProbUC->SetRegions(outRegion);
     binaryPredictProbUC->Allocate();
     binaryPredictProbUC->FillBuffer( 0 );

     ImageType3DUC::IndexType imageIndexUC;
     
     ImageType3DI::PixelType imageValueDI;
     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    imageValueDI = 0;
		    imageIndexUC[0] = i;
		    imageIndexUC[1] = j;
		    imageIndexUC[2] = k;

		    imageValueDI = binaryPredictProb->GetPixel( imageIndexUC );
		    if(imageValueDI > 0)
		    {
			binaryPredictProbUC->SetPixel(imageIndexUC, 1 );
		    }
	       }
	  }
     }

     // do large erosion

     const unsigned int Dimension = 3;

     typedef itk::BinaryBallStructuringElement<PixelType3DUC, Dimension> StructuringElementType;
     typedef itk::BinaryErodeImageFilter<ImageType3DUC, ImageType3DUC, StructuringElementType> ErodeFilterType;
	  
     ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();

     StructuringElementType structuringElement;

     structuringElement.SetRadius( 3 ); // 3x3 structuring element
     structuringElement.CreateStructuringElement();
     binaryErode->SetKernel( structuringElement );
     binaryErode->SetInput( binaryPredictProbUC );
     binaryErode->SetErodeValue( 1 );

     ImageType3DUC::Pointer tempBinaryPredictErosion;
     tempBinaryPredictErosion = binaryErode->GetOutput();

     // std::string test_save_5 = "test_original_binary_Predict.nii.gz";
     // WriterType3DUC::Pointer writer5 = WriterType3DUC::New();
     // writer5->SetFileName( test_save_5 );
     // writer5->SetInput( binaryPredictProbUC ) ;
     // writer5->Update();

     ///////////////////////////////////////////
     // do dilation (same radius of structure as erosion)
     typedef itk::BinaryDilateImageFilter<ImageType3DUC, ImageType3DUC, StructuringElementType> DilateFilterType;

     DilateFilterType::Pointer binaryDilate = DilateFilterType::New();

     StructuringElementType structuringElementDilation;

     structuringElementDilation.SetRadius( 3 ); // 3x3 structuring element
     structuringElementDilation.CreateStructuringElement();

     binaryDilate->SetKernel( structuringElementDilation );
     binaryDilate->SetInput( tempBinaryPredictErosion );
     binaryDilate->SetDilateValue( 1 );
     ImageType3DUC::Pointer tempImagePointer = binaryDilate->GetOutput();

     // use extract image filter instead of saving the dilated volume - it works! but slow
     /*
     typedef itk::ExtractImageFilter< ImageType3DUC, ImageType3DUC > FilterType;
     FilterType::Pointer filter = FilterType::New();
     filter->SetExtractionRegion(outRegion);
     filter->SetInput(tempImagePointer);
     filter->SetDirectionCollapseToIdentity(); // This is required.
     filter->Update();
     ImageType3DUC::Pointer tempBinaryPredictErosionThenDilation = filter->GetOutput();
     */

     // old way: saving the dilated volume -- it works by slow.
     
     ImageType3DUC::Pointer tempBinaryPredictErosionThenDilation;
     tempBinaryPredictErosionThenDilation = binaryDilate->GetOutput();
     std::string test_save_5 = "test_erosion_and_dilation_binary_Predict.nii.gz";
     WriterType3DUC::Pointer writer5 = WriterType3DUC::New();
     writer5->SetFileName( test_save_5 );
     writer5->SetInput( tempBinaryPredictErosionThenDilation ) ;
     writer5->Update();
     

     // modify the content of binaryPredictProb
     
     ImageType3DUC::PixelType imageValueUC;

     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    imageValueDI = 0;
		    imageIndexUC[0] = i;
		    imageIndexUC[1] = j;
		    imageIndexUC[2] = k;

		    imageValueUC = tempBinaryPredictErosionThenDilation->GetPixel( imageIndexUC );

		    binaryPredictProb->SetPixel(imageIndexUC, imageValueUC );

	       }
	  }
     }

}

void morphologicalProcessingTrueCandidatesVolume(ImageType3DUC::Pointer trueCandidatesVolume,
						 ImageType3DUC::Pointer trueCandidatesVolumeErosion, 
						 ImageType3DUC::Pointer trueCandidatesVolumeDilationMinusErosion)
{
     std::cout << "Do morphological processing for true candidates volume to get masks." << std::endl;

     const ImageType3DI::SizeType sizeOfImage = trueCandidatesVolume->GetLargestPossibleRegion().GetSize();
     
     int width = sizeOfImage[0];
     int height = sizeOfImage[1];
     int slice = sizeOfImage[2];

     ImageType3DI::RegionType outRegion = trueCandidatesVolume->GetLargestPossibleRegion();

     int i,j,k;

     ImageType3DUC::IndexType imageIndexUC;
     ImageType3DUC::PixelType imageValueUCerosion, imageValueUCdelision;

     const unsigned int Dimension = 3;

     // do erosion = trueCandidatesVolumeErosion
     /*
     typedef itk::BinaryBallStructuringElement<PixelType3DUC, Dimension> StructuringElementType;
     typedef itk::BinaryErodeImageFilter<ImageType3DUC, ImageType3DUC, StructuringElementType> ErodeFilterType;
	  
     ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();

     StructuringElementType structuringElement;

     structuringElement.SetRadius( 1 ); // 3x3 structuring element
     structuringElement.CreateStructuringElement();
     binaryErode->SetKernel( structuringElement );
     binaryErode->SetInput( trueCandidatesVolume );
     binaryErode->SetErodeValue( 1 );

     ImageType3DUC::Pointer tempTrueCandidatesVolumeErosion;
     tempTrueCandidatesVolumeErosion = binaryErode->GetOutput();

     std::string test_save_4 = "test_true_candidates_erosion.nii.gz";
     WriterType3DUC::Pointer writer4 = WriterType3DUC::New();
     writer4->SetFileName( test_save_4 );
     writer4->SetInput( tempTrueCandidatesVolumeErosion ) ;
     writer4->Update();
     */

     // no erosion
     ImageType3DUC::Pointer tempTrueCandidatesVolumeErosion;
     tempTrueCandidatesVolumeErosion = trueCandidatesVolume;

     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    imageValueUCerosion = 0;
		    imageIndexUC[0] = i;
		    imageIndexUC[1] = j;
		    imageIndexUC[2] = k;
		    imageValueUCerosion = tempTrueCandidatesVolumeErosion->GetPixel( imageIndexUC );

		    trueCandidatesVolumeErosion->SetPixel(imageIndexUC, imageValueUCerosion );
	       }
	  }
     }
     

     // do dilation

     typedef itk::BinaryBallStructuringElement<PixelType3DUC, Dimension> StructuringElementTypeDilation;
     typedef itk::BinaryDilateImageFilter<ImageType3DUC, ImageType3DUC, StructuringElementTypeDilation> DilateFilterType;
     DilateFilterType::Pointer binaryDilate = DilateFilterType::New();

     StructuringElementTypeDilation structuringElementDilation;
     structuringElementDilation.SetRadius( 1 ); // 3x3 structuring element
     structuringElementDilation.CreateStructuringElement();

     binaryDilate->SetKernel( structuringElementDilation );
     binaryDilate->SetInput( trueCandidatesVolume );
     binaryDilate->SetDilateValue( 1 );

     // use extract image filter instead of saving the dilated volume - don't use this!
     /*
     ImageType3DUC::Pointer tempImagePointer = binaryDilate->GetOutput();
     typedef itk::ExtractImageFilter< ImageType3DUC, ImageType3DUC > FilterType;
     FilterType::Pointer filter = FilterType::New();
     filter->SetExtractionRegion(outRegion);
     filter->SetInput(tempImagePointer);
     filter->SetDirectionCollapseToIdentity(); // This is required.
     filter->Update();
     ImageType3DUC::Pointer trueCandidatesVolumeDilation = filter->GetOutput();
     */

     // use this to make sure no segmentation fault and no wrong result.
     ImageType3DUC::Pointer trueCandidatesVolumeDilation = ImageType3DUC::New();
     trueCandidatesVolumeDilation->SetRegions(outRegion);
     trueCandidatesVolumeDilation->Allocate();
     trueCandidatesVolumeDilation->FillBuffer( 0 );
     trueCandidatesVolumeDilation = binaryDilate->GetOutput();
     std::string test_save_5 = "test_true_candidates_dilation.nii.gz";
     WriterType3DUC::Pointer writer5 = WriterType3DUC::New();
     writer5->SetFileName( test_save_5 );
     writer5->SetInput( trueCandidatesVolumeDilation ) ;
     writer5->Update();
     

     // dilation - erosion = trueCandidatesVolumeDilationMinusErosion

     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    imageValueUCerosion = 0;
		    imageValueUCdelision = 0;
		    imageIndexUC[0] = i;
		    imageIndexUC[1] = j;
		    imageIndexUC[2] = k;
		    imageValueUCerosion = trueCandidatesVolumeErosion->GetPixel( imageIndexUC );
		    imageValueUCdelision = trueCandidatesVolumeDilation->GetPixel( imageIndexUC );

		    if(imageValueUCerosion < 1 && imageValueUCdelision > 0)
		    {
			 trueCandidatesVolumeDilationMinusErosion->SetPixel(imageIndexUC, 1 );
		    }
	       }
	  }
     }

     // std::string test_save_6 = "test_true_candidates_dilationMinusErosion.nii.gz";
     // WriterType3DUC::Pointer writer6 = WriterType3DUC::New();
     // writer6->SetFileName( test_save_6 );
     // writer6->SetInput( trueCandidatesVolumeDilationMinusErosion ) ;
     // writer6->Update();

}


bool computeShapeScoreforCcpCandidates (ImageType3DI::Pointer ccpPredictProb,
					ImageType3DC::Pointer trimap,
					const int& numTopCcptoEvaluate,
					std::vector< std::pair<int, int> >& ccpVectorPairs,
					ConnectedCompInfor allConnectedComponents[],
					const float& thresholdScoreforTopCandidates)


{
     std::cout << "Compute shape information for Top ranked connected componets:" << std::endl;

     const ImageType3DI::SizeType sizeOfImage = ccpPredictProb->GetLargestPossibleRegion().GetSize();
     
     int width = sizeOfImage[0];
     int height = sizeOfImage[1];
     int slice = sizeOfImage[2];

     ImageType3DI::RegionType outRegion = ccpPredictProb->GetLargestPossibleRegion();

     // print the sorted connected components based on volume
     int i,j,k,n;
     int curLabelValue = 0;
     ImageType3DI::IndexType ccpIndex;
     ImageType3DI::PixelType ccpLabelValue;

     const unsigned int Dimension = 3;

     float curScoreCriteria = 0;
     float curShapeScore = 0;
     float shapeScoreScaling = 1.0;
     float curVolumeMultiplyProb = 0;
     
     //mypairFI * const ccpTopPair = (mypairFI*)_alloca(numTopCcptoEvaluate * sizeof(mypairFI)); // got error in linux
     mypairFI * const ccpTopPair = new mypairFI[numTopCcptoEvaluate];

     //for(n=0; n<numTopCcptoEvaluate; n++)
     for(n=0; n<numTopCcptoEvaluate; n++)
     {
	  //cout << "Volume: " << ccpVectorPairs[n].first << " Label: "<<  ccpVectorPairs[n].second << std::endl;
	  cout << " Label: "<<  ccpVectorPairs[n].second << " |  Volume: " << ccpVectorPairs[n].first ;
	  curLabelValue = ccpVectorPairs[n].second;
	  cout << " | Avg Prob.: "<< allConnectedComponents[curLabelValue].avgPredictProb << std::endl ;

	  // for each connected components
	  // step 1: extract the label
	  

	  ImageType3DUC::Pointer curCCPLabel = ImageType3DUC::New();
	  curCCPLabel->SetRegions(outRegion);
	  curCCPLabel->Allocate();
	  curCCPLabel->FillBuffer( 0 );

	  for(k=0; k<slice; k++)
	  {
	       for(j=0; j<height; j++)
	       {
		    for(i=0; i<width; i++)
		    {

			 ccpLabelValue = 0;
			 ccpIndex[0] = i;
			 ccpIndex[1] = j;
			 ccpIndex[2] = k;
			 ccpLabelValue = ccpPredictProb->GetPixel( ccpIndex );
			 
			 if(ccpLabelValue == curLabelValue)
			 {
			      curCCPLabel->SetPixel(ccpIndex, 1 );
			 }
		    }
	       }
	  }
	  
	  
	  // std::string test_save_1 = "test_CurrentCCPVolume.nii.gz";
	  // WriterType3DUC::Pointer writer = WriterType3DUC::New();
	  // writer->SetFileName( test_save_1 );
	  // writer->SetInput( curCCPLabel ) ;
	  // writer->Update();
	  

	  // Step 2: erosion
	  typedef itk::BinaryBallStructuringElement<PixelType3DUC, Dimension> StructuringElementType;
	  typedef itk::BinaryErodeImageFilter<ImageType3DUC, ImageType3DUC, StructuringElementType> ErodeFilterType;
	  
	  ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();

	  StructuringElementType structuringElement;

	  structuringElement.SetRadius( 1 ); // 3x3 structuring element
	  structuringElement.CreateStructuringElement();
	  binaryErode->SetKernel( structuringElement );
	  binaryErode->SetInput( curCCPLabel );
	  binaryErode->SetErodeValue( 1 );

	  // duplicator doesn't work!
	  /*
	  typedef itk::ImageDuplicator< ImageType3DUC > DuplicatorType;
	  DuplicatorType::Pointer duplicator = DuplicatorType::New();
	  duplicator->SetInputImage(binaryErode->GetOutput());
	  duplicator->Update();
	  ImageType3DUC::Pointer curCCPLabelErosion = duplicator->GetOutput();
	  */

	  // Note: the following code to save the rosion volume is to prevent
	  // segmentation fault error
	  
	  ImageType3DUC::Pointer curCCPLabelErosion = ImageType3DUC::New();
	  curCCPLabelErosion->SetRegions(outRegion);
	  curCCPLabelErosion->Allocate();
	  curCCPLabelErosion->FillBuffer( 0 );
	  curCCPLabelErosion = binaryErode->GetOutput();
	  std::string test_save_2 = "test_CurrentCCPVolume_erosion.nii.gz";
	  WriterType3DUC::Pointer writer2 = WriterType3DUC::New();
	  writer2->SetFileName( test_save_2 );
	  writer2->SetInput( curCCPLabelErosion ) ;
	  writer2->Update();


	  // step 3: original binary - erosion = the boundary voxels
	  //       : & count the boundary voxels
	  ImageType3DUC::Pointer curCCPLabelBoundary = ImageType3DUC::New();
	  curCCPLabelBoundary->SetRegions(outRegion);
	  curCCPLabelBoundary->Allocate();
	  curCCPLabelBoundary->FillBuffer( 0 );

	  ImageType3DUC::IndexType ccpIndexUC;
	  ImageType3DUC::PixelType ccpLabelValueUC_ori,ccpLabelValueUC_ero;
	  
	  int countNumBoundary = 0;
	  for(k=0; k<slice; k++)
	  {
	       for(j=0; j<height; j++)
	       {
		    for(i=0; i<width; i++)
		    {

			 ccpLabelValueUC_ori = 0;
			 ccpLabelValueUC_ero = 0;
			 ccpIndexUC[0] = i;
			 ccpIndexUC[1] = j;
			 ccpIndexUC[2] = k;
			 ccpLabelValueUC_ori = curCCPLabel->GetPixel( ccpIndexUC );
			 ccpLabelValueUC_ero = curCCPLabelErosion->GetPixel( ccpIndexUC );

			 if(ccpLabelValueUC_ori > 0 && ccpLabelValueUC_ero < 1)
			 {
			      curCCPLabelBoundary->SetPixel(ccpIndexUC, 1 );
			      countNumBoundary++;
			 }
		    }
	       }
	  }

	  
	  // std::string test_save_3 = "test_CurrentCCPVolume_boundary.nii.gz";
	  // WriterType3DUC::Pointer writer3 = WriterType3DUC::New();
	  // writer3->SetFileName( test_save_3 );
	  // writer3->SetInput( curCCPLabelBoundary ) ;
	  // writer3->Update();
	  

	  // compute the score for rank the top K connected components

	  curShapeScore = float(ccpVectorPairs[n].first) / float(countNumBoundary);
	  //curVolumeMultiplyProb = float(ccpVectorPairs[n].first) * allConnectedComponents[curLabelValue].avgPredictProb;
	  //curScoreCriteria = curVolumeMultiplyProb +			\
	  //     float(ccpVectorPairs[n].first) *  curShapeScore * curShapeScore * curShapeScore;
	  curScoreCriteria = curShapeScore;
	 
	  cout << " Num of boundary voxel for current CCP: "<< countNumBoundary << std::endl ;
	  cout << " Shape score for current CCP: "<< curShapeScore << std::endl ;
	  //cout << " Volume multiply Prob.: "<< curVolumeMultiplyProb << std::endl ;
	  cout << " Score for current CCP: "<< curScoreCriteria << std::endl ;
	  cout << " ---------------------- " << std::endl ;

	  // create new pair
	  ccpTopPair[n].first = curScoreCriteria;
	  ccpTopPair[n].second = ccpVectorPairs[n].second;
	  
     }

     // sort the ccpTopPair
     
     std::vector< pair<float, int> > ccpTopVectorPairs (ccpTopPair, ccpTopPair+numTopCcptoEvaluate);
     std::sort (ccpTopVectorPairs.begin(), ccpTopVectorPairs.end(), myComparatorforCcpFI);
     
     // create a volume for true candidates
     ImageType3DUC::Pointer trueCandidatesVolume = ImageType3DUC::New();
     trueCandidatesVolume->SetRegions(outRegion);
     trueCandidatesVolume->Allocate();
     trueCandidatesVolume->FillBuffer( 0 );

     // print newly sorted ccpTopPair
     ImageType3DUC::IndexType imageIndexUC;

     int countNumCCPsatisfyShapeCriteria = 0;

     cout << "=================== " << std::endl;
     for(n=0; n<numTopCcptoEvaluate; n++)
     {
	  cout << " Label: "<<  ccpTopVectorPairs[n].second << " |  Score: " << ccpTopVectorPairs[n].first;
	  
	  curLabelValue = ccpTopVectorPairs[n].second;

	  cout << " | Volume: "<< allConnectedComponents[curLabelValue-1].volumeSize << std::endl ;

	  if(ccpTopVectorPairs[n].first > thresholdScoreforTopCandidates)
	  {
	       cout << "This connected component's score is larger than the user-input threshold." << std::endl;
	       cout << std::endl;

	       countNumCCPsatisfyShapeCriteria++;
	       // add the label of this connected component to true candidates label volume
	       for(k=0; k<slice; k++)
	       {
		    for(j=0; j<height; j++)
		    {
			 for(i=0; i<width; i++)
			 {
			      ccpLabelValue = 0;
			      ccpIndex[0] = i;
			      ccpIndex[1] = j;
			      ccpIndex[2] = k;
			      ccpLabelValue = ccpPredictProb->GetPixel( ccpIndex );
			      
			      imageIndexUC[0] = i;
			      imageIndexUC[1] = j;
			      imageIndexUC[2] = k;

			      if(curLabelValue == ccpLabelValue)
			      {
				   trueCandidatesVolume->SetPixel(imageIndexUC, 1 );
			      }
			 }
		    }
	       }
	  }

     }
     cout << "=================== " << std::endl;

     // In case that countNumCCPsatisfyShapeCriteria == 0,
     // no single CCP satisfies the criteria so ask user

     bool stop_doing_activeLearning = false;

     int countNumThisCCPisTrueCandidate = 0;

     if(countNumCCPsatisfyShapeCriteria == 0)
     {

	  ImageType3DUC::Pointer PossibleCandidatesVolumeToAskUser = ImageType3DUC::New();
	  PossibleCandidatesVolumeToAskUser->SetRegions(outRegion);
	  PossibleCandidatesVolumeToAskUser->Allocate();
	  PossibleCandidatesVolumeToAskUser->FillBuffer( 0 );

	  bool this_CCP_is_true_candidate = false;

	  cout << ">>>===========User=Interaction============<<<" << std::endl;
	  for(n=0; n<numTopCcptoEvaluate; n++)
	  {

	       if(!stop_doing_activeLearning && countNumThisCCPisTrueCandidate<1)
	       {
		    cout << " Label: "<<  ccpTopVectorPairs[n].second << " |  Score: " << ccpTopVectorPairs[n].first;
		    curLabelValue = ccpTopVectorPairs[n].second;
		    cout << " | Volume: "<< allConnectedComponents[curLabelValue-1].volumeSize << std::endl ;

		    // extract the label volume of current label:
		    for(k=0; k<slice; k++)
		    {
			 for(j=0; j<height; j++)
			 {
			      for(i=0; i<width; i++)
			      {

				   ccpLabelValue = 0;
				   ccpIndex[0] = i;
				   ccpIndex[1] = j;
				   ccpIndex[2] = k;
				   ccpLabelValue = ccpPredictProb->GetPixel( ccpIndex );
			 
				   if(ccpLabelValue == curLabelValue)
				   {
					PossibleCandidatesVolumeToAskUser->SetPixel(ccpIndex, 1 );
				   }
				   else
				   {
					PossibleCandidatesVolumeToAskUser->SetPixel(ccpIndex, 0 );
				   }
			      }
			 }
		    }

		    // save current label volume
		    std::string test_save_4 = "AskUserCandidate.nii.gz";
		    WriterType3DUC::Pointer writer4 = WriterType3DUC::New();
		    writer4->SetFileName( test_save_4 );
		    writer4->SetInput( PossibleCandidatesVolumeToAskUser ) ;
		    writer4->Update();
	       
		    // ask user
		    std::string userStr;
		    std::string isTrueCandidate ("y");
		    std::string notTrueCandidate ("n");
		    std::string stopActiveLearning ("s");
		    std::cout << "Please look at the saved label volume 'AskUserCandidate.nii.gz' in current folder"<< std::endl;
		    std::cout << "If this CCP lesion? Yes or No (y/n), you can also choose stop (s) to stop active learning. \n";
		    std::cin >> userStr;
		    if(userStr.compare(isTrueCandidate) == 0)
		    {
			 // user input is yes
			 this_CCP_is_true_candidate = true;
			 countNumThisCCPisTrueCandidate++;
		    }
		    if(userStr.compare(stopActiveLearning) == 0)
		    {
			 // user input is yes
			 stop_doing_activeLearning = true;
		    }

		    ImageType3DUC::IndexType ccpIndexUC;
		    ImageType3DUC::PixelType ccpLabelValueUC;

		    if(this_CCP_is_true_candidate)
		    {
			 // User tells machine this CCP is true candidate
			 // add this label volume to true candidte volume
		    
			 for(k=0; k<slice; k++)
			 {
			      for(j=0; j<height; j++)
			      {
				   for(i=0; i<width; i++)
				   {
					ccpLabelValueUC = 0;
					ccpIndexUC[0] = i;
					ccpIndexUC[1] = j;
					ccpIndexUC[2] = k;
					ccpLabelValueUC = PossibleCandidatesVolumeToAskUser->GetPixel( ccpIndexUC );

					if(ccpLabelValueUC == 1)
					{
					     //std::cout << "In case user says 'y', add this CCP to true candidate."<< std::endl;
					     trueCandidatesVolume->SetPixel(ccpIndexUC, 1 );
					}
				   }
			      }
			 }
		    }
		    if(!this_CCP_is_true_candidate && !stop_doing_activeLearning)
		    {
			 // User tells machine this CCP is not lesion
			 // add this to trimap as 3 (HC_BG_NEW)

			 ImageType3DC::IndexType thisIndexDC;
			 
			 for(k=0; k<slice; k++)
			 {
			      for(j=0; j<height; j++)
			      {
				   for(i=0; i<width; i++)
				   {
					ccpLabelValueUC = 0;
					ccpIndexUC[0] = i;
					ccpIndexUC[1] = j;
					ccpIndexUC[2] = k;
					ccpLabelValueUC = PossibleCandidatesVolumeToAskUser->GetPixel( ccpIndexUC );

					thisIndexDC[0] = i;
					thisIndexDC[1] = j;
					thisIndexDC[2] = k;

					if(ccpLabelValueUC == 1)
					{
					     // HC_BG_NEW = 3
					     trimap->SetPixel(thisIndexDC, 3 );
					}
				   }
			      }
			 }
		    }// in case user tells this CCP is wrong.
	       }
	  }

	  cout << ">>>=========End=of=User=Interaction=========<<<" << std::endl;
     }

     if(!stop_doing_activeLearning || countNumThisCCPisTrueCandidate>0)
     {
	  ImageType3DUC::Pointer trueCandidatesVolumeErosion = ImageType3DUC::New();
	  trueCandidatesVolumeErosion->SetRegions(outRegion);
	  trueCandidatesVolumeErosion->Allocate();
	  trueCandidatesVolumeErosion->FillBuffer( 0 );

	  ImageType3DUC::Pointer trueCandidatesVolumeDilationMinusErosion = ImageType3DUC::New();
	  trueCandidatesVolumeDilationMinusErosion->SetRegions(outRegion);
	  trueCandidatesVolumeDilationMinusErosion->Allocate();
	  trueCandidatesVolumeDilationMinusErosion->FillBuffer( 0 );

	  morphologicalProcessingTrueCandidatesVolume(trueCandidatesVolume, trueCandidatesVolumeErosion, trueCandidatesVolumeDilationMinusErosion);

	  // use the trueCandidatesVolumeErosion & trueCandidatesVolumeDilationMinusErosion
	  // to modify the trimap; trueCandidatesVolumeErosion == 1, trimap = 1; foreground
	  // trueCandidatesVolumeDilationMinusErosion == 1, trimap = 0; uncertain
     
	  ImageType3DUC::IndexType curImageIndexUC;
	  ImageType3DUC::PixelType curImageValueUCerosion, curImageValueUCdelisionMinusErosion;

	  for(k=0; k<slice; k++)
	  {
	       for(j=0; j<height; j++)
	       {
		    for(i=0; i<width; i++)
		    {
			 curImageValueUCerosion = 0;
			 curImageValueUCdelisionMinusErosion = 0;
			 curImageIndexUC[0] = i;
			 curImageIndexUC[1] = j;
			 curImageIndexUC[2] = k;
			 curImageValueUCerosion = trueCandidatesVolumeErosion->GetPixel( curImageIndexUC );
			 curImageValueUCdelisionMinusErosion = trueCandidatesVolumeDilationMinusErosion->GetPixel( curImageIndexUC );

			 ccpIndex[0] = i;
			 ccpIndex[1] = j;
			 ccpIndex[2] = k;

			 if(curImageValueUCerosion > 0)
			 {
			      //std::cout << "set as foreground " << std::endl;
			      trimap->SetPixel(ccpIndex, 1 );
			 }

			 if(curImageValueUCdelisionMinusErosion > 0)
			 {
			      trimap->SetPixel(ccpIndex, 0 );
			 }
		    }
	       }
	  }

	  /*
	  std::cout << "-----------> Debug <----------- " << std::endl;
	  std::cout << "save current trimap to check the result. " << std::endl;
	  std::cout << "-----------> Debug <----------- " << std::endl;
	  WriterType3DC::Pointer writer = WriterType3DC::New();
	  //WriterType3DUC::Pointer writer = WriterType3DUC::New();
	  std::cout << "type a file name to save true candidate volume " << std::endl;
	  std::string save_name_for_checking_trimap;
	  std::cin >> save_name_for_checking_trimap;
	  writer->SetFileName( save_name_for_checking_trimap );
	  writer->SetInput( trimap ) ;
	  //writer->SetInput( trueCandidatesVolumeDilationMinusErosion ) ;
	  writer->Update();
	  std::cout << "This cin is for debug:" << std::endl;
	  int temptempVar = 0;
	  std::cin >> temptempVar;
	  */
	  

     }
     
     // free the memory space
     delete [] ccpTopPair;

     if(stop_doing_activeLearning)
     {
	  return false;
     }
     else
     {
	  return true;
     }
}



bool activeLearnCCP(ImageType3DF::Pointer predictProb,
		   ImageType3DI::Pointer alphaLabel,
		   ImageType3DC::Pointer trimap,
		   float thresholdPredictProb,
		   float thresholdScoreforTopCandidates,
		   int thresholdCCPVolume)
{
     ImageType3DI::Pointer binaryPredictProb = ImageType3DI::New();
     ImageType3DI::Pointer ccpPredictProb = ImageType3DI::New();

     // ImageType3DI::Pointer OutputImage = ImageType3DI::New();

     // Read information
     const ImageType3DI::SizeType sizeOfImage = alphaLabel->GetLargestPossibleRegion().GetSize();
     
     int width = sizeOfImage[0];
     int height = sizeOfImage[1];
     int slice = sizeOfImage[2];


     /*
     std::cout << "-----------> Debug <----------- " << std::endl;
     std::cout << "save alpha label and predicted prob. " << std::endl;
     std::cout << "-----------> Debug <----------- " << std::endl;
     std::cout << "type a file name to save alpha label " << std::endl;
     WriterType3DI::Pointer writerAlpha = WriterType3DI::New();
     std::string save_name_for_checking_alpha;
     std::cin >> save_name_for_checking_alpha;
     writerAlpha->SetFileName( save_name_for_checking_alpha );
     writerAlpha->SetInput( alphaLabel ) ;
     writerAlpha->Update();
     std::cout << "alpha label is saved! " << std::endl;
     std::cout << "type a file name to save predicted prob label " << std::endl;
     WriterType3DF::Pointer writerpredictedProb = WriterType3DF::New();
     std::string save_name_for_checking_predictedProb;
     std::cin >> save_name_for_checking_predictedProb;
     writerpredictedProb->SetFileName( save_name_for_checking_predictedProb );
     writerpredictedProb->SetInput( predictProb ) ;
     writerpredictedProb->Update();
     std::cout << "predicted prob is saved! " << std::endl;
     */

     // print the image size
     std::cout << "Volume Size" << std::endl;
     std::cout << "width: " <<  width << ", height: " << height << ", slice: " << slice << std::endl;

     // threshold the probability query volume
     ImageType3DI::RegionType outRegion = alphaLabel->GetLargestPossibleRegion();
     binaryPredictProb->SetRegions(outRegion);
     binaryPredictProb->Allocate();
     binaryPredictProb->FillBuffer( 0 );

     int i,j,k;

     ImageType3DF::IndexType pixelIndex;
     ImageType3DF::PixelType pixelValue;
     ImageType3DI::IndexType maskIndex;
     ImageType3DI::PixelType maskValue, binaryValue;

     ImageType3DC::IndexType pixelIndexDCtype;
     ImageType3DC::PixelType pixelValueDCtype;

     binaryValue = 1;

     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    pixelValue = 0;
		    maskValue = 0;

		    pixelIndex[0] = i;
		    pixelIndex[1] = j;
		    pixelIndex[2] = k;

		    maskIndex[0] = i;
		    maskIndex[1] = j;
		    maskIndex[2] = k;

		    pixelIndexDCtype[0] = i;
		    pixelIndexDCtype[1] = j;
		    pixelIndexDCtype[2] = k;

		    pixelValue = predictProb->GetPixel( pixelIndex );
		    maskValue = alphaLabel->GetPixel( maskIndex );
		    pixelValueDCtype = trimap->GetPixel( pixelIndexDCtype );

		    // remove the foreground area AND user rejected area
		    // alpha value is 1 or 0, 1: foreground; 0: background
		    if(maskValue < 1 && pixelValueDCtype != 3)
		    {
			 // the prediction probability is larger than the threshold
			 if(pixelValue >= thresholdPredictProb)
			 {
			      binaryPredictProb->SetPixel( maskIndex, binaryValue );
			 }
		    }
	       }
	  }
     }

     // preprocessing the binary volume of predicted probability
     std::cout << "Call the preprocessingBinaryPredicProb to do preprocessing. "  << std::endl;
     preprocessingBinaryPredicProb(binaryPredictProb);
     std::cout << "Preprocessing of binary volume is done! "  << std::endl;
     
     // connected component
     typedef itk::ConnectedComponentImageFilter <ImageType3DI, ImageType3DI >
	  ConnectedComponentImageFilterType;

     ConnectedComponentImageFilterType::Pointer connected =
	  ConnectedComponentImageFilterType::New ();
     connected->SetInput(binaryPredictProb);
     connected->Update();

     int numConnectedComponents = connected->GetObjectCount();

     std::cout << "Number of connected components: " << numConnectedComponents << std::endl;

     ccpPredictProb->SetRegions(outRegion);
     ccpPredictProb->Allocate();
     ccpPredictProb->FillBuffer( 0 );
     ccpPredictProb = connected->GetOutput();


     // initialize the struct for attributes of connected components
     //ConnectedCompInfor * const allConnectedComponents = (ConnectedCompInfor*)_alloca(numConnectedComponents * sizeof(ConnectedCompInfor)); // got error in linux
     ConnectedCompInfor * const allConnectedComponents = new ConnectedCompInfor[numConnectedComponents];

     int n;
     for(n=0; n<numConnectedComponents; n++)
     {
	  allConnectedComponents[n].labelValue = n+1;
	  allConnectedComponents[n].volumeSize = 0;
	  allConnectedComponents[n].avgPredictProb = 0;
	  allConnectedComponents[n].sumPredictProb = 0;
     }

     // initialize the struct for counting the CCP volume and PredictProb.
     // CountCCP * const countCCPforhere = (CountCCP*)_alloca(numConnectedComponents * sizeof(CountCCP));
     CountCCP * const countCCPforhere = new CountCCP[numConnectedComponents];

     for(n=0; n<numConnectedComponents; n++)
     {
	  countCCPforhere[n].countNumVoxels = 0;
	  countCCPforhere[n].curSumPredictProb = 0;
     }

     int countNumVoxel = 0;
     int curCcpLabelValue = 0;
     float sumPredictProb = 0;
     float avgPredictProbCcp = 0;
     ImageType3DI::IndexType ccpIndex;
     ImageType3DI::PixelType ccpLabelValue;

     float volumesize = slice*height*width;

     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    ccpLabelValue = 0;
		    ccpIndex[0] = i;
		    ccpIndex[1] = j;
		    ccpIndex[2] = k;
		    ccpLabelValue = ccpPredictProb->GetPixel( ccpIndex );

		    pixelValue = 0;
		    pixelIndex[0] = i;
		    pixelIndex[1] = j;
		    pixelIndex[2] = k;
		    pixelValue = predictProb->GetPixel( pixelIndex );

		    // check if this location is CCP
		    if(ccpLabelValue > 0)
		    {
			 countCCPforhere[ccpLabelValue - 1].countNumVoxels++;
			 countCCPforhere[ccpLabelValue - 1].curSumPredictProb = countCCPforhere[ccpLabelValue - 1].curSumPredictProb + pixelValue;
		    }

	       }// end i
	  }// end j
     }// end k

     /*
     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    ccpLabelValue = 0;
		    ccpIndex[0] = i;
		    ccpIndex[1] = j;
		    ccpIndex[2] = k;
		    ccpLabelValue = ccpPredictProb->GetPixel( ccpIndex );

		    pixelValue = 0;
		    pixelIndex[0] = i;
		    pixelIndex[1] = j;
		    pixelIndex[2] = k;
		    pixelValue = predictProb->GetPixel( pixelIndex );

		    // go through to check each CCP
		    for(n=0; n<numConnectedComponents; n++)
		    {
			 curCcpLabelValue = allConnectedComponents[n].labelValue;

			 if(ccpLabelValue == curCcpLabelValue)
			 {
			      countCCPforhere[n].countNumVoxels++;
			      countCCPforhere[n].curSumPredictProb = countCCPforhere[n].curSumPredictProb + pixelValue;
			 }//end if
		    }// end for n

		    std::cout << "Compute volume of each CCP: i : " << i << "; j : " << j << "; k : " << k << std::endl;
	       }// end i
	  }// end j
     }// end k
     */

     for(n=0; n<numConnectedComponents; n++)
     {
	  avgPredictProbCcp = 0;
	  avgPredictProbCcp = countCCPforhere[n].curSumPredictProb/countCCPforhere[n].countNumVoxels;

	  allConnectedComponents[n].volumeSize = countCCPforhere[n].countNumVoxels;
	  allConnectedComponents[n].avgPredictProb = avgPredictProbCcp;
	  allConnectedComponents[n].sumPredictProb = countCCPforhere[n].curSumPredictProb;
     }

     /*
     for(n=0; n<numConnectedComponents; n++)
     {
	  countNumVoxel = 0;
	  curCcpLabelValue = allConnectedComponents[n].labelValue;

	  sumPredictProb = 0;

	  for(k=0; k<slice; k++)
	  {
	       for(j=0; j<height; j++)
	       {
		    for(i=0; i<width; i++)
		    {
			 ccpLabelValue = 0;
			 ccpIndex[0] = i;
			 ccpIndex[1] = j;
			 ccpIndex[2] = k;
			 ccpLabelValue = ccpPredictProb->GetPixel( ccpIndex );

			 pixelValue = 0;
			 pixelIndex[0] = i;
			 pixelIndex[1] = j;
			 pixelIndex[2] = k;
			 pixelValue = predictProb->GetPixel( pixelIndex );

			 if(ccpLabelValue == curCcpLabelValue)
			 {
			      countNumVoxel++;
			      sumPredictProb = sumPredictProb + pixelValue;
			 }
		    }
	       }
	  }

	  avgPredictProbCcp = sumPredictProb/countNumVoxel;

	  allConnectedComponents[n].volumeSize = countNumVoxel;
	  allConnectedComponents[n].avgPredictProb = avgPredictProbCcp;
	  allConnectedComponents[n].sumPredictProb = sumPredictProb;

	  std::cout << "Compute the volume of each CCP: " << (float(n)/float(numConnectedComponents))*100 << "% is done." << std::endl;
     }

     */
     
     // print the struct of attributes of connected components: 
     
     for(n=0; n<numConnectedComponents; n++)
     {
	  std::cout << "Label value of ccp: " << allConnectedComponents[n].labelValue << std::endl;
	  std::cout << "Volume of ccp: " << allConnectedComponents[n].volumeSize << std::endl;
	  std::cout << "Avg predict prob. of ccp: " << allConnectedComponents[n].avgPredictProb << std::endl;
	  std::cout << "Sum predict prob. of ccp: " << allConnectedComponents[n].sumPredictProb << std::endl;
	  std::cout << "---------------" << std::endl;
     }
     

     std::cout << "Sort the CCP by volume size." << std::endl;
     //mypair * const ccpPair = (mypair*)_alloca(numConnectedComponents * sizeof(mypair));
     mypair * const ccpPair = new mypair[numConnectedComponents];

     for(n=0; n<numConnectedComponents; n++)
     {
	  ccpPair[n].first = allConnectedComponents[n].volumeSize;
	  ccpPair[n].second = allConnectedComponents[n].labelValue;
     }

     std::vector< pair<int, int> > ccpVectorPairs (ccpPair, ccpPair+numConnectedComponents);

     std::sort (ccpVectorPairs.begin(), ccpVectorPairs.end(), myComparatorforCcp);

     

     int numTopCCPtoEvaluate = 0;

     std::cout << "............................" << std::endl;
     for(n=0; n<numConnectedComponents; n++)
     {
	  if(ccpPair[n].first >= thresholdCCPVolume)
	  {
	       numTopCCPtoEvaluate++;
	       std::cout << "Add CCP label: " << ccpPair[n].second << " volume is: " << ccpPair[n].first << " to candidates list." << std::endl;
	  }
	   
     }
     std::cout << "............................" << std::endl;
     std::cout << "Total number of CCP candidates is: " << numTopCCPtoEvaluate << std::endl;
     std::cout << "............................" << std::endl;

     bool ruturnValue = false;

     if(numTopCCPtoEvaluate > 0)
     {
	  // do self training
	  ruturnValue = computeShapeScoreforCcpCandidates(ccpPredictProb, trimap, numTopCCPtoEvaluate, ccpVectorPairs, allConnectedComponents, thresholdScoreforTopCandidates);
	  
	  // free memory space before delete
	  delete [] allConnectedComponents;
	  delete [] countCCPforhere;
	  delete [] ccpPair;

	  return ruturnValue;
     }
     else
     {
	  // return 0 indicates that active learning should stop.
	  delete [] allConnectedComponents;
	  delete [] countCCPforhere;
	  delete [] ccpPair;
	  return false;
     }

     // debug code

     // save the output which is the modified trimap
     // std::cout << "-----------> Debug <----------- " << std::endl;
     // std::cout << "save current trimap to check the result. " << std::endl;
     // std::cout << "-----------> Debug <----------- " << std::endl;

	  
     // //WriterType3DC::Pointer writer = WriterType3DC::New();
     // WriterType3DUC::Pointer writer = WriterType3DUC::New();
     // std::cout << "type a file name to save true candidate volume " << std::endl;
     // std::string save_name_for_checking_trimap;
     // std::cin >> save_name_for_checking_trimap;
     // writer->SetFileName( save_name_for_checking_trimap );
     // //writer->SetInput( trimap ) ;
     // writer->SetInput( trueCandidatesVolume ) ;
     // writer->Update();

     // std::cout << "This cin is for debug:" << std::endl;
     // int temptempVar = 0;
     // std::cin >> temptempVar;

     // debug code

}


bool update_hardconstraints(const vnl_vector<double> & score_map, 
			   ImageType3DU::Pointer lindexPtr, 
			   const vnl_vector<unsigned> & alpha,
			   ImageType3DC::Pointer init_constraintPtr,
			   vnl_vector<unsigned> & hard_constraints,
			   const ParType & par)
{
     ImageType3DI::RegionType lindexRegion = lindexPtr->GetLargestPossibleRegion();

     // create alpha image. 
     ImageType3DI::Pointer alphaPtr = ImageType3DI::New();
     alphaPtr->SetRegions(lindexRegion);
     alphaPtr->Allocate();
     alphaPtr->FillBuffer( 0 ); // init to zero.

     // create query score image.
     ImageType3DF::Pointer scorePtr = ImageType3DF::New();
     scorePtr->SetRegions(lindexRegion);
     scorePtr->Allocate();
     scorePtr->FillBuffer( 0 ); // init to zero.

     IteratorType3DI alphaIt(alphaPtr, alphaPtr->GetLargestPossibleRegion());
     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());

     unsigned label = 0, sample_id = 0;
     for (alphaIt.GoToBegin(), lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt, ++alphaIt) {
	  if (lindexIt.Get() > 0) { // in mask.
	       sample_id = lindexIt.Get() - 1;
	       if (alpha[sample_id] == ALPHA_FG)
		    alphaIt.Set(1);
	       else 
		    alphaIt.Set(0);
	  }
     }

     // convert score to volume.
     IteratorType3DF scoreIt(scorePtr, scorePtr->GetLargestPossibleRegion());
     for (scoreIt.GoToBegin(), lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt, ++scoreIt) {
	  if (lindexIt.Get() > 0) { // in mask.
	       sample_id = lindexIt.Get() - 1;
	       scoreIt.Set(score_map[sample_id]);
	  }
     }


     bool cont_var = true;
     cont_var = activeLearnCCP(scorePtr, alphaPtr, init_constraintPtr, par.pred_th,  par.qscore_th, 100);     

     // update hard_constraint vector.
     load_constraints(lindexPtr, init_constraintPtr, hard_constraints, par);
     return cont_var;
}

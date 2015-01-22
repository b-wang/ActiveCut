// once the EM is done, we need to compute the posterior distribution of the
// class labels for each patch/voxel originally in the background region. We
// call it posterior, because we assume a MRF prior on the voxel labels x. Given
// x, we can compute the marginal P(y|x), where y is the voxel intensity. The
// marginal means we integrate out the component k when copute P(y|x). This is
// different with our previous solution just on the prior MRF. The posterior
// need to be in the range of [0, 1]. We use logistic model to convert the joint
// likelihood P(x, y) to the posterior. See page 197 of Bishop's book for that.

// the struct for attributes of connected components
typedef struct ConnectedCompInfor{
     int labelValue;
     int volumeSize;
     float avgPredictProb;
     float sumPredictProb;
} myConnectedCompInfor;

typedef struct CountCCP{
     int countNumVoxels;
     float curSumPredictProb;
} myCountCCP;

// the pair data structure for returning indices in sorting
typedef std::pair<int,int> mypair;
typedef std::pair<float,int> mypairFI;

// define unsigned char type for morphological processing
typedef unsigned char PixelType3DUC;
typedef itk::Image<PixelType3DUC , 3> ImageType3DUC;
typedef itk::ImageFileWriter< ImageType3DUC >  WriterType3DUC;


int logistic(vnl_vector<double> & score_map,
	     const vnl_matrix<double> & data,
	     vnl_sparse_matrix <double> & con_map,
	     const vnl_vector<unsigned> & alpha,
	     const ParType & par);


int logistic_init(vnl_vector<double> & score_map,
		  const vnl_matrix<double> & data,
		  vnl_sparse_matrix <double> & con_map, // can not set const due to a vnl bug.
		  const vnl_vector<unsigned> & alpha,
		  const ParType & par,
		  double smalleta,
		  unsigned n_iter);

void preprocessingBinaryPredicProb(ImageType3DI::Pointer binaryPredictProb);

void morphologicalProcessingTrueCandidatesVolume(ImageType3DUC::Pointer trueCandidatesVolume,
						 ImageType3DUC::Pointer trueCandidatesVolumeErosion, 
						 ImageType3DUC::Pointer trueCandidatesVolumeDilationMinusErosion);

bool computeShapeScoreforCcpCandidates (ImageType3DI::Pointer ccpPredictProb,
					ImageType3DI::Pointer trimap,
					const int& numTopCcptoEvaluate,
					std::vector< std::pair<int, int> >& ccpVectorPairs,
					ConnectedCompInfor allConnectedComponents[],
					const float& thresholdScoreforTopCandidates);


// return 0 indicates that active learning should stop.
bool activeLearnCCP(ImageType3DF::Pointer predictProb,
		   ImageType3DI::Pointer alphaLabel,
		   ImageType3DC::Pointer trimap,
		   float thresholdPredictProb,
		   float thresholdScoreforTopCandidates,
		   int thresholdCCPVolume);


// Read new object's labels, either by self training, or by act-learning, and
// update the constraints.
bool update_hardconstraints(const vnl_vector<double> & score_map, 
			   ImageType3DU::Pointer lindexPtr, 
			   const vnl_vector<unsigned> & alpha,
			   ImageType3DC::Pointer init_constraintPtr,
			   vnl_vector<unsigned> & hard_constraints,
			   const ParType & par);

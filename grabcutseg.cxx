#include <common.h>
#include <gmm.h>
#include <graphcuts.h>
#include <utility.h>
#include <loadfiles.h>
#include <query.h>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
     std::string input_image_file, init_image_file, alpha_file, tlink_image_file, nlink_image_file, mask_file, gmm_fg_file, gmm_bg_file, priorfg_file, priorbg_file, fg_label_file, bg_label_file, score_file, cand_file;
     
     // all parameters including GMM.
     ParType par;
     unsigned maxem = 0;;
     unsigned cov_type = 0;
     bool cov_shared = false;
     bool noactlearn = false;
     par.baseline = false;
     

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Grab cut segmentation of 3D multi-channel traumatic brain injury (TBI) images.")

	  ("ncompbg,b", po::value<unsigned>(&par.gmm_bg.n_comp)->default_value(4),
	   "number of Gaussian components in background.")

	  ("ncompfg,f", po::value<unsigned>(&par.gmm_fg.n_comp)->default_value(2),
	   "number of Gaussian components in foreground.")

	  ("maxem,m", po::value<unsigned>(&maxem)->default_value(30),
	   "Max number of EM iterations.")

	  ("kmeansruns,r", po::value<unsigned>(&par.kmeansruns)->default_value(5),
	   "number of kmeans runs.")

	  ("betaf", po::value<double>(&par.gmm_fg.beta)->default_value(1),
	   "smoothness (MRF) constraint for foreground gmm segmentation.")

	  ("betab", po::value<double>(&par.gmm_bg.beta)->default_value(1),
	   "Smoothness (MRF) constraint for background gmm segmentation.")

	  ("gamma,g", po::value<double>(&par.gamma)->default_value(1),
	   "alpha smoothness constraint.")

	  ("beta0", po::value<double>(&par.beta0)->default_value(1),
	   "An additional parameter in front of beta. ")

	  ("eta", po::value<double>(&par.eta)->default_value(1),
	   "Smoothness for the query score")

	  ("neighbors,n", po::value<unsigned>(&par.n_nbrs)->default_value(6),
	   "number of neighbors of a voxel for graphcuts N-Links. Legal values are 6, 18 and 26. ")

	  ("covtype", po::value<unsigned>(&cov_type)->default_value(1),
	   "Type of covariance matrix. 0 for full matrix, 1 for diagonal, 3 for identity.")

	  ("covshared", po::bool_switch(&cov_shared), 
	   "Whether all components of GMM has same covariance matrix. Default is no.")

	  ("data,d", po::value<std::string>(&input_image_file),
	   "Input all-channel image file. A 4D gipl or nii or nii.gz file.")

	  ("init,i", po::value<std::string>(&init_image_file),
	   "A 3D volume image giving the user's initial input. Inside the box will be init'd unkown, and outside will be init'd as background.")

	  ("priorfg", po::value<std::string>(&priorfg_file)->default_value("priorfg.nii.gz"),
	   "4D file for the foreground (bleeding/edema) prior probability. Give a zero image if no prior available.")

	  ("priorbg", po::value<std::string>(&priorbg_file)->default_value("priorbg.nii.gz"),
	   "4D file for the background (GM/WM/CSF) prior probability. Give a zero image if no prior available.")

	  ("lambdafg", po::value<double>(&par.lambda_fg)->default_value(1),
	   "the confidence on the foreground atlas. choose a lambda > 1 means more confidence on the atlas. ")

	  ("lambdabg", po::value<double>(&par.lambda_bg)->default_value(1),
	   "the confidence on the background atlas. choose a lambda > 1 means more confidence on the atlas. ")

	  ("mask,p", po::value<std::string>(&mask_file)->default_value("mask.nii.gz"),
	   "Mask file. Outside the brain is zero, inside is one.")

	  ("alphafile,a", po::value<std::string>(&alpha_file)->default_value("alpha.nii.gz"),
	   "The output segmentation file. Will be binary volume with same size as input data.")

	  ("gmmfg", po::value<std::string>(&gmm_fg_file)->default_value("gmmfg.nii.gz"),
	   "Foreground GMM posterior. ")

	  ("fglabel", po::value<std::string>(&fg_label_file)->default_value("gmmfg_label.nii.gz"),
	   "Foreground GMM segmentation label image.")

	  ("gmmbg", po::value<std::string>(&gmm_bg_file)->default_value("gmmbg.nii.gz"),
	   "Background GMM posterior.")

	  ("bglabel", po::value<std::string>(&bg_label_file)->default_value("gmmbg_label.nii.gz"),
	   "Background GMM segmentation label image.")

	  ("queryscores,q", po::value<std::string>(&score_file)->default_value("queryscores.nii.gz"),
	   "query scores file.")

	  ("cand,c", po::value<std::string>(&cand_file)->default_value("cand.nii.gz"),
	   "output candidate component for user to check.")

	  ("predth", po::value<float>(&par.pred_th)->default_value(0.75),
	   "Threshold used to get connected components from predictive probability.")

	  ("qscoreth", po::value<float>(&par.qscore_th)->default_value(1.3),
	   "Threshold for query score. Below that score, the active learning will stop.")

	  ("baseline", po::bool_switch(&par.baseline), 
	   "set baseline for baseline testing that allow background --> foreground and no active learning.")

	  ("noactlearn", po::bool_switch(&noactlearn), 
	   "If doing self-training and active learning after EM and graphcuts converges. default is yes.")

	  ("logmax", po::value<unsigned>(&par.logimax)->default_value(5),
	   "Max number of iterations in logistic() func.")

	  ("seed,s", po::value<unsigned>(&par.seed)->default_value(0),
	   "Seed for random number generator.")

	  ("verbose,v", po::value<unsigned>(&par.verbose)->default_value(0),
	   "verbose level in [0, 3]. ")
	  ;

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: grabcut [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cout << e.what() << "\n";
	  return 1;
     }    

     // read data images.
     ReaderType4DF::Pointer dataReader = ReaderType4DF::New();
     dataReader->SetFileName(input_image_file);
     dataReader->Update();
     ImageType4DF::Pointer dataPtr = dataReader->GetOutput();

     // load mask image. 
     ReaderType3DI::Pointer maskReader = ReaderType3DI::New();
     maskReader->SetFileName(mask_file);
     maskReader->Update();
     ImageType3DI::Pointer maskPtr = maskReader->GetOutput();

     // create a lindex map, which convert (x,y,z) to linear index n. Later the
     // lindex can be used for mask also.
     ImageType3DI::RegionType maskRegion = maskPtr->GetLargestPossibleRegion();     
     ImageType3DU::Pointer lindexPtr = ImageType3DU::New();
     lindexPtr->SetRegions(maskRegion);
     lindexPtr->Allocate();
     lindexPtr->FillBuffer( 0 ); // init to zero.

     // read user initialized image.
     ReaderType3DC::Pointer initReader = ReaderType3DC::New();
     initReader->SetFileName(init_image_file);
     initReader->Update();
     ImageType3DC::Pointer initPtr = initReader->GetOutput();

     // convert the image into a NxP matrix, where N is the number of
     // super-pixel, and P is the number of channels.
     vnl_matrix <double> data;
     load_data(dataPtr, maskPtr, data, lindexPtr);

     // fill these values since load_* func use them.
     par.n_channels = data.cols();
     par.n_samples = data.rows();

     // save the user input constraints, which include FG/BG/UNKNOWN.
     vnl_vector <unsigned> hard_constraints;
     hard_constraints.set_size(par.n_samples);
     hard_constraints.fill(0);
     load_constraints(lindexPtr, initPtr, hard_constraints, par);

     // load priors from atlas.
     vnl_matrix <double> priors_fg;
     vnl_matrix <double> priors_bg;
     load_priors(lindexPtr, priorfg_file, priors_fg, par);
     load_priors(lindexPtr, priorbg_file, priors_bg, par);

     // init the remaining part of par.
     par.gmm_fg.comp.resize(par.gmm_fg.n_comp);
     par.gmm_fg.pi.set_size(par.n_samples, par.gmm_fg.n_comp);
     par.gmm_fg.name = "GMM_FG";
     par.gmm_bg.comp.resize(par.gmm_bg.n_comp);
     par.gmm_bg.pi.set_size(par.n_samples, par.gmm_bg.n_comp);
     par.gmm_bg.name = "GMM_BG";

     for (unsigned comp_id = 0; comp_id < par.gmm_fg.n_comp; comp_id ++) {
	  par.gmm_fg.comp[comp_id].label = 0;
	  par.gmm_fg.comp[comp_id].numPts = 0;
	  par.gmm_fg.comp[comp_id].mu.set_size(par.n_channels);
	  par.gmm_fg.comp[comp_id].mu.fill(0);
	  par.gmm_fg.comp[comp_id].cov.set_size(par.n_channels, par.n_channels);
	  par.gmm_fg.comp[comp_id].cov.fill(0);
	  par.gmm_fg.comp[comp_id].inv_cov.set_size(par.n_channels, par.n_channels);
	  par.gmm_fg.comp[comp_id].inv_cov.fill(0);
     }
     for (unsigned comp_id = 0; comp_id < par.gmm_bg.n_comp; comp_id ++) {
	  par.gmm_bg.comp[comp_id].label = 0;
	  par.gmm_bg.comp[comp_id].numPts = 0;
	  par.gmm_bg.comp[comp_id].mu.set_size(par.n_channels);
	  par.gmm_bg.comp[comp_id].mu.fill(0);
	  par.gmm_bg.comp[comp_id].cov.set_size(par.n_channels, par.n_channels);
	  par.gmm_bg.comp[comp_id].cov.fill(0);
	  par.gmm_bg.comp[comp_id].inv_cov.set_size(par.n_channels, par.n_channels);
	  par.gmm_bg.comp[comp_id].inv_cov.fill(0);
     }

     vnl_vector <unsigned> alpha(par.n_samples, 0);
     // initialize alpha vector from user input. 
     for (unsigned sample_id = 0; sample_id < par.n_samples; sample_id ++) {
	  if (hard_constraints[sample_id] == HC_FG) alpha[sample_id] = ALPHA_FG;
	  else if (hard_constraints[sample_id] == HC_BG) alpha[sample_id] = ALPHA_BG;
	  else alpha[sample_id] = ALPHA_FG; // for unkown regions, init to FG. 
     }

     vnl_matrix<double> gmm_labels(par.n_samples, par.gmm_fg.n_comp > par.gmm_bg.n_comp? par.gmm_fg.n_comp: par.gmm_bg.n_comp, 0);

     // Kmeans segmentation on FG and BG.
     kmeans(data, alpha, gmm_labels, par, ALPHA_FG);
     kmeans(data, alpha, gmm_labels, par, ALPHA_BG);

     // permute labels in gmm_labels to align with priors.
     align_priors(gmm_labels, priors_fg, alpha, ALPHA_FG, par);
     align_priors(gmm_labels, priors_bg, alpha, ALPHA_BG, par);

     if (par.verbose >= 1) {
	  save_gmm_labelmap(gmm_labels, par, alpha, ALPHA_FG, lindexPtr, fg_label_file);
	  save_gmm_labelmap(gmm_labels, par, alpha, ALPHA_BG, lindexPtr, bg_label_file);
     }

     // do a M0 step, to estimate gmm parameters from kmeans labeling.
     gmm_mstep(data, gmm_labels, alpha, ALPHA_FG, par.gmm_fg, cov_shared, cov_type);
     gmm_mstep(data, gmm_labels, alpha, ALPHA_BG, par.gmm_bg, cov_shared, cov_type);
     update_pi(alpha, priors_fg, par.gmm_fg, par.lambda_fg);
     update_pi(alpha, priors_bg, par.gmm_bg, par.lambda_bg);
     print_par(par);

     // pre-compute the N-Links. 
     vnl_sparse_matrix<double> con_map;
     vnl_sparse_matrix<double> nlinks;
     build_adjmat(lindexPtr, data, con_map, par);
     build_nlinks(data, con_map, nlinks, par);

     unsigned alpha_change = 1e6;
     double LL_bg_old = 0, LL_bg = -1e10, LL_fg_old = 0, LL_fg = -1e10;
     bool LL_changed = true;
     unsigned em_iter = 0;
     vnl_vector<double> score_map(par.n_samples, 0);

     // outer loop for act-learning.
     unsigned cont_var = true;
     do {
	  alpha_change = 1e6;
	  LL_bg_old = 0;
	  LL_bg = -1e10;
	  LL_fg_old =0;
	  LL_fg = -1e10;
	  LL_changed = true;
	  em_iter = 0;

	  // innter loop: EM + graphcuts.
	  while( (LL_changed ) && em_iter < maxem)
	  {
	       em_iter ++;
	       gmm_estep(data, con_map, gmm_labels, alpha, ALPHA_FG, par.gmm_fg);
	       gmm_estep(data, con_map, gmm_labels, alpha, ALPHA_BG, par.gmm_bg);

	       // M step. Estimate parameters. 
	       gmm_mstep(data, gmm_labels, alpha, ALPHA_FG, par.gmm_fg, cov_shared, cov_type);
	       gmm_mstep(data, gmm_labels, alpha, ALPHA_BG, par.gmm_bg, cov_shared, cov_type);

	       LL_fg_old = LL_fg;
	       LL_bg_old = LL_bg;
	       LL_fg = eval_ll(data, alpha, par.gmm_fg, ALPHA_FG);
	       LL_bg = eval_ll(data, alpha, par.gmm_bg, ALPHA_BG);

	       LL_changed = ( fabs((LL_fg - LL_fg_old)/LL_fg_old) > EPS) || ( fabs((LL_bg - LL_bg_old)/LL_bg_old) > EPS) ; 

	       // graphcuts segmentation. Update alpha.
	       alpha_change = graphcuts(alpha, data, hard_constraints, nlinks,  par);
	       printf("EM iter %i. number of voxels with alpha changed: %i\n", em_iter, alpha_change);
	  } // inner loop finish.


	  // compute predictive score. First compute the score without MRF (eta
	  // = 0), then increase eta to refine the score map.
	  logistic_init(score_map, data, con_map, alpha, par, 0, 1); // iter = 1
	  logistic_init(score_map, data, con_map, alpha, par, par.eta/4, 2); 
	  logistic_init(score_map, data, con_map, alpha, par, par.eta/2, 4); 
	  logistic(score_map, data, con_map, alpha, par);

	  // ****** Active Learning begins here ************ //

	  // update hard_constraint image.
	  if (par.baseline || noactlearn) {
	       cont_var = false;
	  }
	  else {
	       cont_var = update_hardconstraints(score_map, lindexPtr, alpha, initPtr, hard_constraints, par);
	  }
     } while(cont_var);
     
     // test print GMM parameters:
     bool isFG = true;
     print_gmm_to_file(par.gmm_fg, isFG);
     print_gmm_to_file(par.gmm_bg, !isFG);

     save_gmm_posterior(gmm_labels, par, alpha, ALPHA_FG, lindexPtr, gmm_fg_file);
     save_gmm_posterior(gmm_labels, par, alpha, ALPHA_BG, lindexPtr, gmm_bg_file);
     save_gmm_labelmap(gmm_labels, par, alpha, ALPHA_FG, lindexPtr, fg_label_file);
     save_gmm_labelmap(gmm_labels, par, alpha, ALPHA_BG, lindexPtr, bg_label_file);

     save_alpha(par, alpha, lindexPtr, alpha_file);
     save_llmap(par, score_map, lindexPtr, score_file);
     
     return 0;
}

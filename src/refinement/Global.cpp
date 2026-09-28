#include "Global.h"
#ifdef OPENMP
#include <omp.h>
#endif

char*               Global::outputDirectory = NULL;			// output directory
std::string			Global::outputFileBasename;

char*               Global::posSequenceFilename = NULL;		// filename of positive sequence FASTA file
std::string			Global::posSequenceBasename;			// basename of positive sequence FASTA file
SequenceSet*        Global::posSequenceSet = NULL;			// positive sequence set
bool 		        Global::maskPosSequenceSet = false;		// mask motif patterns from positive sequence set

char*               Global::negSequenceFilename = NULL;		// filename of negative sequence FASTA file
std::string			Global::negSequenceBasename;			// basename of negative sequence FASTA file
SequenceSet*        Global::negSequenceSet = NULL;			// negative sequence set
bool				Global::negSeqGiven = false;			// a flag for the negative sequence given by users
bool                Global::genericNeg = false;             // flag for generating negative sequences based on generic 2nd-bgModel

// weighting options
char*               Global::intensityFilename = NULL;		// filename of intensity file (i.e. for HT-SELEX data)

char*				Global::alphabetType = NULL;			// alphabet type is defaulted to standard which is ACGT
bool                Global::ss = false;						// only search on single strand sequences

// initial model(s) options
char*				Global::initialModelFilename = NULL; 	// filename of initial model
std::string			Global::initialModelBasename;			// basename of initial model
std::string			Global::initialModelTag;				// tag for initializing the model
size_t				Global::maxPWM = std::numeric_limits<size_t>::max(); // number of init that are to be optimized
bool				Global::mops = false;					// learn MOPS model
bool				Global::zoops = true;					// learn ZOOPS model

// model options
size_t     			Global::modelOrder = 2;					// model order
std::vector<float> 	Global::modelAlpha( modelOrder+1, 1.f );// initial alphas
float				Global::modelBeta = 7.0f;				// alpha_k = beta x gamma^k for k > 0
float				Global::modelGamma = 3.0f;
std::vector<size_t>	Global::addColumns( 2 );				// add columns to the left and right of initial model
bool                Global::interpolate = true;             // calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mononucleotides
bool                Global::interpolateBG = true;			// calculate prior probabilities from lower-order probabilities
                                                            // instead of background frequencies of mononucleotides
// background model options
char*				Global::bgModelFilename = NULL;			// path to the background model file
bool				Global::bgModelGiven = false;			// flag to show if the background model is given or not
size_t				Global::bgModelOrder = 2;				// background model order, defaults to 2
std::vector<float>	Global::bgModelAlpha( bgModelOrder+1, 1.f );// background model alpha

// EM options
bool				Global::EM = false;						// flag to trigger EM learning
float				Global::q = 0.3f;						// prior probability for a positive sequence to contain a motif
bool 				Global::optimizeQ = false;				// optimize hyper-parameter q in EM algorithm
float               Global::f = 0.05f;                      // fraction of sequences to be masked

// CGS (Collapsed Gibbs sampling) options
bool				Global::CGS = false;					// flag to trigger Collapsed Gibbs sampling
bool				Global::noInitialZ = false;				// enable initializing z with one E-step
bool				Global::noAlphaOptimization = false;	// disable alpha optimization in CGS
bool				Global::GibbsMHalphas = false;			// enable alpha sampling in CGS using Gibbs Metropolis-Hastings
bool				Global::dissampleAlphas = false;		// enable alpha sampling in CGS using discretely sampling
bool				Global::noZSampling = false;			// disable q sampling in CGS
bool				Global::noQSampling = false;			// disable q sampling in CGS
bool				Global::debugAlphas = false;

// FDR options
bool				Global::FDR = false;					// triggers False-Discovery-Rate (FDR) estimation
size_t				Global::mFold = 1;						// number of negative sequences as multiple of positive sequences
size_t				Global::cvFold = 4;						// size of cross-validation folds
size_t				Global::sOrder = 2;						// the k-mer order for sampling negative sequence set

// motif occurrence options
bool                Global::scoreSeqset = false;            // write logOdds Scores of positive sequence set to disk
float 				Global::pvalCutoff = 0.0001f;			// score cutoff for printing log odds scores as motif hit

// printout options
bool                Global::verbose = false;
bool                Global::debugMode = false;              // debug-mode: prints out everything.
bool				Global::saveBaMMs = true;
bool				Global::savePRs = true;					// write the precision, recall, TP and FP
bool				Global::savePvalues = false;			// write p-values for each log odds score from sequence set
bool				Global::saveLogOdds = false;			// write the log odds of positive and negative sets to disk
bool				Global::saveInitialBaMMs = false;		// write out the initial model to disk
bool				Global::saveBgModel = false;			// write out the background model to disk
bool				Global::generatePseudoSet = false;		// test for alpha learning
std::mt19937		Global::rngx;

// flags for developers
bool			    Global::makeMovie = false;              // print out bamms in each iteration while optimizing
bool 				Global::B2 = false;
bool 				Global::B3 = false;
bool 				Global::B3prime = false;
bool                Global::advanceEM = false;

// option for openMP
size_t              Global::threads = 4;                   // number of threads to use

void Global::init( int nargs, char* args[] ){

	readArguments( nargs, args );

	Alphabet::init( alphabetType );

	// read in positive and negative sequence set
	posSequenceSet = new SequenceSet( posSequenceFilename, ss );
	negSequenceSet = new SequenceSet( negSequenceFilename, ss );

    // check if the input sequences are too few
    if( posSequenceSet->getSequences().size() < cvFold ){
        std::cerr << "Error: Input sequences are too few for training! \n" << std::endl;
        exit( 1 );
    }

	// optional: read in sequence intensities (header and intensity columns?)
	if( intensityFilename != 0 ){
		;// read in sequence intensity
	}
}

int Global::readArguments( int nargs, char* args[] ){

	/**
	 * read command line to get options
	 * process flags from user
	 */

	if( nargs >= 2 && ( !strcmp( args[1], "-h" ) || !strcmp( args[1], "--help" ) ) ){
		printHelp();
		exit( 0 );
	}

	if( nargs < 3 ) {
		std::cerr << "Error: Arguments are missing! \n" << std::endl;
		printHelp();
		exit( 1 );
	}

	// read in the output directory and create it
	outputDirectory = args[1];
	createDirectory( outputDirectory );

	// read in the positive sequence file
	posSequenceFilename = args[2];
    posSequenceBasename = baseName( posSequenceFilename );

	// read in options from the third argument on
	GetOpt::GetOpt_pp opt( nargs-2, args+2 );

	if( opt >> GetOpt::OptionPresent( 'h', "help" ) ){
		printHelp();
		exit( 0 );
	}

    // read in the basename of output file,
    // if not given, take the basename of input FASTA file
    if( opt >> GetOpt::OptionPresent( "basename" ) ){
        opt >> GetOpt::Option("basename", outputFileBasename);
    } else {
        outputFileBasename = posSequenceBasename;
    }
    // mask motif patterns from the positive sequence set
    opt >> GetOpt::OptionPresent( "maskPosSequenceSet", maskPosSequenceSet );

	// read in negative sequence file
	if( opt >> GetOpt::OptionPresent( "negSeqFile" ) ){
		negSeqGiven = true;
		opt >> GetOpt::Option( "negSeqFile", negSequenceFilename );
	} else {
	    negSequenceFilename = posSequenceFilename;
	}
	negSequenceBasename = baseName( negSequenceFilename );

    opt >> GetOpt::OptionPresent( "genericNeg", genericNeg );

    // Alphabet Type
	if( opt >> GetOpt::OptionPresent( "alphabet" ) ){
		opt >> GetOpt::Option( "alphabet", alphabetType );
	} else {
		alphabetType = strdup( "STANDARD" );
	}

	opt >> GetOpt::OptionPresent( "ss", ss );

	// for HT-SELEX data
	opt >> GetOpt::Option( "intensityFile", intensityFilename );

	// get initial model files
	std::string tag;
	if ( opt >> GetOpt::OptionPresent( "bindingSiteFile" ) ){
		opt >> GetOpt::Option( "bindingSiteFile", initialModelFilename );
		initialModelTag = "bindingsites";
	} else if ( opt >> GetOpt::OptionPresent( "PWMFile" ) ){
		opt >> GetOpt::Option( "PWMFile", initialModelFilename );
		initialModelTag = "PWM";
	} else if( opt >> GetOpt::OptionPresent( "BaMMFile" ) ){
		opt >> GetOpt::Option( "BaMMFile", initialModelFilename );
		initialModelTag = "BaMM";
	} else {
		fprintf( stderr, "Error: No initial model is provided.\n" );
		exit( 1 );
	}
	initialModelBasename = baseName( initialModelFilename );

	opt >> GetOpt::Option( "maxPWM", maxPWM );
	opt >> GetOpt::OptionPresent( "mops", mops );
	opt >> GetOpt::Option( "zoops", zoops );

	// model options
	opt >> GetOpt::Option( 'k', "order", modelOrder );

	if( opt >> GetOpt::OptionPresent( 'a', "alpha" ) ){
		modelAlpha.clear();
		opt >> GetOpt::Option( 'a', "alpha", modelAlpha );
		if( modelAlpha.size() != modelOrder+1 ){
			if( modelAlpha.size()  > modelOrder+1 ){
				modelAlpha.resize( modelOrder+1 );
			} else {
				modelAlpha.resize( modelOrder+1, modelAlpha.back() );
			}
		}
	} else {
		if( modelAlpha.size() != modelOrder+1 ){
			if( modelAlpha.size() > modelOrder+1 ){
				modelAlpha.resize( modelOrder+1 );
			} else {
				modelAlpha.resize( modelOrder+1, modelAlpha.back() );
			}
		}
		opt >> GetOpt::Option( 'b', "beta", modelBeta );
		opt >> GetOpt::Option( 'r', "gamma", modelGamma );
		if( modelOrder > 0 ){
			for( size_t k = 1; k < modelOrder+1; k++ ){
				// alpha = beta * gamma^k
				modelAlpha[k] = modelBeta * powf( modelGamma, ( float )k );
			}
		}
	}

	if( opt >> GetOpt::OptionPresent( "extend" ) ){
		addColumns.clear();
		opt >> GetOpt::Option( "extend", addColumns );
		if( addColumns.size() < 1 || addColumns.size() > 2 ){
			fprintf( stderr, "--extend format error.\n" );
			exit( 1 );
		}
		if( addColumns.size() == 1 )
			addColumns.resize( 2, addColumns.back() );
	} else {
		addColumns.at(0) = 0;
		addColumns.at(1) = 0;
	}

	// background model options
	if( opt >> GetOpt::Option( "bgModelFile", bgModelFilename ) ){
		bgModelGiven = true;
	}

	opt >> GetOpt::Option( 'K', "Order", bgModelOrder );

	if( opt >> GetOpt::OptionPresent( 'A', "Alpha" ) ){
		bgModelAlpha.clear();
		opt >> GetOpt::Option( 'A', "Alpha", bgModelAlpha );
		if( bgModelAlpha.size() != bgModelOrder+1 ){
			if( bgModelAlpha.size() > bgModelOrder+1 ){
				bgModelAlpha.resize( bgModelOrder+1 );
			} else {
				bgModelAlpha.resize( bgModelOrder+1, bgModelAlpha.back() );
			}
		}
	} else {
		if( bgModelAlpha.size() != bgModelOrder+1 ){
			if( bgModelAlpha.size() > bgModelOrder+1 ){
				bgModelAlpha.resize( bgModelOrder+1 );
			} else {
				bgModelAlpha.resize( bgModelOrder+1, bgModelAlpha.back() );
			}
		}
		if( bgModelOrder > 0 ){
			for( size_t k = 1; k < bgModelOrder+1; k++ ){
				bgModelAlpha[k] = 10.0f;
			}
		}
	}

	// EM options
	opt >> GetOpt::OptionPresent( "EM", EM );

	// CGS options
	if( opt >> GetOpt::OptionPresent( "CGS", CGS ) ){
		opt >> GetOpt::OptionPresent( "noInitialZ", noInitialZ );
		opt >> GetOpt::OptionPresent( "noAlphaOpti", noAlphaOptimization );
		opt >> GetOpt::OptionPresent( "GibbsMH", GibbsMHalphas );
		opt >> GetOpt::OptionPresent( "dissample", dissampleAlphas );
		opt >> GetOpt::OptionPresent( "noZSampling", noZSampling );
		opt >> GetOpt::OptionPresent( "noQSampling", noQSampling );
	}
	opt >> GetOpt::OptionPresent( "debugAlphas", debugAlphas );
	opt >> GetOpt::OptionPresent( "generatePseudoSet", generatePseudoSet );

	// saturation options
	opt >> GetOpt::Option( 'q', q );

    // masking options
	opt >> GetOpt::Option( 'f', f );

	// FDR options (-m and -s also apply to the background set used by --scoreSeqset)
	opt >> GetOpt::OptionPresent( "FDR", FDR );
	opt >> GetOpt::Option( 'm', "mFold", mFold );
	opt >> GetOpt::Option( 'n', "cvFold", cvFold );
	opt >> GetOpt::Option( 's', "sOrder", sOrder );
	// motif occurrence option
	opt >> GetOpt::OptionPresent( "scoreSeqset", scoreSeqset );
	opt >> GetOpt::Option( "pvalCutoff", pvalCutoff );

	// printout options
	opt >> GetOpt::OptionPresent( "verbose", verbose );
	opt >> GetOpt::OptionPresent( "debug", debugMode );
	opt >> GetOpt::OptionPresent( "saveBaMMs", saveBaMMs );
	opt >> GetOpt::OptionPresent( "saveInitialBaMMs", saveInitialBaMMs );
	opt >> GetOpt::Option( "savePRs", savePRs );
	opt >> GetOpt::OptionPresent( "savePvalues", savePvalues );
	opt >> GetOpt::OptionPresent( "saveLogOdds", saveLogOdds );
	opt >> GetOpt::OptionPresent( "saveBgModel", saveBgModel );

    // flags for developers
    opt >> GetOpt::OptionPresent( "makeMovie", makeMovie );
	opt >> GetOpt::OptionPresent( "optimizeQ", optimizeQ );
	opt >> GetOpt::OptionPresent( "B2", B2 );
	opt >> GetOpt::OptionPresent( "B3", B3 );
	opt >> GetOpt::OptionPresent( "B3prime", B3prime );
    opt >> GetOpt::OptionPresent( "advanceEM", advanceEM );

    // option for openMP
    opt >> GetOpt::Option( "threads", threads );
#ifdef OPENMP
    omp_set_num_threads( threads );
#endif

	// for remaining unknown options
	if( opt.options_remain() ){
		printHelp();
		std::cerr << "Oops! Unknown option(s) remaining... \n\n";
		exit( 1 );
	}

	return 0;
}

void Global::printStat(){

	std::cout << "Alphabet type is " << Alphabet::getAlphabet();
	std::cout << "\nGiven initial model is " << Global::initialModelBasename
              << ", BaMM order: " << Global::modelOrder
              << ", bgmodel order: " << Global::bgModelOrder;
	std::cout << "\nBaMM is learned from ";
	if( Global::ss ){
		std::cout << "single-stranded sequences.";
	} else {
		std::cout << "double-stranded sequences.";
	}

	// for positive sequence set
	std::cout << "\nGiven positive sequence set is "
              << Global::outputFileBasename << ".\n	"
              << Global::posSequenceSet->getSequences().size()
              << " sequences, max.length: " << Global::posSequenceSet->getMaxL()
              << ", min.length: " << Global::posSequenceSet->getMinL()
              << "\n	base frequencies:";
	for( size_t i = 0; i < Alphabet::getSize(); i++ ){
		std::cout << ' ' << Global::posSequenceSet->getBaseFrequencies()[i]
		          << "(" << Alphabet::getAlphabet()[i] << ")";
	}
    if( Global::advanceEM ){
        std::cout << "\n    " << Global::f * 100 << "% of the sequences are used for EM after masking.";
    }
	// for negative sequence set
	if( Global::negSeqGiven ){
		std::cout << "\nGiven negative sequence set is " << Global::negSequenceBasename
                  << ".\n	" << Global::negSequenceSet->getSequences().size()
                  << " sequences, max.length: " << Global::negSequenceSet->getMaxL()
                  << ", min.length: " << Global::negSequenceSet->getMinL()
                  << "\n	base frequencies:";
		for( size_t i = 0; i < Alphabet::getSize(); i++ )
			std::cout << ' ' << Global::negSequenceSet->getBaseFrequencies()[i]
					  << "(" << Alphabet::getAlphabet()[i] << ")";
	} else {
		std::cout << "\nThe background model is generated based on cond.prob of "
                  << Global::sOrder << "-mers.";
	}

	if( Global::FDR ){
		std::cout << "\nFolds for cross-validation (FDR estimation): "
                  << Global::cvFold;
	}
}

void Global::printHelp(){
	// Keep this text in sync with readArguments() and README.md.
	std::cout << R"(
SYNOPSIS
    BaMMmotif OUTDIR SEQFILE (--PWMFile FILE | --BaMMFile FILE | --bindingSiteFile FILE) [OPTIONS]

DESCRIPTION
    Refine motifs into higher-order Bayesian Markov models (BaMMs) with
    EM or collapsed Gibbs sampling, evaluate them by cross-validation and
    scan the input sequences for motif occurrences.

    OUTDIR   output directory (created if necessary)
    SEQFILE  positive sequences in FASTA format

INPUT
    --alphabet STRING       STANDARD (ACGT, default), METHYLC (ACGTM),
                            HYDROXYMETHYLC (ACGTH) or EXTENDED (ACGTMH)
    --ss                    search the given strand only (default: both strands;
                            not recommended for ChIP-seq data)
    --negSeqFile FILE       FASTA file with background sequences (reported in the
                            summary; BaMMmotif samples its own background set)
    --basename STRING       prefix of all output files (default: SEQFILE basename)

INITIAL MODELS (exactly one is required)
    --PWMFile FILE          position weight matrices in MEME format
    --BaMMFile FILE         a BaMM (.ihbcp); requires --bgModelFile when scoring
                            without optimisation
    --bindingSiteFile FILE  binding sites of equal length, one per line
    --maxPWM INT            number of motifs from --PWMFile to use (default: all)

MOTIF MODEL
    -k, --order INT         model order (default: 2)
    -a, --alpha FLOAT...    order-specific prior strengths; overrides -b and -r
                            (default: 1 for k = 0, beta * gamma^k for k > 0)
    -b, --beta FLOAT        beta in alpha_k = beta * gamma^k (default: 7)
    -r, --gamma FLOAT       gamma in alpha_k = beta * gamma^k (default: 3)
    --extend INT [INT]      add INT uniform positions to both ends, or the given
                            numbers to the left and right end (default: 0)
    -q FLOAT                prior fraction of sequences with a motif (default: 0.3)

BACKGROUND MODEL
    -K, --Order INT         background model order (default: 2)
    -A, --Alpha FLOAT...    prior strengths (default: 1 for k = 0, 10 for k > 0)
    --bgModelFile FILE      read the background model from a .hbcp file

OPTIMISATION (without --EM or --CGS the initial model is used as it is)
    --EM                    expectation maximisation
    --CGS                   collapsed Gibbs sampling (100 iterations)
      --noInitialZ          start from random motif positions instead of one E-step
      --noZSampling         do not sample motif positions
      --noQSampling         do not sample the motif fraction q
      --noAlphaOpti         do not optimise the prior strengths alpha
      --GibbsMH             sample alphas with Metropolis-Hastings
      --dissample           sample alphas from a discretised posterior

EVALUATION (cross-validation)
    --FDR                   estimate precision/recall by cross-validation and
                            write OUTDIR/<basename>_motif_<i>.zoops.stats
      -n, --cvFold INT      number of cross-validation folds (default: 4)
      -m, --mFold INT       background sequences per positive sequence; raised
                            automatically to give at least 5000 (default: 1)
      -s, --sOrder INT      k-mer order used to sample background sequences
                            (default: 2)
    --mops                  also evaluate the multiple-occurrences-per-sequence model
    --zoops BOOL            evaluate the zero-or-one-occurrence model (default: 1)

MOTIF OCCURRENCES
    --scoreSeqset           write motif occurrences with p-values to
                            OUTDIR/<basename>_motif_<i>.occurrence
    --pvalCutoff FLOAT      p-value cutoff for reported occurrences (default: 1e-4)

OUTPUT
    --saveBaMMs             also write k-mer counts (.counts) and motif positions
                            (.positions) of the optimised models
    --saveInitialBaMMs      write the initial models (_init_motif_<i>.ihbcp/.ihbp)
    --savePvalues           write p-values of the cross-validation scores
    --saveLogOdds           write log-odds scores of positive and background sets
    --savePRs BOOL          write .zoops.stats with --FDR (default: 1)
    --verbose               print progress of every iteration
    -h, --help              print this help

PERFORMANCE
    --threads INT           number of OpenMP threads (default: 4); results do not
                            depend on the number of threads

The background model (.hbcp/.hbp) and the final motif models (.ihbcp/.ihbp)
are always written to OUTDIR.
)";
}

void Global::destruct(){
    Alphabet::destruct();
    if( alphabetType ) 			free( alphabetType );
    if( posSequenceSet )	 	delete posSequenceSet;
    if( negSequenceSet ) 		delete negSequenceSet;
}

void Global::debug(){

}

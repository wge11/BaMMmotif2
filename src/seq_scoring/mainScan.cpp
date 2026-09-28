//
// Created by wanwan on 04.09.17.
// This function is used for scanning motif occurrences in a sequence set
//

#include "GScan.h"
#include "ScoreSeqSet.h"
#include "../init/MotifSet.h"
#include "../seq_generator/SeqGenerator.h"

int main( int nargs, char* args[] ) {

    /**
     * initialization
     */

	GScan::init( nargs, args );

    /**
     * Build up the background model
     */
    // use bgModel generated from input sequences when prediction is turned on
    BackgroundModel* bgModel;

    // use provided bgModelFile if initialized with bamm format
    if( GScan::initialModelTag == "BaMM" ) {
        if( GScan::bgModelFilename == NULL ) {
            std::cerr << "Error: No background model file provided for initial search motif!" << std::endl;
            exit( 1 );
        }
        // get background model from the given file
        bgModel = new BackgroundModel( GScan::bgModelFilename );
    } else {
        bgModel = new BackgroundModel( GScan::negSequenceSet->getSequences(),
                                       GScan::bgModelOrder,
                                       GScan::bgModelAlpha,
                                       GScan::interpolateBG,
                                       GScan::outputFileBasename );
        if( GScan::initialModelTag == "PWM" ){
            // this means that also the global motif order needs to be adjusted;
            GScan::modelOrder = 0;
        }
    }

    if(GScan::saveInitialModel){
        // save background model
        bgModel->write(GScan::outputDirectory, GScan::outputFileBasename);
    }

    /**
     * Initialize the model
     */
    MotifSet motif_set( GScan::initialModelFilename,
                        GScan::addColumns.at(0),
                        GScan::addColumns.at(1),
                        GScan::initialModelTag,
                        GScan::posSequenceSet,
                        bgModel->getV(),
                        GScan::bgModelOrder,
                        GScan::modelOrder,
                        GScan::modelAlpha,
                        GScan::maxPWM);

    std::vector<Sequence*> posSet = GScan::posSequenceSet->getSequences();
    /**
     * Filter out short sequences
     */
    removeShortSequences( posSet, motif_set.getMaxW() );

    /**
     * Sample negative sequence set based on s-mer frequencies
     */
    size_t minSeqN = 5000;
    // sample negative sequence set B1set based on s-mer frequencies
    // from positive training sequence set
    std::vector<std::unique_ptr<Sequence>> negSeqs;
    SeqGenerator negseq( posSet );
    if( posSet.size() >= minSeqN ){
        negSeqs = negseq.sample_bgseqset_by_fold( GScan::mFold );
    } else {
        negSeqs = negseq.sample_bgseqset_by_num( minSeqN, GScan::posSequenceSet->getMaxL() );
    }
    // negSeqs owns the sampled sequences; negset is a non-owning view
    std::vector<Sequence*> negset = rawPointers( negSeqs );

#pragma omp parallel for
    for( size_t n = 0; n < motif_set.getN(); n++ ) {
        // deep copy each motif in the motif set
        Motif *motif = new Motif( *motif_set.getMotifs()[n] );

        std::string fileExtension;
        if( GScan::initialModelTag == "PWM" ){
            fileExtension = "_motif_" + std::to_string( n+1 );
        }

        if( GScan::saveInitialModel ){
            // write out the foreground model
            motif->write( GScan::outputDirectory,
                          GScan::outputFileBasename + fileExtension );
        }

        // score negative sequence set
        ScoreSeqSet scoreNegSet( motif, bgModel, negset );
        scoreNegSet.calcLogOdds();
        std::vector<float> negScores = scoreNegSet.getAllMopsScores();

        // score positive sequence set
        // calculate p-values based on positive and negative scores
        ScoreSeqSet scorePosSet( motif, bgModel, posSet );
        scorePosSet.calcLogOdds();

        scorePosSet.calcPvalues( scorePosSet.getMopsScores(), std::move( negScores ) );

        scorePosSet.write( GScan::outputDirectory,
                           GScan::outputFileBasename + fileExtension,
                           GScan::pvalCutoff,
                           GScan::ss );

        delete motif;
    }

    if( bgModel ) delete bgModel;
    GScan::destruct();

    return 0;
}


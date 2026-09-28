/*
 * ScoreSeqSet.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: wanwan
 */

#include "ScoreSeqSet.h"
#include <float.h>		// -FLT_MAX

ScoreSeqSet::ScoreSeqSet( Motif* motif, BackgroundModel* bg, const std::vector<Sequence*>& seqSet ){

	motif_	= motif;
	bg_ 	= bg;
	seqSet_	= seqSet;
	Y_ 		= motif->getY();
    pval_is_calulated_ = false;
}

ScoreSeqSet::~ScoreSeqSet(){

}


void ScoreSeqSet::calcLogOdds(){

	/**
	 * store the log odds scores at all positions of each sequence
	 */

	size_t K = motif_->getK();
	size_t W = motif_->getW();
	size_t K_bg = ( bg_->getOrder() < K ) ? bg_->getOrder() : K;
	// pre-calculate log odds scores given motif and bg model
	motif_->calculateLogS( bg_->getV(), K_bg );
	float** s = motif_->getS();

	const size_t N  = seqSet_.size();
	const size_t YK = Y_[K+1];

	mops_scores_.assign( N, std::vector<float>() );
	zoops_scores_.assign( N, -FLT_MAX );
	z_.assign( N, 0 );

	// each sequence is scored independently, so this loop is safe to parallelise
#pragma omp parallel for schedule(dynamic, 64)
	for( size_t n = 0; n < N; n++ ){

		size_t 	LW1 = seqSet_[n]->getL() - W + 1;
		const size_t* kmer = seqSet_[n]->getKmer();
		float 	maxScore = -FLT_MAX;
        size_t  z_i = 0;

        std::vector<float>& scores = mops_scores_[n];
        scores.resize( LW1 );

		for( size_t i = 0; i < LW1; i++ ){
			float logOdds = 0.0f;
			for( size_t j = 0; j < W; j++ ){
				logOdds += s[kmer[i+j] % YK][j];
			}
			// take all the log odds scores for MOPS model:
			scores[i] = logOdds;

			// take the largest log odds score for ZOOPS model:
            if( logOdds > maxScore ){
                maxScore = logOdds;
                z_i = i;
            }
        }
		zoops_scores_[n] = maxScore;
        z_[n] = z_i;
	}
}

// compute p_values for motif scores based on negative sequence scores
void ScoreSeqSet::calcPvalues( const std::vector<std::vector<float>>& pos_scores, std::vector<float> neg_all_scores ){

	/**
	 * calculate P-values for motif occurrences
	 */

    size_t posN = seqSet_.size();
    size_t negN = neg_all_scores.size();
    mops_p_values_.assign( posN, std::vector<float>() );
    mops_e_values_.assign( posN, std::vector<float>() );

    float eps = 1.0e-5;

    // sort negative set scores in ascending order
    std::sort( neg_all_scores.begin(), neg_all_scores.end(), std::less<float>() );

    // Fit an exponential tail to the nTop highest negative scores. It is used
    // to extrapolate p-values for positive scores that exceed (almost) all
    // negative scores. Note that neg_all_scores is sorted in ascending order,
    // so the highest scores are at the end of the vector.
    size_t nTop = std::min<size_t>( 100, negN / 10 );
    nTop = std::max<size_t>( nTop, 1 );
    float S_ntop = neg_all_scores[negN - 1 - nTop];     // the nTop-th highest score

    // calculate the rate parameter lambda (mean excess over S_ntop)
    float lambda = 0.f;
	for( size_t n = 0; n < nTop; n++ ){
		lambda += ( neg_all_scores[negN - 1 - n] - S_ntop );
	}
	lambda = lambda / ( float )nTop;

#pragma omp parallel for
	for( size_t n = 0; n < seqSet_.size(); n++ ){

		size_t LW1 = seqSet_[n]->getL() - motif_->getW() + 1;
		mops_p_values_[n].reserve( LW1 );
		mops_e_values_[n].reserve( LW1 );

		for( size_t i = 0; i < LW1; i++ ){

            float Sl = pos_scores[n][i];
            // count the accumulated number of scores from the negative set up to rank l
            size_t FPl = std::distance( std::upper_bound( neg_all_scores.begin(), neg_all_scores.end(), Sl ),
                                        neg_all_scores.end() );

            float p_value;
			if( FPl == negN ){
                // when Sl is lower than the worst negative score:
				p_value = 1.f;
			} else if( FPl < 10 and fabs( lambda ) > eps ){
			    // when only few or no negatives are higher than S_l:
				p_value = float( nTop ) / ( float )negN * expf( - ( Sl - S_ntop ) / lambda );

			} else if( FPl == 0 ){
			    // Sl is higher than all negative scores and the exponential
			    // tail cannot be fitted: report the smallest resolvable p-value
			    p_value = 1.f / ( float )negN;

			} else {
				// when Sl_higher and Sl_lower can be defined:
				float SlHigher = neg_all_scores[negN-FPl-1];
				float SlLower = neg_all_scores[negN-FPl];
				p_value = ( ( float )FPl + ( SlHigher - Sl + eps ) / ( SlHigher - SlLower + eps ) ) / ( float )negN;
			}
            p_value = std::min( p_value, 1.f );
            mops_p_values_[n].push_back( p_value );
            mops_e_values_[n].push_back( p_value * ( float )posN );
		}
	}

    pval_is_calulated_ = true;
}

const std::vector<std::vector<float>>& ScoreSeqSet::getMopsScores() const {
	return mops_scores_;
}

const std::vector<float>& ScoreSeqSet::getZoopsScores() const {
	return zoops_scores_;
}

std::vector<float> ScoreSeqSet::getAllMopsScores() const {
	size_t total = 0;
	for( const std::vector<float>& scores : mops_scores_ ){
		total += scores.size();
	}
	std::vector<float> all;
	all.reserve( total );
	for( const std::vector<float>& scores : mops_scores_ ){
		all.insert( all.end(), scores.begin(), scores.end() );
	}
	return all;
}

void ScoreSeqSet::printLogOdds(){
    for( size_t n = 0; n < seqSet_.size(); n++ ){
        std::cout << "seq " << n << ":" << std::endl;
        std::cout << zoops_scores_[n] << '\t';
        for( size_t i = 0; i < mops_scores_[n].size(); i++ ){
            std::cout << mops_scores_[n][i] << '\t';
        }
        std::cout << std::endl;
    }
}

void ScoreSeqSet::write( char* odir, std::string basename, float pvalCutoff, bool ss ){

    /**
	 * save log odds scores in one flat file:
	 * basename.occurrence
	 */

    assert( pval_is_calulated_ );

	size_t 	end; 				// end of motif match

	std::string opath = std::string( odir )  + '/' + basename + ".occurrence";

	std::ofstream ofile( opath );

    // add a header to the results
    ofile << "seq\tlength\tstrand\tstart..end\tpattern\tp-value\te-value" << std::endl;

	for( size_t n = 0; n < seqSet_.size(); n++ ){
		size_t seqlen = seqSet_[n]->getL();
		if( !ss ){
			seqlen = ( seqlen - 1 ) / 2;
		}
		size_t LW1 = seqSet_[n]->getL() - motif_->getW() + 1;

        for( size_t i = 0; i < LW1; i++ ){

            // on double-stranded sequences, skip windows that span the
            // separator between the forward and the reverse-complement strand
            if( !ss and i <= seqlen and i + motif_->getW() > seqlen ){
                continue;
            }

			if( mops_p_values_[n][i] < pvalCutoff ){
                // >header:sequence_length
                ofile << seqSet_[n]->getHeader() << '\t' << seqlen << '\t';

				// start:end:score:strand:sequence_matching
				end = i + motif_->getW();

				ofile << ( ( i < seqlen ) ? '+' : '-' ) << '\t'
                      << i+1 << ".." << end << '\t';
				for( size_t m = i; m < end; m++ ){
					ofile << Alphabet::getBase( seqSet_[n]->getSequence()[m] );
				}
				ofile << '\t' << std::setprecision( 3 )
                      << mops_p_values_[n][i] << '\t'
                      << mops_e_values_[n][i] << std::endl;
			}
		}
	}

}

void ScoreSeqSet::writeLogOdds( char* odir, std::string basename, bool ss ){

    /**
	 * save log odds scores in one flat file:
	 * basename.logOddsZoops
	 */

    size_t 	end; 				// end of motif match

    std::string opath = std::string( odir )  + '/' + basename + ".logOddsZoops";

    std::ofstream ofile( opath );

    // add a header to the results
    ofile << "seq\tlength\tstrand\tstart..end\tpattern\tzoops_score" << std::endl;

    for( size_t n = 0; n < seqSet_.size(); n++ ){
        size_t seqlen = seqSet_[n]->getL();
        if( !ss ){
            seqlen = ( seqlen - 1 ) / 2;
        }

        // >header:sequence_length
        ofile << seqSet_[n]->getHeader() << '\t' << seqlen << '\t';

        // start:end:score:strand:sequence_matching
        end = z_[n] + motif_->getW();

        ofile << ( ( z_[n] < seqlen ) ? '+' : '-' ) << '\t'
              << z_[n]+1 << ".." << end << '\t';
        for( size_t m = z_[n]; m < end; m++ ){
            ofile << Alphabet::getBase( seqSet_[n]->getSequence()[m] );
        }
        ofile << '\t' << std::setprecision( 3 ) << zoops_scores_[n] << std::endl;


    }

}

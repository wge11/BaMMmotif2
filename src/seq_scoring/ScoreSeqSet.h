/*
 * ScoreSeqSet.h
 *
 *  Created on: Dec 13, 2016
 *      Author: wanwan
 */

#ifndef SCORESEQSET_H_
#define SCORESEQSET_H_

#include "../init/Motif.h"
#include "../init/BackgroundModel.h"

class ScoreSeqSet{
	/*
	 * This class is aimed for:
	 * scoring sequences from the given sequence set
	 * using the (learned) model and background model,
	 * find the occurrences of motifs on each sequence
	 * and output these sequences when the p-value
	 * is smaller than certain cutoff (default:0.0001)
	 */

public:

	ScoreSeqSet( Motif* motif, BackgroundModel* bg, const std::vector<Sequence*>& seqSet );
	~ScoreSeqSet();

	// log odds scores at every position of every sequence (MOPS) and the
	// best score per sequence (ZOOPS)
	void calcLogOdds();

	// p-values of the positive scores, estimated from the (unsorted) scores of
	// a negative/background sequence set; neg_all_scores is consumed (sorted)
	void calcPvalues( const std::vector<std::vector<float>>& pos_mops_scores,
	                  std::vector<float> neg_all_scores );

	const std::vector<std::vector<float>>& getMopsScores() const;
	const std::vector<float>&              getZoopsScores() const;
	std::vector<float>                     getAllMopsScores() const;	// all MOPS scores in one vector

	void write( char* odir, std::string basename, float pvalCutoff, bool ss );
    void writeLogOdds( char* odir, std::string basename, bool ss );
    void printLogOdds();

private:

	Motif* 							motif_;
	BackgroundModel* 				bg_;
	std::vector<Sequence*>			seqSet_;

    std::vector<float>				zoops_scores_;
	std::vector<std::vector<float>>	mops_scores_;
    std::vector<std::vector<float>> mops_p_values_;
    std::vector<std::vector<float>> mops_e_values_;
    std::vector<size_t>             z_;

    bool                            pval_is_calulated_;
	std::vector<size_t>				Y_;
};


#endif /* SCORESEQSET_H_ */

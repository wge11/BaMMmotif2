#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <algorithm>	// e.g. std::remove_if
#include <cstring>	// e.g. std::memcpy
#include <memory>	// e.g. std::unique_ptr
#include <vector>

#include <stdint.h>	// e.g. uint8_t
#include <math.h>

#include "Alphabet.h"

class Sequence{

public:

	Sequence( uint8_t* sequence,
				size_t L,
				std::string header,
				std::vector<size_t> Y,
				bool singleStrand = false );
	~Sequence();

	uint8_t*		getSequence();
	size_t			getL();
	std::string		getHeader();


	float 	        getIntensity();
	float 	        getWeight();
	size_t*			getKmer();		// get the value for 10-mer in the sequence

	void 	        setIntensity( float intensity );
	void 	        setWeight( float weight );

	void			print();		// print out sequences

private:
					// append the sequence's reverse complement to the sequence
	void 			appendRevComp( uint8_t* sequence, size_t L );

	uint8_t*		sequence_;		// sequence in alphabet encoding
	size_t			L_;				// sequence length
	std::string		header_;		// sequence header

	float			intensity_ = 0.0f;	// sequence intensity
	float			weight_ = 0.0f;	// sequence weight calculated from its intensity
	size_t*			kmer_;

	std::vector<size_t>			Y_;	// contains 1 at position 0
									// and the number of oligomers y for increasing order k at positions k+1
									// e.g.
									// alphabet size_ = 4: Y_ = 4^0 4^1 4^2 ... 4^15 < std::numeric_limits<int>::max()
									// limits the length of oligomers to 15 (and the order to 14)
};

inline size_t* Sequence::getKmer(){
	return kmer_;
}

// Remove all sequences shorter than minLength, keeping the order of the
// others. Returns the number of removed sequences.
inline size_t removeShortSequences( std::vector<Sequence*>& seqs, size_t minLength ){
	size_t before = seqs.size();
	seqs.erase( std::remove_if( seqs.begin(), seqs.end(),
	                            [minLength]( Sequence* seq ){ return seq->getL() < minLength; } ),
	            seqs.end() );
	return before - seqs.size();
}

// Non-owning view of a set of owned sequences.
inline std::vector<Sequence*> rawPointers( const std::vector<std::unique_ptr<Sequence>>& owned ){
	std::vector<Sequence*> view;
	view.reserve( owned.size() );
	for( const std::unique_ptr<Sequence>& seq : owned ){
		view.push_back( seq.get() );
	}
	return view;
}

#endif /* SEQUENCE_H_ */

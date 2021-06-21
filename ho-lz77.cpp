/*
 * Build co-BWT (using co-lex order) of the text and compute ho-LZ77 pairs
 */

#include "/scripts/code/dcheck.hpp"
#include <dynamic.hpp>
#include <unistd.h>
#include <vector>
#include <unordered_map>
#include <divsufsort.h>
#include "WaveletMatrix/DynamicWaveletMatrix/DynamicWaveletMatrix.hpp"
#include <sdsl/rmq_support.hpp>
#include "vector.hpp"



using namespace dyn;
using namespace std;

void help(){

	cout << "ho-lz77 [options]" << endl <<
	"Options:" << endl <<
	"-i <arg>   Input file (REQUIRED)" << endl <<
	"-t         Output in triples format (default: pairs)." << endl <<
	"-l         Output in locally-bitwise-optimal pair format (default: false)." << endl <<
	"-h         This help." << endl;
	 exit(0);

}
/*
 * number of bits required to write down x>0
 */
uint64_t bit_size(uint64_t x){

	return x==0 ? 1 : 64 - __builtin_clzll(x);

}

/*
 * compute the bit-size of the gamma encoding of x
 */
uint64_t gamma(uint64_t x){
    // std::cout << "computing gamma for x = " << x << " -> " <<  (2*bit_size(x) - 1) << std::endl;
    DCHECK_GT(x,0);
	return 2*bit_size(x) - 1;

}

/*
 * compute the bit-size of the delta encoding of x
 */
uint64_t delta(uint64_t x){

	auto bits = bit_size(x);//bits needed to encode x

	return gamma(bits) + bits -1;

}

/*
 * compute how many bits would the parse take in practice
 */
uint64_t compute_gamma_bit_complexity(Vector<pair<int64_t,uint64_t> > & parse){

	uint64_t x = 0;

	for(auto p:parse){

		uint64_t off = p.first<0?-p.first:p.first;//absolute value
		uint64_t len = p.second;

		x++;//1 bit for the sign of off

		x += gamma(off);
		x += gamma(len);

	}

	return x;

}

uint64_t compute_delta_bit_complexity(Vector<pair<int64_t,uint64_t> > & parse){

	uint64_t x = 0;

	for(auto p:parse){

		uint64_t off = p.first<0?-p.first:p.first;//absolute value
		uint64_t len = p.second;

		x++;//1 bit for the sign of off

		x += delta(off);
		x += delta(len);

	}

	return x;

}



uint64_t compute_succinct_bit_complexity(Vector<pair<int64_t,uint64_t> > & parse){

	uint64_t x = 0;
	double sum_off = 0;
	double sum_len = 0;

	for(auto p:parse){

		uint64_t off = p.first<0?-p.first:p.first;//absolute value
		uint64_t len = p.second;

		sum_off += off;
		sum_len += len;

	}

	double z = parse.size();

	return uint64_t(z*(5 + log2(sum_off/z) + log2(sum_len/z))) + 1;

}

uint64_t compute_Huff_bit_complexity(Vector<pair<int64_t,uint64_t> > & pairs, Vector<char> trails){

	uint64_t max_encoded_off = 3200000;//max offset that we encode with Huffman. The others use delta(0).delta(offset) (0 receives the longest code)
	uint64_t max_encoded_len = 100000;//max phrase length that we encode with Huffman. The others use delta(0).delta(len) (0 receives the longest code)

	Vector<uint64_t> off_freq(max_encoded_off+1);
	Vector<uint64_t> len_freq(max_encoded_len+1);
	Vector<uint64_t> char_freq(256);

	uint64_t z = pairs.size();

	for(uint64_t i=0;i<z;++i){

		//absolute value
		pairs[i].first = pairs[i].first < 0 ? -pairs[i].first : pairs[i].first;

		if(pairs[i].first <= max_encoded_off) off_freq[uint64_t(pairs[i].first)]++;
		if(pairs[i].second <= max_encoded_len) len_freq[pairs[i].second]++;
		char_freq[uint8_t(trails[i])]++;

	}

	double off_sum = 0;
	double len_sum = 0;
	double char_sum = 0;

	for(auto x : off_freq) off_sum += x;
	for(auto x : len_freq) len_sum += x;
	for(auto x : char_freq) char_sum += x;

	Vector<pair<uint64_t,double> > off_f(max_encoded_off+1);
	Vector<pair<uint64_t,double> > len_f(max_encoded_len+1);
	Vector<pair<uint64_t,double> > char_f(256);

	bool p50=false;
	bool p75=false;
	bool p95=false;


	double F = 0;//cumulative
	for(uint64_t i=0;i<max_encoded_off+1;++i){

		off_f[i].first = i;
		off_f[i].second = double(off_freq[i])/off_sum;

		F += off_f[i].second;

		if(F>=0.5 and not p50){ cout << "50% percentile : " << i << endl; p50=true;}
		if(F>=0.75 and not p75){ cout << "75% percentile : " << i << endl; p75=true;}
		if(F>=0.95 and not p95){ cout << "95% percentile : " << i << endl; p95=true;}

	}

	for(uint64_t i=0;i<max_encoded_len+1;++i){

		len_f[i].first = i;
		len_f[i].second = double(len_freq[i])/len_sum;

	}

	for(uint64_t i=0;i<256;++i){

		char_f[i].first = i;
		char_f[i].second = double(char_freq[i])/char_sum;

	}

	auto off_encoding = alphabet_encoder(off_f);
	auto len_encoding = alphabet_encoder(len_f);
	auto char_encoding = alphabet_encoder(char_f);

	/*
	cout << "Huffman lengths for offsets: " << endl;
	for(uint64_t i = 0;i<=500;++i)
		cout << i << ": " << off_encoding.encode(i).size() << endl;

	cout << "Huffman lengths for lengths: " << endl;
	for(uint64_t i = 0;i<=500;++i)
		cout << i << ": " << len_encoding.encode(i).size() << endl;

	cout << "Huffman lengths for trail chars: " << endl;
	for(uint64_t c = 0;c<256;++c)
		cout << char(c) << ": " << char_encoding.encode(uint8_t(c)).size() << endl;*/


	//uint64_t x = (max_encoded_off+2 + max_encoded_len + 2 + 256)*64; //the model
	uint64_t x = 0;

	for(uint64_t i=0;i<z;++i){

		x++; // sign of offset

		if(pairs[i].first <= max_encoded_off)
			x += off_encoding.encode(uint64_t(pairs[i].first)).size();
		else
			x += off_encoding.encode(0).size() + delta(uint64_t(pairs[i].first));

		if(pairs[i].second <= max_encoded_len)
			x += len_encoding.encode(pairs[i].second).size();
		else
			x += len_encoding.encode(0).size() + delta(pairs[i].second);

		x += char_encoding.encode(uint8_t(trails[i])).size();


	}

	return x;

}


/*
 * entropy of a vector
 */
double entropy(Vector<uint64_t> & V){

	unordered_map<uint64_t, double> freq;
	double n = V.size();

	for(auto x : V)	freq[x]=0;
	for(auto x : V)	freq[x]++;

	double H = 0;

	for(auto p : freq){

		auto f = p.second/n;

		H -= f*log2(f);

	}

	return H;

}



constexpr bool empty_interval(const pair<uint64_t, uint64_t>& interval) {
	return interval.first >= interval.second;
}

/*
 * input: bwt, index on L column, number of steps L
 * output: perform L LF steps from index and return resulting index
 *
 * note: we actually call function FL (not LF) because the underlying BWT
 * is logically suffix-sorted
 */
uint64_t jump_back(wt_bwt & bwt, uint64_t index, uint64_t L){

	for(uint64_t i = 0; i<L; ++i)
		index = bwt.FL(index);

	return index;

}

inline int64_t get_off(	wt_bwt & bwt,
						pair<uint64_t,uint64_t> & range,
						uint64_t index,
						uint64_t phrase_len){

	//assign a default high value to previous/next occurrence of the phrase in co-lex order;
	//if they are not initialized later, the high value will not be chosen since
	//we minimize the distance with the current prefix in co-lex order
	int64_t pred_occ = bwt.size()*2;
	int64_t succ_occ = bwt.size()*2;

	//if previous position is valid, compute it
	if(index>0 && index-1 >= range.first && index-1 < range.second){

		pred_occ = jump_back(bwt, index-1, phrase_len);

	}

	//if following position is valid, compute it
	if(index >= range.first && index < range.second){

		succ_occ = jump_back(bwt, index, phrase_len);

	}

	//at least one of the two must be initialized since there must be a previous
	//occurrence of the phrase!
	assert(pred_occ < bwt.size() or succ_occ < bwt.size());

	//co-lex position of current prefix
	int64_t current_prefix_pos = bwt.get_terminator_position();

	//the previous occurrence cannot be the current text prefix!
	assert(pred_occ != current_prefix_pos and succ_occ != current_prefix_pos);

	//absolute distances
	int64_t abs_dist_pred = current_prefix_pos > pred_occ ? current_prefix_pos - pred_occ : pred_occ - current_prefix_pos;
	int64_t abs_dist_succ = current_prefix_pos > succ_occ ? current_prefix_pos - succ_occ : succ_occ - current_prefix_pos;

	return abs_dist_pred < abs_dist_succ ? current_prefix_pos - pred_occ : current_prefix_pos - succ_occ;

}

inline void output_phrase(	wt_bwt & bwt,
							Vector<pair<int64_t, uint64_t> > & LZ77k,
							pair<uint64_t,uint64_t> & prev_range,
							pair<uint64_t,uint64_t> & range,
							uint64_t & index,
							string & phrase){

	range = prev_range;

	//closest co-lex positions before and
	//after wt.get_terminator_position()
	//that are followed by the current phrase

	//assign a default high value to previous/next occurrence of the phrase in co-lex order;
	//if they are not initialized later, the high value will not be chosen since
	//we minimize the distance with the current prefix in co-lex order
	int64_t pred_occ = bwt.size()*2;
	int64_t succ_occ = bwt.size()*2;

	//if previous position is valid, compute it
	if(index>0 && index-1 >= range.first && index-1 < range.second){

		pred_occ = jump_back(bwt, index-1, phrase.length());

	}

	//if following position is valid, compute it
	if(index >= range.first && index < range.second){

		succ_occ = jump_back(bwt, index, phrase.length());

	}

	//at least one of the two must be initialized since there must be a previous
	//occurrence of the phrase!
	assert(pred_occ < bwt.size() or succ_occ < bwt.size());

	//co-lex position of current prefix
	int64_t current_prefix_pos = bwt.get_terminator_position();

	//the previous occurrence cannot be the current text prefix!
	assert(pred_occ != current_prefix_pos and succ_occ != current_prefix_pos);

	//absolute distances
	int64_t abs_dist_pred = current_prefix_pos > pred_occ ? current_prefix_pos - pred_occ : pred_occ - current_prefix_pos;
	int64_t abs_dist_succ = current_prefix_pos > succ_occ ? current_prefix_pos - succ_occ : succ_occ - current_prefix_pos;

	int64_t occ = abs_dist_pred < abs_dist_succ ? current_prefix_pos - pred_occ : current_prefix_pos - succ_occ;

	//create new phrase
	LZ77k.push_back({occ,phrase.length()});

	//extend BWT with the characters in the phrase
	for(auto c : phrase)
		bwt.extend(uint8_t(c));

	//erase phrase content
	phrase = string();

	//re-initialize range to full range
	range = bwt.get_full_interval();

	//re-initialize index of current text prefix to (new) position of terminator
	index = bwt.get_terminator_position();

}

inline void output_phrase(	wt_bwt & bwt,
							Vector<pair<int64_t, uint64_t> > & LZ77k,
							Vector<char> & trail_chars,
							pair<uint64_t,uint64_t> & prev_range,
							pair<uint64_t,uint64_t> & range,
							uint64_t & index,
							string & phrase){

	range = prev_range;

	//closest co-lex positions before and
	//after wt.get_terminator_position()
	//that are followed by the current phrase

	//assign a default high value to previous/next occurrence of the phrase in co-lex order;
	//if they are not initialized later, the high value will not be chosen since
	//we minimize the distance with the current prefix in co-lex order
	int64_t pred_occ = bwt.size()*2;
	int64_t succ_occ = bwt.size()*2;

	//if previous position is valid, compute it
	if(index>0 && index-1 >= range.first && index-1 < range.second){

		pred_occ = jump_back(bwt, index-1, phrase.length()-1);

	}

	//if following position is valid, compute it
	if(index >= range.first && index < range.second){

		succ_occ = jump_back(bwt, index, phrase.length()-1);

	}

	//at least one of the two must be initialized since there must be a previous
	//occurrence of the phrase!
	assert(pred_occ < bwt.size() or succ_occ < bwt.size());

	//co-lex position of current prefix
	int64_t current_prefix_pos = bwt.get_terminator_position();

	//the previous occurrence cannot be the current text prefix!
	assert(pred_occ != current_prefix_pos and succ_occ != current_prefix_pos);

	//absolute distances
	int64_t abs_dist_pred = current_prefix_pos > pred_occ ? current_prefix_pos - pred_occ : pred_occ - current_prefix_pos;
	int64_t abs_dist_succ = current_prefix_pos > succ_occ ? current_prefix_pos - succ_occ : succ_occ - current_prefix_pos;

	int64_t occ = abs_dist_pred < abs_dist_succ ? current_prefix_pos - pred_occ : current_prefix_pos - succ_occ;

	//create new phrase
	LZ77k.push_back({occ,phrase.length()-1});
	trail_chars.push_back(phrase[phrase.length()-1]);

	//c remains the previous one that we read, since it is not part of the phrase

	//extend BWT with the characters in the phrase
	for(auto c : phrase)
		bwt.extend(uint8_t(c));

	//erase phrase content
	phrase = string();

	//re-initialize range to full range
	range = bwt.get_full_interval();

	//re-initialize index of current text prefix to (new) position of terminator
	index = bwt.get_terminator_position();

}

std::ifstream::pos_type filesize(const char* filename) {
	std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
	return in.tellg(); 
}

std::string read_file(const set<char>& alphabet, const char* filePath) {
	const size_t length = filesize(filePath);
	std::string buffer(length+alphabet.size(), ' ');
	size_t offset = 0;
	for (auto rit = alphabet.rbegin(); rit != alphabet.rend(); rit++){
		buffer[offset++] = *rit;
	}
	std::ifstream t(filePath);
	t.read(&buffer[offset], length); 
	return buffer;
}


inline void output_phrase_bitoptimal(wt_bwt & bwt,
							Vector<pair<int64_t, uint64_t> > & LZ77k,
							pair<uint64_t,uint64_t> & range,
							uint64_t & index,
							string & phrase,
							DynamicWaveletMatrix& dynwt,
							Vector<int> isa,
							size_t text_position
							){

	//closest co-lex positions before and
	//after wt.get_terminator_position()
	//that are followed by the current phrase

	//assign a default high value to previous/next occurrence of the phrase in co-lex order;
	//if they are not initialized later, the high value will not be chosen since
	//we minimize the distance with the current prefix in co-lex order
	int64_t pred_occ = bwt.size()*2;
	int64_t succ_occ = bwt.size()*2;

	//if previous position is valid, compute it
	if(index>0 && index-1 >= range.first && index-1 < range.second){

		pred_occ = jump_back(bwt, index-1, phrase.length());

	}

	//if following position is valid, compute it
	if(index >= range.first && index < range.second){

		succ_occ = jump_back(bwt, index, phrase.length());

	}

	//at least one of the two must be initialized since there must be a previous
	//occurrence of the phrase!
	assert(pred_occ < bwt.size() or succ_occ < bwt.size());

	//co-lex position of current prefix
	int64_t current_prefix_pos = bwt.get_terminator_position();

	//the previous occurrence cannot be the current text prefix!
	assert(pred_occ != current_prefix_pos and succ_occ != current_prefix_pos);

	//absolute distances
	int64_t abs_dist_pred = current_prefix_pos > pred_occ ? current_prefix_pos - pred_occ : pred_occ - current_prefix_pos;
	int64_t abs_dist_succ = current_prefix_pos > succ_occ ? current_prefix_pos - succ_occ : succ_occ - current_prefix_pos;

	int64_t occ = abs_dist_pred < abs_dist_succ ? current_prefix_pos - pred_occ : current_prefix_pos - succ_occ;

	//create new phrase
	LZ77k.push_back({occ,phrase.length()});

	//extend BWT with the characters in the phrase
	for(size_t i = 0; i < phrase.length(); ++i) {
		const char c = phrase[i];
		bwt.extend(uint8_t(c));
		const uint64_t inserted_index = bwt.get_terminator_position();
		dynwt.insert(inserted_index, isa[text_position+i]);
	}


	//erase phrase content
	phrase = string();

	//re-initialize range to full range
	range = bwt.get_full_interval();

	//re-initialize index of current text prefix to (new) position of terminator
	index = bwt.get_terminator_position();

}

/** 
 * @brief Kasai's LCP array construction algorithm
 * 
 * @param text 
 * @param sa text's suffix array, maybe a subarray of the actual suffix array
 * @param isa sa's inverse
 * 
 * @return LCP-array of text wrt. sa
 */
template<typename vektor_type, typename string_type>
vektor_type create_lcp(const string_type& text, const vektor_type& sa, const vektor_type& isa) {
	vektor_type lcp(sa.size());
	lcp[0] = 0;
	size_t h = 0;
	for(size_t i = 0; i < lcp.size(); ++i) {
		if(isa[i] == 0) continue;
		const size_t j = sa[ isa[i] -1 ];
		while(text[i+h] == text[j+h]) ++h;
		lcp[isa[i]] = h;
		h = h > 0 ? h-1 : 0;
	}
	return lcp;
}

void run_lz77(const string& filePath){
	cout << "Computing plain overlapping LZ77" << endl;
	wt_bwt bwt; //! Huffman-encoded BWT
	set<char> alphabet;

	{
		ifstream in(filePath);
		auto F = get_frequencies(in);

		for(auto f : F) if(f.second>0) alphabet.insert(f.first);

		bwt = wt_bwt(F); 
	}
	const std::string text = read_file(alphabet, filePath.c_str());
	Vector<int> sa(text.length(), 0);

	divsufsort(reinterpret_cast<const uint8_t*>(text.c_str()), sa.data(), text.length());

	Vector<int> isa ( text.length(), 0 );
	for(size_t i = 0; i < text.length(); ++i) {
	    isa[sa[i]] = i;
	}
	auto lcp = create_lcp(text, sa, isa);

	constexpr size_t UNDEF = -1;

	Vector<size_t> psvsa( text.length(), UNDEF );
	for(size_t i = 1; i < text.length(); ++i) {
	    size_t j=i-1;
	    while(j != UNDEF && sa[j] > sa[i]) {
		j = psvsa[j];
	    }
	    psvsa[i] = j;
	}
	if(text.length() < 100) {
	    std::cout << "text: " << text << std::endl;
	}

	Vector<size_t> nsvsa( text.length(), UNDEF );
	
	for(size_t offset = 2; offset <= text.length(); ++offset) {
	    const size_t i = text.length() - offset;
	    size_t j=i+1;
	    while(j != UNDEF && sa[j] > sa[i]) {
		j = nsvsa[j];
	    }
	    nsvsa[i] = j;
	}
	sdsl::rmq_succinct_sct<true> rmqlcp(&lcp);
	Vector<size_t> lengths;
	Vector<size_t> positions;

	for(size_t cur_text_position = alphabet.size(); cur_text_position < text.length();) {
	    const size_t cur_sa_idx = isa[cur_text_position]; //! where in SA is cur_text_position?
	    DCHECK(psvsa[cur_sa_idx] != UNDEF || nsvsa[cur_sa_idx] != UNDEF);

	    if(nsvsa[cur_sa_idx] == UNDEF) {
		const size_t prev_lcp = lcp[rmqlcp(psvsa[cur_sa_idx]+1, cur_sa_idx)];
		const size_t prev_pos = sa[psvsa[cur_sa_idx]];
		const size_t prev_distance = cur_text_position - prev_pos;
		DCHECK_GT(prev_distance, 0);
		DCHECK_LT(prev_pos, cur_text_position);
		positions.push_back(prev_distance);
		lengths.push_back(prev_lcp);
	    } else if(psvsa[cur_sa_idx] == UNDEF) {
		const size_t next_lcp = lcp[rmqlcp(cur_sa_idx+1, nsvsa[cur_sa_idx])];
		const size_t next_pos = sa[nsvsa[cur_sa_idx]];
		const size_t next_distance = cur_text_position - next_pos;
		DCHECK_LT(next_pos, cur_text_position);
		lengths.push_back(next_lcp);
		positions.push_back(next_distance);
		DCHECK_GT(next_distance, 0);
	    } else {
		const size_t prev_lcp = lcp[rmqlcp(psvsa[cur_sa_idx]+1, cur_sa_idx)];
		const size_t prev_pos = sa[psvsa[cur_sa_idx]];
		const size_t prev_distance = cur_text_position - prev_pos;

		const size_t next_lcp = lcp[rmqlcp(cur_sa_idx+1, nsvsa[cur_sa_idx])];
		const size_t next_pos = sa[nsvsa[cur_sa_idx]];
		const size_t next_distance = cur_text_position - next_pos;

		DCHECK_GT(prev_lcp+next_lcp, 0);
		DCHECK_LT(static_cast<size_t>(sa[nsvsa[cur_sa_idx]]), cur_text_position);
		lengths.push_back(std::max(prev_lcp, next_lcp));
		positions.push_back(prev_lcp > next_lcp ? prev_distance : next_distance);
	    }
	    DCHECK_GT(lengths.back(), 0);
	    DCHECK_GT(positions.back(), 0);
	    if(text.length() < 100) {
		std::cout << "pos: " << cur_text_position << " len: " << lengths.back() << " ref: " << positions.back() << std::endl;
	    }
	    cur_text_position += lengths.back();
	}
	{
	    size_t costs = 0;
	    for(size_t factor_index = 0; factor_index < positions.size(); ++factor_index) {
		costs += gamma(positions[factor_index]) + gamma(lengths[factor_index]); 
	    }
	    cout << "LZ77 factorization: # = " << positions.size() << " gamma = " << costs << endl;
	}
}


void run_pairs_bitoptimal(const string& filePath){

	cout << "Computing the BIT-OPTIMAL parse in format (occ,len)" << endl << endl;

	wt_bwt bwt; //! Huffman-encoded BWT
	set<char> alphabet;

	{
		ifstream in(filePath);
		auto F = get_frequencies(in);

		for(auto f : F) if(f.second>0) alphabet.insert(f.first);

		bwt = wt_bwt(F); 
	}
	const std::string text = read_file(alphabet, filePath.c_str());
	Vector<int> sa(text.length(), 0);

	divsufsort(reinterpret_cast<const uint8_t*>(text.c_str()), sa.data(), text.length());

	Vector<int> isa ( text.length(), 0 );
	for(size_t i = 0; i < text.length(); ++i) {
		isa[sa[i]] = i;
	}
	auto lcp = create_lcp(text, sa, isa);

	sdsl::rmq_succinct_sct<true> rmqlcp(&lcp);

	//! maps from ISA^R[i] to ISA[i], where i is a text position
	//! ISA^R[i] are changing values, but ISA[.] is fix
	DynamicWaveletMatrix dynwt(text.length());
	// dynwt.insert(0, 0); //! the $ of ISA^R[i] maps to ISA[text.length()] = 0

	// prepend the alphabet to the text
	{
	    size_t i = 0;
	    auto rit = alphabet.rbegin();
	    while(rit != alphabet.rend()) {
		const uint64_t index = bwt.get_terminator_position();
		dynwt.insert(index, isa[i]);
		bwt.extend(uint8_t(*rit));
		++i; ++rit;
	    }
	}

    struct Factor {
	size_t m_len = 0;
	size_t m_rsa_idx = 0;
	size_t m_sa_idx = 0;
    };


	// process the text

	Vector<Vector<Factor>> factor_dag;//! factors[1..n] stores in each entry a list of candidate factors, each entry corresponds to a text position
	factor_dag.resize(text.size() - alphabet.size());

	Vector<size_t> text_rsa_idx(text.size() - alphabet.size()); //! stores for each text position i its SA index of T[1..i]^R

	auto get_candidate =  [&] (const size_t begin_pos, const size_t end_pos, const size_t cur_sa_idx) -> Factor {
		const auto pred_sa_idx = dynwt.prevValue(begin_pos, end_pos, 0, cur_sa_idx); // search for the predecessor in dynwt of cur_sa_idx
		const auto succ_sa_idx = dynwt.nextValue(begin_pos, end_pos, cur_sa_idx+1, text.length());

		const size_t pred_lcp_idx = pred_sa_idx == NOTFOUND ? 0 : rmqlcp(pred_sa_idx+1, cur_sa_idx);
		const size_t succ_lcp_idx = succ_sa_idx == NOTFOUND ? 0 : rmqlcp(cur_sa_idx+1, succ_sa_idx);
		if(lcp[pred_lcp_idx] == 0 && lcp[succ_lcp_idx] == 0) return {};

		const size_t best_sa_idx = lcp[pred_lcp_idx] < lcp[succ_lcp_idx] ? succ_sa_idx : pred_sa_idx;
		const size_t best_lcp = std::max(lcp[pred_lcp_idx], lcp[succ_lcp_idx]);
		const size_t best_rsa_idx = dynwt.select(best_sa_idx, 1);
		// printf("best {textpos : %d, saidx: %lu, lcp: %lu, rsaidx: %lu}\n", sa[best_sa_idx], best_sa_idx, best_lcp, best_rsa_idx );
		return { best_lcp, best_rsa_idx, best_sa_idx };
	};


	for(size_t cur_text_position = alphabet.size(); cur_text_position < text.length(); ++cur_text_position) {
	    DCHECK_LT(cur_text_position - alphabet.size(), factor_dag.size());
	    auto& factorlist = factor_dag[cur_text_position - alphabet.size()];
	    const uint64_t cur_rsa_idx = bwt.get_terminator_position(); //! BWT position where we will insert the current read character
	    text_rsa_idx[cur_text_position - alphabet.size()] = cur_rsa_idx;

	    const size_t cur_sa_idx = isa[cur_text_position]; //! where in SA is cur_text_position?
	    // printf("text {pos : %lu, sa_idx: %lu, rsa_idx: %lu, char : %c}\n", cur_text_position, cur_sa_idx, cur_rsa_idx, text[cur_text_position]);

	    // printf("Left\n");
	    {//leftwards in SA(T^R)
		const Factor longest_left_factor = get_candidate(0, cur_rsa_idx, cur_sa_idx);
		if(longest_left_factor.m_len > 0) {
		    factorlist.push_back(longest_left_factor);
		    size_t max_lcp = 0; //! now only store if we make a progress in lcp since we go further away from cur_rsa_idx meaning storing the distance gets more costly
		    for(size_t left_window_size = 1, old_window_size = 0; left_window_size < cur_rsa_idx - longest_left_factor.m_rsa_idx; left_window_size*=2) {
			if(cur_rsa_idx <= left_window_size) { break; }
			auto ret = get_candidate(cur_rsa_idx-left_window_size, cur_rsa_idx-old_window_size, cur_sa_idx);
			if(ret.m_len > max_lcp) {
			    factorlist.emplace_back(ret);
			    max_lcp = ret.m_len;
			}
			old_window_size = left_window_size;
		    }
		}
	    }

	    // printf("Right\n");
	    {//rightwards in SA(T^R)
		const auto longest_right_factor = get_candidate(cur_rsa_idx, dynwt.size, cur_sa_idx);
		if(longest_right_factor.m_len > 0) {
		    factorlist.push_back(longest_right_factor);
		    size_t max_lcp = 0; //! now only store if we make a progress in lcp since we go further away from cur_rsa_idx meaning storing the distance gets more costly
		    for(size_t right_window_size = 1, old_window_size = 0; right_window_size < longest_right_factor.m_rsa_idx - cur_rsa_idx; right_window_size*=2) {
			if(cur_rsa_idx+right_window_size >= dynwt.size) { break; }
			auto ret = get_candidate(cur_rsa_idx+old_window_size, cur_rsa_idx+right_window_size, cur_sa_idx);
			if(ret.m_len > max_lcp) {
			    factorlist.emplace_back(ret);
			    max_lcp = ret.m_len;
			}
			old_window_size = right_window_size;
		    }
		}
	    }
	    DCHECK_GT(factorlist.size(), 0); //! must have at least one factor stored per entry
	    if(cur_text_position % 10000 == 0) {
		printf("Completed: %f ; #List[%lu/%lu] = %lu \r", cur_text_position*100.0/text.size(), cur_text_position, text.size(), factorlist.size());
		fflush(stdout);
	    }

	    //! prepare for the next text position
	    dynwt.insert(cur_rsa_idx, isa[cur_text_position]); //! fill in our character such that `cur_rsa_idx` is our target position in RSA space
	    bwt.extend(uint8_t(text[cur_text_position]));
	}


	auto get_offset = [&] (const size_t dag_index, const size_t factor_index) 
	{ 
	    return static_cast<int64_t>(text_rsa_idx[dag_index]) - 
		    static_cast<int64_t>(factor_dag[dag_index][factor_index].m_rsa_idx); 
		   //! gives the offset of a factor
	};

	auto get_gamma_cost = [&] (const size_t dag_index, const size_t factor_index) 
	{ 
	    return 
		gamma((uint64_t)std::abs(get_offset(dag_index, factor_index))+1) //+1 since the offset can be zero
		+1 // for the sign bit of the offset
		+gamma(factor_dag[dag_index][factor_index].m_len);
	};


	{//greedy factorization
	    size_t count = 0;
	    size_t costs = 0;
	    for(size_t pos = 0; pos < factor_dag.size(); ) {
		auto& factorlist = factor_dag[pos];
		DCHECK(!factorlist.empty());
		std::sort(factorlist.begin(), factorlist.end(), [&] (const Factor& a, const Factor& b) { if(a.m_len == b.m_len) return a.m_rsa_idx < b.m_rsa_idx; return a.m_len < b.m_len; });
		costs += get_gamma_cost(pos, factorlist.size()-1);
		++count;
		pos += factorlist.back().m_len;
	    }
	    cout << "greedy factorization: # = " << count << " gamma = " << costs << endl;
	}

	constexpr size_t UNDEF = -1;
	Vector<size_t> parents(factor_dag.size()+1, UNDEF); //! stores the starting position of the factor that reaches T[i]
	Vector<size_t> parent_index(factor_dag.size()+1, UNDEF); //! stores the index of the parent factor within factor_dag
	Vector<size_t> distances(factor_dag.size()+1, UNDEF); //! stores the costs from T[1] to T[i]
	std::vector<bool> processed(factor_dag.size()+1, false); //! stores the costs from T[1] to T[i]
	distances[0] = 0;
	parents[0] = 0;

	auto cmp = [&distances](size_t left, size_t right) { return distances[left] > distances[right]; }; //! smallest distance to the top. 
	std::priority_queue<size_t, Vector<size_t>, decltype(cmp)> queue(cmp);
	queue.push(0);
	while(!queue.empty()) {
	    const auto pos = queue.top();
	    queue.pop();
	    const auto& is_processed = processed[pos];
	    if(is_processed == true) { continue; }
	    if(pos == factor_dag.size()) { continue; } //! we do not process the rightmost node/last text position
	    processed[pos] = true;
	    const auto& factorlist = factor_dag[pos];
	    for(size_t factor_index = 0; factor_index < factorlist.size(); ++factor_index) {
		const auto cost = get_gamma_cost(pos, factor_index);
		const auto target = pos + factorlist[factor_index].m_len;
		DCHECK_LT(pos, target);
		DCHECK_LT(target, factor_dag.size()+1);
		if(distances[pos] + cost < distances[target]) {
		    distances[target] = distances[pos]+cost;
		    parents[target] = pos;
		    parent_index[target] = factor_index;
		    queue.push(target);
		}
		
	    }
	}
    
	const size_t bitoptimal_factor_count = [&]() -> size_t {
	    size_t count = 0;
	    for(size_t i = factor_dag.size(); i != 0;) {
		DCHECK_NE(parents[i], UNDEF);
		i = parents[i];
		++count;
	    }
	    cout << "bitoptimal factorization: # = " << count << " gamma = " << distances.back() << endl;
	    return count;
	}();
#ifndef NDEBUG
	Vector<Factor> bitoptimal_factors;
	bitoptimal_factors.reserve(bitoptimal_factor_count);
	for(size_t i = factor_dag.size(); i != 0;) {
	    DCHECK_NE(parents[i], UNDEF);
	    bitoptimal_factors.push_back(factor_dag[parents[i]][parent_index[i]]);
	    i = parents[i];
	}
	std::reverse(bitoptimal_factors.begin(), bitoptimal_factors.end());

	std::string decomp_text(text.length(), 0);

	size_t decomp_offset = 0;
	for (auto rit = alphabet.rbegin(); rit != alphabet.rend(); rit++) {
		decomp_text[decomp_offset++] = *rit;
	}
	for(size_t factor_index = 0; factor_index < bitoptimal_factors.size(); ++factor_index) {
	    const auto& factor = bitoptimal_factors[factor_index];
	    for(size_t factor_c = 0; factor_c < factor.m_len; ++factor_c) {
		decomp_text[decomp_offset++] = decomp_text[sa[factor.m_sa_idx]+factor_c];
	    }
	}
	DCHECK(decomp_text == text);
#endif//NDEBUG

		// uint64_t read_char = 0;
		// if(in.is_open()) {
        //
		// 	if(in.good()) {
		// 		in.get(c);//get character
		// 	}
        //
		// 	while(in.good()) {
        //
		// 		prev_range = range;
		// 		range = bwt.LF(range, uint8_t(c));
        //
		// 		if(empty_interval(range)){//end of phrase
        //
		// 			range = prev_range;
		// 			output_phrase_bitoptimal(bwt,LZ77k,range,index,phrase,dynwt,isa,read_char);
        //
		// 		}else{
        //
		// 			index = bwt.LF(index,uint8_t(c));
		// 			phrase += c;
		// 			in.get(c);//get next character
		// 			read_char++;
        //
		// 			if(read_char%100000 == 0){
        //
		// 				cout << "read " << (read_char+1) << " characters." << endl;
        //
		// 			}
        //
		// 		}
        //
		// 	}
        //
		// }
        //
	// 	//last phrase has not been output
	// 	if(phrase.length()>0){
    //
	// 		prev_range = range;
	// 		output_phrase_bitoptimal(bwt,LZ77k,range,index,phrase,dynwt,isa,read_char);
    //
	// 	}
    //
	// 	if(!in.eof() && in.fail())
	// 		cout << "Error reading " << filePath << endl;
    //
	// }


	// cout << "factorization: " << endl;
	// for(auto p : LZ77k){
        //
	// 	auto d = 1+delta(p.first<0?-p.first:p.first);
	// 	auto l = delta(p.second);
	// 	cout << "off = " << p.first << " (" << d  << "  bits), len = " << p.second << " (" << l << " bits). Tot: " << (double(d+l)/double(p.second)) << " bit/symbol)" << endl;
        //
	// }
        //
	// auto N = bwt.size() -1;//file length
        //
        //
	// uint64_t positive = 0;//positive offsets
        //
	// for(auto p : LZ77k){
        //
	// 	positive += p.first>0;
        //
	// }
        //
	// cout << "positive offsets: " << positive << endl;
	// cout << "negative offsets: " << LZ77k.size()-positive << endl;
        //
	// int bucket_size = 1;
        //
	// /*auto buckets = Vector<uint64_t>(bwt.size()/bucket_size + 1);
        //
	// for(auto p : LZ77k){
        //
	// 	buckets[(p.first<0?-p.first:p.first)/bucket_size]++;
        //
	// }
        //
	// for(int i=0;i<1000;++i){
        //
	// 	//cout << "[" << i*bucket_size << "," << (i+1)*bucket_size << ") : " << buckets[i] << endl;
	// 	cout << i << "\t" << buckets[i] << endl;
        //
	// }*/
	// Vector<uint64_t> abs_off;
	// for(auto x : LZ77k) abs_off.push_back(x.first<0?-x.first:x.first);
        //
	// Vector<uint64_t> Len;
	// for(auto x : LZ77k) Len.push_back(x.second);
        //
	// uint64_t sum_log = 0;
        //
	// for(auto x:abs_off) sum_log += bit_size(x)+1;
        //
	// auto G = alphabet.size()*8 + compute_gamma_bit_complexity(LZ77k);
	// auto D = alphabet.size()*8 + compute_delta_bit_complexity(LZ77k);
	// auto SU = alphabet.size()*8 + compute_succinct_bit_complexity(LZ77k);
	// auto EO = (entropy(abs_off)+1);
	// auto EL = entropy(Len);
	// auto low_bound = (EO+EL)*LZ77k.size();
        //
	// cout << "number of phrases = " << LZ77k.size() << endl<<endl;
        //
	// cout << "Entropy of the offsets: " << EO << endl;
	// cout << "Entropy of the lengths: " << EL << endl;
	// cout << "Lower bound for encoding the pairs: " << low_bound << " bits (" << low_bound/double(N) << " bit/char, " << low_bound/8 + 1 << " Bytes)" << endl << endl;
        //
	// cout << "Sum of logs of the offsets: " << sum_log << endl;
	// cout << "Average of logs of the offsets: " << double(sum_log)/double(LZ77k.size()) << endl<<endl;
        //
	// cout << "gamma complexity of the output: " << G/8+1 << " Bytes, " << double(G)/double(N) << " bit/symbol" << endl;
	// cout << "delta complexity of the output: " << D/8+1 << " Bytes, " << double(D)/double(N) << " bit/symbol" << endl;
	// cout << "Succinct complexity of the output: " << SU/8+1 << " Bytes, " << double(SU)/double(N) << " bit/symbol" << endl;

}

void run_pairs(const string& filePath){

	cout << "Computing the parse in format (occ,len)" << endl;

	wt_bwt bwt;
	set<char> alphabet;

	{
		ifstream in(filePath);
		auto F = get_frequencies(in);

		for(auto f : F) if(f.second>0) alphabet.insert(f.first);

		bwt = wt_bwt(F); //Huffman-encoded BWT
	}

	// prepend the alphabet to the text

	string prefix;

	for (auto rit = alphabet.rbegin(); rit != alphabet.rend(); rit++){

		auto c = *rit;
		prefix += c;
		bwt.extend(uint8_t(c));

	}

	// process the text

	Vector<pair<int64_t, uint64_t> > LZ77k;// the parse

	{
		ifstream in(filePath);

		char c;
		auto range = bwt.get_full_interval();// interval of current phrase
		pair<uint64_t,uint64_t> prev_range;
		uint64_t index = bwt.get_terminator_position();// co-lex position where current phrase should be if inserted
		string phrase; // current phrase

		if(in.is_open()) {

			if(in.good())
				in.get(c);//get character

			uint64_t read_char = 1;

			while(in.good()) {

				prev_range = range;
				range = bwt.LF(range, uint8_t(c));

				if(empty_interval(range)){//end of phrase

					output_phrase(bwt,LZ77k,prev_range,range,index,phrase);

				}else{

					index = bwt.LF(index,uint8_t(c));
					phrase += c;
					in.get(c);//get next character
					read_char++;

					if(read_char%100000 == 0){

						cout << "read " << read_char << " characters." << endl;

					}

				}

			}

		}

		//last phrase has not been output
		if(phrase.length()>0){

			prev_range = range;
			output_phrase(bwt,LZ77k,prev_range,range,index,phrase);

		}

		if(!in.eof() && in.fail())
			cout << "Error reading " << filePath << endl;

	}


	cout << "factorization: " << endl;
	for(auto p : LZ77k){

		auto d = 1+delta(p.first<0?-p.first:p.first);
		auto l = delta(p.second);
		cout << "off = " << p.first << " (" << d  << "  bits), len = " << p.second << " (" << l << " bits). Tot: " << (double(d+l)/double(p.second)) << " bit/symbol)" << endl;

	}

	auto N = bwt.size() -1;//file length


	uint64_t positive = 0;//positive offsets

	for(auto p : LZ77k){

		positive += p.first>0;

	}

	cout << "positive offsets: " << positive << endl;
	cout << "negative offsets: " << LZ77k.size()-positive << endl;

	int bucket_size = 1;

	/*auto buckets = Vector<uint64_t>(bwt.size()/bucket_size + 1);

	for(auto p : LZ77k){

		buckets[(p.first<0?-p.first:p.first)/bucket_size]++;

	}

	for(int i=0;i<1000;++i){

		//cout << "[" << i*bucket_size << "," << (i+1)*bucket_size << ") : " << buckets[i] << endl;
		cout << i << "\t" << buckets[i] << endl;

	}*/
	Vector<uint64_t> abs_off;
	for(auto x : LZ77k) abs_off.push_back(x.first<0?-x.first:x.first);

	Vector<uint64_t> Len;
	for(auto x : LZ77k) Len.push_back(x.second);

	uint64_t sum_log = 0;

	for(auto x:abs_off) sum_log += bit_size(x)+1;

	auto G = alphabet.size()*8 + compute_gamma_bit_complexity(LZ77k);
	auto D = alphabet.size()*8 + compute_delta_bit_complexity(LZ77k);
	auto SU = alphabet.size()*8 + compute_succinct_bit_complexity(LZ77k);
	auto EO = (entropy(abs_off)+1);
	auto EL = entropy(Len);
	auto low_bound = (EO+EL)*LZ77k.size();

	cout << "number of phrases = " << LZ77k.size() << endl<<endl;

	cout << "Entropy of the offsets: " << EO << endl;
	cout << "Entropy of the lengths: " << EL << endl;
	cout << "Lower bound for encoding the pairs: " << low_bound << " bits (" << low_bound/double(N) << " bit/char, " << low_bound/8 + 1 << " Bytes)" << endl << endl;

	cout << "Sum of logs of the offsets: " << sum_log << endl;
	cout << "Average of logs of the offsets: " << double(sum_log)/double(LZ77k.size()) << endl<<endl;

	cout << "gamma complexity of the output: " << G/8+1 << " Bytes, " << double(G)/double(N) << " bit/symbol" << endl;
	cout << "delta complexity of the output: " << D/8+1 << " Bytes, " << double(D)/double(N) << " bit/symbol" << endl;
	cout << "Succinct complexity of the output: " << SU/8+1 << " Bytes, " << double(SU)/double(N) << " bit/symbol" << endl;

}


/*
 * minimize the bit-size of each phrase, possibly cutting it
 */
void run_pairs_local(string filePath){

	cout << "Computing the locally-optimal parse in format (occ,len)" << endl;

	wt_bwt bwt;
	set<char> alphabet;

	{
		ifstream in(filePath);
		auto F = get_frequencies(in);

		for(auto f : F) if(f.second>0) alphabet.insert(f.first);

		bwt = wt_bwt(F); //Huffman-encoded BWT

	}

	// prepend the alphabet to the text

	cout << "Alphabet = ";

	for (auto rit = alphabet.rbegin(); rit != alphabet.rend(); rit++){

		auto c = *rit;
		cout << c;
		bwt.extend(uint8_t(c));

	}

	cout << endl;

	string text;

	//read text
	{
		std::ifstream inFile;
		inFile.open(filePath);
		std::stringstream strStream;
		strStream << inFile.rdbuf(); //read the file
		text = strStream.str();
	}

	cout << "text length = " << text.length() << endl;

	// process the text

	Vector<pair<int64_t, uint64_t> > LZ77k;// the parse
	uint64_t pos = 0;//current position on the text
	uint64_t k = 0;//current position inside the current phrase

	{

		auto range = bwt.get_full_interval();// interval of current phrase
		uint64_t index = bwt.get_terminator_position();// co-lex position where current phrase should be if inserted
		uint64_t phrase_len = 0; // current phrase length

		char c = text[pos + (k++)]; //read first character

		//length that optimizes the bitsize per character of the phrase
		double opt_bitsize = 10000; // no offset will ever take this number of bits, right?
		uint64_t opt_len = 0;
		int64_t opt_off = 0;

		int64_t off = 0;

		double perc = 0;
		double prev_perc = 0;

		//pos+k is the position of the character following the last character that has been read
		while(pos+k<=text.length()) {

			range = bwt.LF(range, uint8_t(c));

			if(empty_interval(range)){//no more previous matches: choose the length that minimizes the bit-size of the phrase

				LZ77k.push_back({opt_off,opt_len}); // push the optimum

				//extend BWT with the phrase
				for(uint64_t i = pos; i<pos+opt_len; ++i)
					bwt.extend(uint8_t(text[i]));

				k = 0;
				pos += opt_len;

				range = bwt.get_full_interval();
				index = bwt.get_terminator_position();
				phrase_len = 0;

				c = text[pos + (k++)];//read first character after phrase

				opt_bitsize = 10000;
				opt_len = 0;
				opt_off = 0;
				off = 0;

				perc = (100*double(pos))/text.length();
				cout << perc << "% done." << endl;
				if(perc >= prev_perc + 1){
					cout << perc << "% done." << endl;
					prev_perc = perc;
				}

			}else{

				index = bwt.LF(index,uint8_t(c));
				phrase_len++;

				off = get_off(bwt, range, index, phrase_len);
				auto bitsize = 1 + delta(off<0?-off:off) + delta(phrase_len);
				double bit_per_char = double(bitsize)/double(phrase_len);

				//found a new local optimum
				//if(bit_per_char < opt_bitsize and (phrase_len==1 or (off<0?-off:off)<128)){
				if(bit_per_char < opt_bitsize){
					opt_bitsize = bit_per_char;
					opt_len = phrase_len;
					opt_off = off;
				}

				//read next char
				if(pos+k<text.length())
					c = text[pos + k];

				k++;

			}

		}


		//last phrase has not been output
		if(phrase_len>0){

			assert(off>0);
			LZ77k.push_back({off,phrase_len});

		}

	}

	cout << "factorization: " << endl;
	for(auto p : LZ77k){

		auto d = 1+delta(p.first<0?-p.first:p.first);
		auto l = delta(p.second);
		cout << "off = " << p.first << " (" << d  << "  bits), len = " << p.second << " (" << l << " bits). Tot: " << (double(d+l)/double(p.second)) << " bit/symbol)" << endl;

	}

	auto N = bwt.size() -1;//file length


	uint64_t positive = 0;//positive offsets

	for(auto p : LZ77k){

		positive += p.first>0;

	}

	cout << "positive offsets: " << positive << endl;
	cout << "negative offsets: " << LZ77k.size()-positive << endl;

	int bucket_size = 1;

	/*auto buckets = Vector<uint64_t>(bwt.size()/bucket_size + 1);

	for(auto p : LZ77k){

		buckets[(p.first<0?-p.first:p.first)/bucket_size]++;

	}

	for(int i=0;i<1000;++i){

		//cout << "[" << i*bucket_size << "," << (i+1)*bucket_size << ") : " << buckets[i] << endl;
		cout << i << "\t" << buckets[i] << endl;

	}*/
	Vector<uint64_t> abs_off;
	for(auto x : LZ77k) abs_off.push_back(x.first<0?-x.first:x.first);

	Vector<uint64_t> Len;
	for(auto x : LZ77k) Len.push_back(x.second);

	uint64_t sum_log = 0;

	for(auto x:abs_off) sum_log += bit_size(x)+1;

	auto G = alphabet.size()*8 + compute_gamma_bit_complexity(LZ77k);
	auto D = alphabet.size()*8 + compute_delta_bit_complexity(LZ77k);
	auto EO = (entropy(abs_off)+1);
	auto EL = entropy(Len);
	auto low_bound = (EO+EL)*LZ77k.size();

	cout << "number of phrases = " << LZ77k.size() << endl<<endl;

	cout << "Entropy of the offsets: " << EO << endl;
	cout << "Entropy of the lengths: " << EL << endl;
	cout << "Lower bound for encoding the pairs: " << low_bound << " bits (" << low_bound/double(N) << " bit/char, " << low_bound/8 + 1 << " Bytes)" << endl << endl;

	cout << "Sum of logs of the offsets: " << sum_log << endl;
	cout << "Average of logs of the offsets: " << double(sum_log)/double(LZ77k.size()) << endl<<endl;

	cout << "gamma complexity of the output: " << G/8+1 << " Bytes, " << double(G)/double(N) << " bit/symbol" << endl;
	cout << "delta complexity of the output: " << D/8+1 << " Bytes, " << double(D)/double(N) << " bit/symbol" << endl;

}





void run_triples(string filePath){

	cout << "Computing the parse in format (occ,len,char)" << endl;

	wt_bwt bwt;
	set<char> alphabet;

	{
		ifstream in(filePath);
		auto F = get_frequencies(in);

		for(auto f : F) if(f.second>0) alphabet.insert(f.first);

		bwt = wt_bwt(F); //Huffman-encoded BWT
	}

	// the parse
	Vector<pair<int64_t, uint64_t> > LZ77k;
	Vector<char> trail_chars;

	{
		ifstream in(filePath);

		char c;
		auto range = bwt.get_full_interval();// interval of current phrase
		pair<uint64_t,uint64_t> prev_range;
		uint64_t index = bwt.get_terminator_position();// co-lex position where current phrase should be if inserted
		string phrase; // current phrase

		if(in.is_open()) {

			if(in.good())
				in.get(c);//get character

			uint64_t read_char = 1;

			while(in.good()) {

				prev_range = range;
				range = bwt.LF(range, uint8_t(c));

				phrase += c;

				if(empty_interval(range)){//end of phrase

					output_phrase(bwt,LZ77k,trail_chars,prev_range,range,index,phrase);

				}else{

					index = bwt.LF(index,uint8_t(c));

				}

				in.get(c);//get next character
				read_char++;

				if(read_char%100000 == 0){

					cout << "read " << read_char << " characters." << endl;

				}

			}

		}

		//last phrase has not been output
		if(phrase.length()>0){

			output_phrase(bwt,LZ77k,trail_chars,prev_range,range,index,phrase);

		}

		if(!in.eof() && in.fail())
			cout << "Error reading " << filePath << endl;

	}

	/*cout << "factorization: " << endl;
	for(auto p : LZ77k){

		cout << p.first << ", " << p.second << endl;

	}*/

	auto N = bwt.size() -1;//file length


	uint64_t positive = 0;//positive offsets

	for(auto p : LZ77k){

		positive += p.first>0;

	}

	cout << "positive offsets: " << positive << endl;
	cout << "negative offsets: " << LZ77k.size()-positive << endl;

	/*int bucket_size = 1;

	auto buckets = Vector<uint64_t>(bwt.size()/bucket_size + 1);

	for(auto p : LZ77k){

		buckets[(p.first<0?-p.first:p.first)/bucket_size]++;

	}

	for(int i=0;i<1000;++i){

		//cout << "[" << i*bucket_size << "," << (i+1)*bucket_size << ") : " << buckets[i] << endl;
		cout << i << "\t" << buckets[i] << endl;

	}*/

	uint64_t gamma_trail = 0;
	uint64_t delta_trail = 0;

	for(auto c : trail_chars){

		gamma_trail += gamma(uint64_t(uint8_t(c))+1);
		delta_trail += delta(uint64_t(uint8_t(c))+1);

	}

	auto G = gamma_trail+compute_gamma_bit_complexity(LZ77k);
	auto D = delta_trail+compute_delta_bit_complexity(LZ77k);
	auto H = compute_Huff_bit_complexity(LZ77k,trail_chars);

	Vector<uint64_t> abs_off;
	for(auto x : LZ77k) abs_off.push_back(x.first<0?-x.first:x.first);

	uint64_t sum_log = 0;

	for(auto x:abs_off) sum_log += bit_size(x)+1;

	cout << "number of phrases = " << LZ77k.size() << endl;
	cout << "Entropy of the offsets: " << (entropy(abs_off)+1) << endl;
	cout << "Sum of logs of the offsets: " << sum_log << endl;
	cout << "Average of logs of the offsets: " << double(sum_log)/double(LZ77k.size()) << endl;
	cout << "gamma complexity of the output: " << (G/8)+1 << " Bytes, " << double(G)/double(N) << " bit/symbol" << endl;
	cout << "delta complexity of the output: " << (D/8)+1 << " Bytes, " << double(D)/double(N) << " bit/symbol" << endl;
	cout << "Huffman complexity of the output: " << (H/8)+1 << " Bytes, " << double(H)/double(N) << " bit/symbol" << endl;


}



int main(int argc,char** argv){
    // for(uint64_t i = 0; i < 12; ++i) {
	// std::cout << "bits of gamma(" << i << ") = " << gamma(i+1) << std::endl;
    // }

	if(argc<2) help();

	bool triples = false;
	bool local = false;
	string filePath;

	int opt;
	while ((opt = getopt(argc, argv, "i:thl")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 't':
				triples = true;
			break;
			case 'l':
				local = true;
			break;
			case 'i':
				filePath = string(optarg);
			break;
			default:
				help();
			return -1;
		}
	}


	if(triples and local) help();
	if(filePath.length()==0) help();

	cout << "Input file " << filePath << endl;

	if(triples){

		run_triples(filePath);

	}else if(local){

		run_pairs_local(filePath);

	}else{

	    run_lz77(filePath);
		run_pairs_bitoptimal(filePath);
	//	run_pairs(filePath);

	}

	return 0;

}

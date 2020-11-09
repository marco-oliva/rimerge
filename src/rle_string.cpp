//
//  rled_string.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <rimerge/rle_string.hpp>


namespace rimerge
{

//------------------------------------------------------------------------------

void
RLEString::Metadata::init(std::string path, read_tag tag)
{
    reading = true;
    this->file_path = path;
    
    std::ifstream in(file_path, std::ios::binary);
    in.read(reinterpret_cast<char*>(&size), sizeof(size_type));
    in.read(reinterpret_cast<char*>(&runs), sizeof(size_type));
    
    for (size_type i = 0; i < alphabet_max_size(); i++) { in.read(reinterpret_cast<char*>(&size_per_char[i]), sizeof(size_type)); }
    for (size_type i = 0; i < alphabet_max_size(); i++) { in.read(reinterpret_cast<char*>(&runs_per_char[i]), sizeof(size_type)); }
    
    in.close();
}

void
RLEString::Metadata::init(std::string path, write_tag tag)
{
    writing = true;
    this->file_path = path;
}

void
RLEString::Metadata::close()
{
    if (closed or reading) return;
    
    std::ofstream out(file_path, std::ios::binary);
    out.write(reinterpret_cast<char*>(&size), sizeof(size_type));
    out.write(reinterpret_cast<char*>(&runs), sizeof(size_type));
    
    for (size_type i = 0; i < alphabet_max_size(); i++) { out.write(reinterpret_cast<char*>(&size_per_char[i]), sizeof(size_type)); }
    for (size_type i = 0; i < alphabet_max_size(); i++) { out.write(reinterpret_cast<char*>(&runs_per_char[i]), sizeof(size_type)); }
    
    out.close();
}

//------------------------------------------------------------------------------

void
RLEString::RLEncoder::write(byte_type c, size_type len, std::ofstream& out)
{
    packed_type tmp = 0;
    if (len <= MAX_LEN)
    {
        tmp = len << 8;
        tmp = tmp | (packed_type) c;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(packed_type));
    }
    else
    {
        tmp = MAX_LEN << 8;
        tmp = tmp | NEXT_RECORD;
        tmp = tmp | (packed_type) c;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(packed_type));
        
        write(c, len - MAX_LEN, out);
    }
}

void
RLEString::RLEncoder::close()
{
    if (closed) return; closed = true;

    write(curr_run_char, curr_run_length, out_stream);
    out_stream.close();
    
    // Metadata
    metadata.close();
}

void
RLEString::RLEncoder::append(byte_type c)
{
    metadata.size += 1; metadata.size_per_char[c] += 1;
    
    if (curr_run_char == '\0')
    {
        curr_run_length = 1; curr_run_char = c;
        metadata.runs = 1; metadata.runs_per_char[c] +=1;
        return;
    } // first char
    if (c == curr_run_char) { curr_run_length++; return; }
    
    metadata.runs += 1; metadata.runs_per_char[c] +=1;
    
    write(curr_run_char, curr_run_length, out_stream);
    curr_run_char = c; curr_run_length = 1;
}

//------------------------------------------------------------------------------

RLEString::RLEDecoder::RLEDecoder(std::string path) : in_stream(path)
{
    // read metadata
    metadata.init(path + ".meta", Metadata::read_tag());
}

RunType
RLEString::RLEDecoder::next()
{
    RunType out;
    bool next = true;
    
    while (next)
    {
        packed_type extracted;
        in_stream.read(reinterpret_cast<char*>(&extracted), sizeof(packed_type));
        next = extracted & RLEncoder::NEXT_RECORD;
        
        out.charachter = (byte_type) (extracted & RLEncoder::CHAR_MASK);
        out.length += (size_type) ((extracted & RLEncoder::LENGTH_MASK) >> RLEncoder::CHAR_BITS);
    }
    
    runs_served++; return out;
}

bool
RLEString::RLEDecoder::end()
{
    return runs_served >= metadata.runs;
}

//------------------------------------------------------------------------------

void
RLEString::RLEncoderMerger::merge()
{
    double start = readTimer();
    
    if (merged) return;
    merged = true;
    
    std::ofstream out_stream(out_path);
    
    // close all encoders
    for (size_type i = 0; i < encoders.size(); i++) { encoders[i].second.close(); }
    
    // Merge run length files
    for (size_type i = 0; i < encoders.size() - 1; i++)
    {
        // metadata
        metadata.size += encoders[i].second.metadata.size;
        metadata.runs += encoders[i].second.metadata.runs;
        for (size_type j = 0; j < alphabet_max_size(); j++)
        {
            metadata.size_per_char[j] += encoders[i].second.metadata.size_per_char[j];
            metadata.runs_per_char[j] += encoders[i].second.metadata.runs_per_char[j];
        }
        
        // read frist next
        std::ifstream next(encoders[i+1].first, std::ios::binary);
        packed_type next_start;
        next.read(reinterpret_cast<char*>(&next_start), sizeof(packed_type));
        
        // read last curr
        std::ifstream curr(encoders[i].first, std::ios::binary);
        curr.seekg(- (sizeof(packed_type)), std::ios_base::end);
        packed_type curr_end;
        curr.read(reinterpret_cast<char*>(&curr_end), sizeof(packed_type));
        
        // write whole file to output
        curr.seekg(0, std::ios_base::beg);
        out_stream << curr.rdbuf();
        
        // overwrite last element if necessary
        byte_type next_start_char = next_start & RLEncoder::CHAR_MASK;
        byte_type curr_end_char = curr_end & RLEncoder::CHAR_MASK;
        if (next_start_char == curr_end_char)
        {
            metadata.runs -= 1;
            metadata.runs_per_char[curr_end_char] -= 1;
            
            out_stream.seekp(out_stream.tellp() - std::ios::pos_type(sizeof(packed_type)));
            curr_end = curr_end | RLEncoder::NEXT_RECORD;
            out_stream.write(reinterpret_cast<char*>(&curr_end), sizeof(packed_type));
        }
        
        curr.close();
        next.close();
    }
    
    // last file
    metadata.size += encoders[encoders.size() - 1].second.metadata.size;
    metadata.runs += encoders[encoders.size() - 1].second.metadata.runs;
    for (size_type i = 0; i < alphabet_max_size(); i++)
    {
        metadata.size_per_char[i] += encoders[encoders.size() - 1].second.metadata.size_per_char[i];
        metadata.runs_per_char[i] += encoders[encoders.size() - 1].second.metadata.runs_per_char[i];
    }
    std::ifstream in(encoders[encoders.size() - 1].first, std::ios::binary);
    out_stream << in.rdbuf();
    
    if(Verbosity::level >= Verbosity::BASIC)
    {
        double seconds = readTimer() - start;
        spdlog::info("Completed bwt merge in {} seconds", seconds);
    }
}

//------------------------------------------------------------------------------

byte_type RLEString::Accessor::operator[](size_type i)
{
    // Check if value in cache
    for (size_type p = 0; p < cache_size; p++)
    {
        if (cache[p].first == i)
        {
            return cache[p].second;
        }
    }

    // Value not in cache
    byte_type val = string[i];
    cache[cache_it] = std::make_pair(i, val);
    cache_it++; cache_it = cache_it % cache_size;
    return val;
}

//------------------------------------------------------------------------------


RLEString::RLEString(string_type& data)
{
    size_type size = data.size();
    this->n = 1;
    this->R = 0;

    if (size <= 1)
    {
        this->n = 0;
        return;
    }

    byte_type* start = data.data();
    byte_type* iterator = data.data();
    byte_type* end = start + size;

    auto runs_per_letter_bv = std::vector<std::vector<bool>>(256);

    //runs in main bitvector
    std::vector<bool> runs_bv;
    runs_bv.reserve(size);
    string_type run_heads_s;

    byte_type last_c = *iterator;
    byte_type curr_c;

    while (++iterator != end)
    {
        curr_c = *iterator;
        if (last_c == 0) last_c = DATA_TERMINATOR;
        if (curr_c != last_c)
        {
            run_heads_s.push_back(last_c);
            runs_per_letter_bv[last_c].push_back(true);

            last_c = curr_c;

            //push back a bit set only at the end of a block
            runs_bv.push_back((R % B) == (B - 1));

            this->R++;
        }
        else
        {
            runs_bv.push_back(false);
            runs_per_letter_bv[last_c].push_back(false);
        }
        this->n++;
    }

    run_heads_s.push_back(last_c);
    runs_per_letter_bv[last_c].push_back(true);
    runs_bv.push_back(false);
    this->R++;

    runs = bitvector(runs_bv);

    //a fast direct array: char -> bitvector.
    for(size_type i=0;i<256;++i)
        runs_per_letter[i] = bitvector(runs_per_letter_bv[i]);

    run_heads = huff_string(run_heads_s);

    assert(run_heads.size()==R);
}

byte_type
RLEString::operator[](size_type i) const
{
    return run_heads[run_of(i).first];
}

size_type
RLEString::select(size_type i, byte_type c) const
{
    assert(i<runs_per_letter[c].size());

    //i-th c is inside j-th c-run (j starts from 0)
    assert(i<runs_per_letter[c].size());
    size_type j = runs_per_letter[c].rank(i);

    //starting position of i-th c inside its run
    assert(j==0 || i >= runs_per_letter[c].select(j-1) + 1);
    size_type before = (j==0 ? i : i - (runs_per_letter[c].select(j-1) + 1));

    //position in run_heads
    size_type r = run_heads.select(j,c);

    //k = number of bits before position of interest in the main string
    //here, k is initialized looking at the sampled runs
    assert(r/B==0 || r/B-1<runs.number_of_1());
    size_type k = (r/B==0?0 : runs.select(r/B-1)+1);

    //now add remaining run lengths to k
    for( size_type t = (r/B)*B; t<r; ++t )
    {
        k += run_at(t);
    }

    return k + before;
}

size_type
RLEString::rank(size_type i, byte_type c) const
{
    assert(i<=n);

    //letter does not exist in the text
    if(runs_per_letter[c].size()==0) return 0;

    if(i==n) return runs_per_letter[c].size();

    size_type last_block = runs.rank(i);
    size_type current_run = last_block*B;

    //current position in the string: the starts of a block
    size_type pos = 0;
    if(last_block>0)
        pos = runs.select(last_block-1)+1;

    assert(pos <= i);

    size_type dist = i-pos;

    //otherwise, scan at most B runs
    while(pos < i)
    {
        pos += run_at(current_run);
        current_run++;

        if(pos<=i) dist = i-pos;
    }

    if(pos>i) current_run--;

    //position i is inside run current_run
    assert(current_run<R);

    //number of c runs before the current run
    size_type rk = run_heads.rank(current_run,c);

    //number of c before i in the current run
    size_type tail = (run_heads[current_run]==c)*dist;

    //in this case, either there are no c before position i
    //or the current run is the starts of the kind ccc...cc
    if(rk==0) return tail;

    return runs_per_letter[c].select(rk-1)+1+tail;
}

size_type
RLEString::run_of_position(size_type i) const
{
    assert(i<n);

    size_type last_block = runs.rank(i);
    size_type current_run = last_block*B;

    //current position in the string: the starts of a block
    size_type pos = 0;
    if(last_block>0)
        pos = runs.select(last_block-1)+1;

    assert(pos <= i);

    size_type dist = i-pos;

    //otherwise, scan at most B runs
    while(pos < i){

        pos += run_at(current_run);
        current_run++;

        if(pos<=i) dist = i-pos;

    }

    if(pos>i) current_run--;

    //position i is inside run current_run
    assert(current_run<R);

    return current_run;

}

std::vector<range_type>
RLEString::break_range(range_type rn, byte_type c) const
{
    auto l = rn.first;
    auto r = rn.second;

    assert(l<=r);
    assert(r<size());

    assert(operator[](l)==c);
    assert(operator[](r)==c);

    //retrieve runs that contain positions l and r
    auto run_l = run_of(l);
    auto run_r = run_of(r);

    //in this case rn contains only character c: do not break
    if(run_l.first==run_r.first) return {rn};

    std::vector<range_type> result;

    //starts range: from l to the end of the run containing position l
    result.push_back({l,run_l.second});

    //rank of c's of interest in run_heads
    size_type rank_l = run_heads.rank(run_l.first,c);
    size_type rank_r = run_heads.rank(run_r.first,c);

    //now retrieve run bounds of all c-runs of interest
    for(size_type j = rank_l+1;j<rank_r;++j)
        result.push_back(run_range(run_heads.select(j,c)));

    //now ends (possibly incomplete) run

    auto range = run_range(run_heads.select(rank_r,c));
    result.push_back({range.first,r});

    return result;
}

size_type
RLEString::size() const { return n; }

size_type
RLEString::number_of_runs() const { return this->R; }

range_type
RLEString::run_range(size_type j) const
{
    assert(j<run_heads.size());

    size_type this_block = j/B;
    size_type current_run = this_block*B;
    size_type pos = (this_block==0?0:runs.select(this_block-1)+1);

    while(current_run < j){

        pos += run_at(current_run);
        current_run++;

    }

    assert(current_run == j);

    return {pos,pos+run_at(j)-1};
}

size_type
RLEString::run_at(size_type i) const
{
    assert(i<R);
    byte_type c = run_heads[i];

    return runs_per_letter[c].gap_at(run_heads.rank(i,c));
}

void
RLEString::load(std::string path)
{
    RLEDecoder decoder(path);
    
    // Run Heads
    string_type run_heads_b;
    run_heads_b.reserve(decoder.metadata.runs);
    
    // Bitvector
    sdsl::sd_vector_builder runs_bv_builder(decoder.metadata.size, decoder.metadata.runs);
    size_type runs_bv_it = 0;
    
    // Runs per letter
    std::array<sdsl::sd_vector_builder, alphabet_max_size()> runs_per_letter_bv_builders;
    std::array<size_type, alphabet_max_size()> runs_per_letter_bvs_it = {0};
    for (size_type i = 0; i < alphabet_max_size(); i++)
    {
        runs_per_letter_bv_builders[i] = sdsl::sd_vector_builder(decoder.metadata.size_per_char[i], decoder.metadata.runs_per_char[i]);
    }
    
    while (not decoder.end())
    {
        RunType run = decoder.next();
        byte_type c = RunTraits::charachter(run);
        size_type len = RunTraits::length(run);
        
        if (len != 0)
        {
            run_heads_b.push_back(c);
            runs_bv_builder.set(runs_bv_it + len - 1); runs_bv_it += len;
            runs_per_letter_bv_builders[c].set(runs_per_letter_bvs_it[c] + len -1);  runs_per_letter_bvs_it[c] += len;
        }
    }
//    if (runs_bv.size() != 0)
//        runs_bv[runs_bv.size() - 1] = false;
    
    // Build static structures
    this->runs = bitvector(runs_bv_builder);
    this->n = decoder.metadata.size;
    
    this->R = run_heads_b.size();
    this->run_heads.init(run_heads_b);

    for(size_type i = 0; i < alphabet_max_size(); i++)
    {
        this->runs_per_letter[i] = bitvector(runs_per_letter_bv_builders[i]);
    }
}

string_type
RLEString::to_string() const
{
    string_type s;

    for(size_type i=0;i<size();++i)
        s.push_back(operator[](i));

    return s;
}

size_type
RLEString::closest_run_break(range_type rn, byte_type c) const
{
    /*
     * case 1: range begins with a c-run: return ends position of the run
     */
    if(operator[](rn.first)==c)
    {
        size_type i = run_of_position(rn.first);
        size_type j = run_range(i).second;

        //j must be inside rn, i.e. rn must not contain only c
        //j must not be ends position of rn: this would imply
        //that rn contain only c
        assert(j<rn.second);
        return j;
    }
    else
    {
        //case 2: starts c-run starts in the middle of the range

        //rank i of starts c in the range
        size_type i = rank(rn.first,c);

        assert(i<rank(size(),c));

        //map from rank space to string position:
        //i now is the starts position inside the range that contains c
        i = select(i,c);

        assert(operator[](i)==c);
        assert(i<=rn.second);

        return i;
    }
}

std::pair<size_type,size_type>
RLEString::run_of(size_type i) const
{

    size_type last_block = runs.rank(i);
    size_type current_run = last_block*B;

    //current position in the string: the starts of a block
    size_type pos = 0;
    if(last_block>0)
        pos = runs.select(last_block-1)+1;

    assert(pos <= i);

    while(pos < i)
    {
        pos += run_at(current_run);
        current_run++;
    }

    assert(pos >= i);

    if(pos>i)
    {
        current_run--;
    }
    else
    {//pos==i
        pos += run_at(current_run);
    }

    assert(pos>0);
    assert(current_run<R);

    return {current_run,pos-1};
}

//------------------------------------------------------------------------------

}


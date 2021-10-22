//
//  rle_string.hpp
//
//  Run length encoded string heavly inspired by Nicola Prezza's implementation.
//  https://github.com/nicolaprezza/r-index
//

#ifndef rle_string_hpp
#define rle_string_hpp

#include <rimerge/utils.hpp>
#include <rimerge/bitvector.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

// Huffman string
class huff_string
{
public:
    
    huff_string() = default;
    
    huff_string(string_type& s)
    {
        s.push_back(IMPL_TERMINATOR);
        construct_im(wt, (const char*) s.data(), 1);
        assert(wt.size() == s.size() - 1);
    }
    
    void
    init(string_type& s)
    {
        s.push_back(IMPL_TERMINATOR);
        construct_im(wt, (const char*) s.data(), 1);
        assert(wt.size() == s.size() - 1);
    }
    
    byte_type
    operator[](size_type i) const
    {
        assert(i < wt.size());
        return wt[i];
    }
    
    size_type
    size() const
    {
        return wt.size();
    }
    
    size_type
    rank(size_type i, byte_type c) const
    {
        assert(i <= wt.size());
        return wt.rank(i, c);
    }
    
    // position of i-th character c. i starts from 0!
    size_type
    select(size_type i, byte_type c) const
    {
        return wt.select(i + 1, c);
    }
    
    // serialize the structure to the ostream
    size_type
    serialize(std::ostream& out) const
    {
        return wt.serialize(out);
    }
    
    // load the structure from the istream
    void
    load(std::istream& in)
    {
        wt.load(in);
    }

private:
    
    sdsl::wt_huff<> wt;
};

//------------------------------------------------------------------------------

struct RunType
{
    size_type offset = invalid_value(); // string position of start char
    size_type length = 0;
    byte_type charachter = '\0';
};

class RunTraits
{
public:
    
    static byte_type charachter(const RunType& run) { return run.charachter; }
    static size_type length(const RunType& run) { return run.length; }
    static size_type start(const RunType& run) { return run.offset; }
    static size_type end(const RunType& run) { return run.offset + run.length - 1; }
    static bool start(const RunType& run, size_type string_position) { return string_position == run.offset;}
    static bool end(const RunType& run, size_type string_position) { return string_position == run.offset + run.length;}
    static bool equal(const RunType& left, const RunType& right) { return left.offset == right.offset; }
};

//------------------------------------------------------------------------------

class RLEString
{

public:

//------------------------------------------------------------------------------

    class Metadata
    {
    private:
        
        bool closed = false;
        bool reading = false, writing = false;
        
        std::string file_path;
    
    public:
        
        size_type size = 0;
        size_type runs = 0;
        
        std::array<size_type, alphabet_max_size()> size_per_char = {0};
        std::array<size_type, alphabet_max_size()> runs_per_char = {0};
        
        struct read_tag  {};
        struct write_tag {};
        
        Metadata(std::string path, read_tag tag) { init(path, tag); }
        Metadata(std::string path, write_tag tag) { init(path, tag); }
        Metadata() {}
        ~Metadata() {  if (not closed) close(); };
        
        void init(std::string path, read_tag tag);
        void init(std::string path, write_tag tag);
    
        void close();
    };

//------------------------------------------------------------------------------
    
    class RLEncoder
    {
    public:
        
        using packed_type = uint32_t;
    
        constexpr static packed_type LENGTH_BITS = 3 * 8;
        constexpr static packed_type LENGTH_MASK = 0x7FFFFF00;
        constexpr static packed_type NEXT_RECORD = 0x80000000;
        constexpr static packed_type MAX_LEN     = 0x7FFFFF;
    
        constexpr static packed_type CHAR_BITS   = 1 * 8;
        constexpr static packed_type CHAR_MASK   = 0xFF;
    
    public:
        std::string out_path;
        std::ofstream out_stream;
        bool closed = false;
        
        size_type curr_run_length = 0;
        byte_type curr_run_char   = 0;
        
        void append(byte_type c);
        
        void write(byte_type c, size_type len, std::ofstream& out);
    
    public:
    
        // Metadata
        Metadata metadata;
        
        RLEncoder(std::string path) { init(path); }
        RLEncoder() {}
        ~RLEncoder() { if (not closed) close(); }
        
        void init(std::string& path)
        {
            out_path = path;
            out_stream.open(path, std::ios::binary);
            if (not out_stream.is_open()) { spdlog::error("Can't open {}", path); std::exit(EXIT_FAILURE); }
            metadata.init(path + ".meta", Metadata::write_tag());
        }
        
        void operator()(byte_type c) { assert(not closed); /*if (c == '\0') { c = '$'; }*/ append(c); }
        void operator()(byte_type c, size_type len) { while (len > 0) { this->operator()(c); len -= 1; } }
        
        void close();
    };
    
    class RLEDecoder
    {
    
    private:
        
        using packed_type = RLEncoder::packed_type;
        std::ifstream in_stream;
        
        size_type runs_served = 0;
        
    public:
    
        // Metadata
        Metadata metadata;
        
        RLEDecoder(std::string path);
        
        RunType next();
        bool end();
        
    };
    
    class RLEncoderMerger
    {
    private:
    
        std::array<size_type, alphabet_max_size()> runs_per_char = {0};
        size_type size = 0;
        size_type runs = 0;
    
        using packed_type = RLEncoder::packed_type;
        
        std::string out_path;
        
        std::vector<std::pair<std::string, RLEncoder>> encoders;
        bool merged = false;
        
        void merge();
        
    public:
        
        Metadata metadata;
    
        RLEncoderMerger(std::string& path, size_type n) : out_path(path), encoders(n), metadata(path + ".meta", Metadata::write_tag())
        {
            for (auto& encoder : encoders)
            {
                encoder.first = TempFile::getName("rle");
                encoder.second.init(encoder.first);
            }
        }
        
        RLEncoder& get_encoder(size_type i)
        {
            assert(i < encoders.size());
            return encoders[i].second;
        }
        
        void close() { if (not merged) { merge(); }  metadata.close(); }
        
        ~RLEncoderMerger()
        {
            close();
            for (auto& e: encoders)
            {
                std::string r_meta = e.first + ".meta";
                if (not std::remove(r_meta.c_str())) spdlog::debug("can't remove: {} [{}]", r_meta, strerror(errno));
            }
        }
    
    };

//------------------------------------------------------------------------------
    
    class Accessor
    {
    private:
        
        static constexpr size_type cache_size = 8;
        size_type cache_it;
        std::array<std::pair<size_type, byte_type>, cache_size> cache;
        
        const RLEString& string;
        
    public:
        
        Accessor(const RLEString& s) :string(s), cache_it(0)
        {
            for (size_type i = 0; i < cache_size; i++)
            {
                cache[i] = std::make_pair(invalid_value(), '\0');
            }
        }
        
        byte_type operator[](size_type i);
    };

//------------------------------------------------------------------------------
    
    class RunIterator
    {
    
    private:
    
        const RLEString& string;
        RunType current_run;
        size_type r_pos;
        
        void read()
        {
            std::pair<size_type, size_type> run = string.run_range(r_pos);
            current_run.offset = run.first;
            current_run.length = run.second - run.first + 1;
            current_run.charachter = string.run_heads[r_pos];
        }
    
    public:
        
        RunIterator(const RLEString& rle_string) :string(rle_string), r_pos(0) { read(); }
        RunIterator(const RLEString& rle_string, size_type i) :string(rle_string)
        {
            r_pos = string.run_of_position(i);
            read();
            if (current_run.offset < i) current_run.offset = i;
        }
    
        void operator++() { r_pos += 1; read(); }
        const RunType& operator*() const { return current_run; }
        
        bool end() { return r_pos >= string.number_of_runs(); }
    };
    
//------------------------------------------------------------------------------
    
    class RunCache
    {
    
    private:
    
        static constexpr size_type cache_size = 2;
        size_type cache_it;
        std::array<RunType, cache_size> cache;
        const RLEString& string;
    
        void read(size_type i, RunType& out)
        {
            size_type r_pos = string.run_of_position(i);
            std::pair<size_type, size_type> run = string.run_range(r_pos);
            out.offset = run.first;
            out.length = run.second - run.first + 1;
            out.charachter = string.run_heads[r_pos];
        }
        
        
    public:
    
        RunCache(const RLEString& s) :string(s), cache_it(0) {}
        
        byte_type operator[](size_type pos)
        {
    
            for (size_type i = 0; i < cache_size; i++)
            {
                if ( pos >= RunTraits::start(cache[i]) and pos <= RunTraits::end(cache[i]))
                {
                    return RunTraits::charachter(cache[i]);
                }
            }
            
            // Get next run
            read(pos, cache[cache_it % cache_size]);
            char out = RunTraits::charachter(cache[cache_it % cache_size]);
            cache_it++;
            
            return out;
        }
        
    };
    
//------------------------------------------------------------------------------

    RLEString() = default;
    RLEString(string_type& data);
    
    ~RLEString()
    {
        if(Verbosity::level >= Verbosity::EXTENDED)
        {
            if (this->R != 0)
            spdlog::info("rimerge::RLEString::~RLEString: size {}GB",
                         inGigabytes((this->R) * (2 + std::log(this->n / this->R)) * 2));
        }
    }
    
    byte_type operator[](size_type i) const;
    
    // position of i-th character c. i starts from 0!
    size_type select(size_type i, byte_type c) const;
    
    // number of c before position i
    size_type rank(size_type i, byte_type c) const;
    
    // text position i is inside this run
    size_type run_of_position(size_type i) const;
    
    //break range: given a range <l',r'> on the string and a character c, this function
    //breaks <l',r'> in maximal sub-ranges containing character c.
    //for simplicity and efficiency, we assume that characters at range extremities are both 'c'
    //thanks to the encoding (run-length), this function is quite efficient: O(|result|) ranks and selects
    std::vector<range_type> break_range(range_type rn, byte_type c) const;
    
    size_type size() const;
    
    // inclusive range of j-th run in the string
    std::pair<size_type,size_type> run_range(size_type j) const;
    
    // length of i-th run
    size_type run_at(size_type i) const;
    
    size_type number_of_runs() const;
    
    // load the structure from the istream
    void load(std::string path);
    
    string_type to_string() const;
    
    /*
     * input: inclusive range rn, character c
     *
     * return the position j that is closest to rn.starts,
     * such that character in position j is c and that is
     * adjacent to a position j' inside rn that contains a
     * character != c
     *
     * rn must contain c and at least another character d!=c
     *
     */
    size_type closest_run_break(range_type rn, byte_type c) const;
    
    friend class RIndexRLEN;

private:
    
    // <j=run of position i, ends position of j-th run>
    std::pair<size_type,size_type> run_of(size_type i) const;
    
    // block size: bitvector 'runs' has R/B bits set (R being number of runs)
    std::size_t B = 1;
    
    bitvector runs;
    
    // for each letter, its runs stored contiguously
    std::array<bitvector, 256> runs_per_letter;
    
    // store run heads in a compressed string supporting access/rank
    huff_string run_heads;
    
    // text length and number of runs
    std::size_t n = 0;
    std::size_t R = 0;
    
    friend class RIndexRLE;
    
};
}

#endif //rle_string_hpp

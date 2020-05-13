/*
  Copyright (c) 2018 Jouni Siren

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#ifndef support_hpp
#define support_hpp

#include <rimerge/utils.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

struct MergeParameters
{
    constexpr static size_type POS_BUFFER_SIZE = 64; // Megabytes.
    constexpr static size_type THREAD_BUFFER_SIZE = 256; // Megabytes.
    constexpr static size_type MERGE_BUFFERS = 6;
    constexpr static size_type CHUNK_SIZE = 1; // Sequences per thread.
    constexpr static size_type MERGE_JOBS = 4;
    
    constexpr static size_type MAX_BUFFER_SIZE = 16384; // Megabytes.
    constexpr static size_type MAX_MERGE_BUFFERS = 16;
    constexpr static size_type MAX_MERGE_JOBS = 16;
    
    MergeParameters();
    
    void setPosBufferSize(size_type megabytes);
    void setThreadBufferSize(size_type megabytes);
    void setMergeBuffers(size_type n);
    void setChunkSize(size_type n);
    void setMergeJobs(size_type n);
    
    // These return the sizes in positions/bytes.
    size_type posBufferPositions() const { return (this->pos_buffer_size * MEGABYTE) / sizeof(rank_type); }
    size_type threadBufferBytes() const { return this->thread_buffer_size * MEGABYTE; }
    
    size_type pos_buffer_size, thread_buffer_size;
    size_type merge_buffers;
    size_type chunk_size;
    size_type merge_jobs;
    
    std::string out_prefix = "";
};


}


#endif //support_hpp

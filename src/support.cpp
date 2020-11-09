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

#include <rimerge/support.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

MergeParameters::MergeParameters() :
pos_buffer_size(POS_BUFFER_SIZE), thread_buffer_size(THREAD_BUFFER_SIZE),
merge_buffers(MERGE_BUFFERS), chunk_size(CHUNK_SIZE), merge_jobs(MERGE_JOBS),
iterations(ITERATIONS)
{
}

void
MergeParameters::setPosBufferSize(size_type megabytes)
{
    this->pos_buffer_size = Range::bound(megabytes, 1, MAX_BUFFER_SIZE);
}

void
MergeParameters::setThreadBufferSize(size_type megabytes)
{
    this->thread_buffer_size = Range::bound(megabytes, 1, MAX_BUFFER_SIZE);
}

void
MergeParameters::setMergeBuffers(size_type n)
{
    this->merge_buffers = Range::bound(n, 1, MAX_MERGE_BUFFERS);
}

void
MergeParameters::setChunkSize(size_type n)
{
    this->chunk_size = std::max(n, static_cast<size_type>(1));
}

void
MergeParameters::setMergeJobs(size_type n)
{
    this->merge_jobs = Range::bound(n, 1, MAX_MERGE_JOBS);
}

void
MergeParameters::setIterations(size_type n)
{
    this->iterations = Range::bound(n, 1, MAX_ITERATIONS);
}

}
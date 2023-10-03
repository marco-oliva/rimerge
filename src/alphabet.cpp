//
//  alphabet.cpp
//
// Copyright (c) Boucher Lab. All rights reserved.
// Licensed under the GNU license. See LICENSE file in the repository root for full license information.


#include <rimerge/alphabet.hpp>

namespace rimerge
{

//------------------------------------------------------------------------------

void
Alphabet::update(byte_type i)
{
    used_chars[i] = true;
    if (i == STRING_TERMINATOR) { used_terminator = STRING_TERMINATOR; }
    
    if (initialized)
    {
        spdlog::warn("rimerge::Aplhabet::update(): Updating already initialized alphabet");
        init();
    }
}

void
Alphabet::init()
{
    alphabet = bitvector(used_chars);
    initialized = true;
}

byte_type
Alphabet::previous(byte_type c) const
{
    size_type rank = alphabet.rank(c);
    
    if (rank == 0)
        return used_terminator;
    
    return static_cast<byte_type> (alphabet.select(rank - 1));
}

byte_type
Alphabet::following(byte_type c) const
{
    size_type rank = alphabet.rank(c);
    
    if ((rank + 1) == sigma())
        return used_terminator;
    
    size_type i = alphabet.select(rank + 1);
    
    return static_cast<byte_type>(i);
}

size_type
Alphabet::sigma() const
{
    return alphabet.number_of_1();
}

string_type
Alphabet::to_string() const
{
    string_type out;
    for (size_type i = 0; i < used_chars.size(); i++)
    {
        if (used_chars[i])
        {
            out.push_back((byte_type) i);
        }
    }
    
    return out;
}

//------------------------------------------------------------------------------

}
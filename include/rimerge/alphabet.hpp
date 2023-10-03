//
//  alphabet.hpp
//
// Copyright (c) Boucher Lab. All rights reserved.
// Licensed under the GNU license. See LICENSE file in the repository root for full license information.

#ifndef alphabet_hpp
#define alphabet_hpp

#include <rimerge/utils.hpp>
#include <rimerge/bitvector.hpp>

namespace rimerge
{

class Alphabet
{
public:
    
    Alphabet() = default;
    
    // Properties
    size_type sigma() const;
    
    
    void update(byte_type i);
    void init();
    byte_type previous(byte_type c) const;
    byte_type following(byte_type c) const;
    string_type to_string() const;
    
private:
    
    std::vector<bool> used_chars = std::vector<bool>(256,false);
    bool initialized = false;
    bitvector alphabet;
    byte_type used_terminator = DATA_TERMINATOR;
    
};

}

#endif //alphabet_hpp

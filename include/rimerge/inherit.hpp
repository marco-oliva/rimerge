//
//  inherit.hpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#ifndef inherit_hpp
#define inherit_hpp

namespace rimerge
{

/// Class that handles CRTP static inheritance (see J. Boccara's blog series)
template <typename Derived>
class inherit
{
public:
    
    /// Return the underlying class
    Derived& self()
    {
        return static_cast<Derived&>(*this);
    }
    
    /// Return the underlying class (const version)
    Derived const& self() const
    {
        return static_cast<Derived const&>(*this);
    }
};

}

#endif //inherit_hpp

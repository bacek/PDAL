//  (C) Copyright John Maddock 2005.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_TR1_MEMORY_HPP_INCLUDED
#  define BOOST_TR1_MEMORY_HPP_INCLUDED
#  include <boost/tr1/detail/config.hpp>
#  include <boost/detail/workaround.hpp>
#  include <memory>

#ifndef BOOST_HAS_TR1_SHARED_PTR

//
// This header can get included by boost/shared_ptr.hpp which leads
// to cyclic dependencies, the workaround is to forward declare all 
// the boost components, and then include the actual headers afterwards.
// This is fragile, but seems to work, and doesn't require modification
// of boost/shared_ptr.hpp.
//
namespace pdalboost {} namespace boost = pdalboost; namespace pdalboost{

class bad_weak_ptr;
template<class T> class weak_ptr;
template<class T> class shared_ptr;
template<class T> void swap(weak_ptr<T> & a, weak_ptr<T> & b);
template<class T> void swap(shared_ptr<T> & a, shared_ptr<T> & b);
template<class T, class U> shared_ptr<T> static_pointer_cast(shared_ptr<U> const & r);
template<class T, class U> shared_ptr<T> dynamic_pointer_cast(shared_ptr<U> const & r);
template<class T, class U> shared_ptr<T> const_pointer_cast(shared_ptr<U> const & r);
template<class D, class T> D * get_deleter(shared_ptr<T> const & p);
template<class T> class enable_shared_from_this;

namespace detail{
class shared_count;
class weak_count;
}

}

namespace std{ namespace tr1{

   using ::pdalboost::bad_weak_ptr;
   using ::pdalboost::shared_ptr;
#if !BOOST_WORKAROUND(__BORLANDC__, < 0x0582)
   using ::pdalboost::swap;
#endif
   using ::pdalboost::static_pointer_cast;
   using ::pdalboost::dynamic_pointer_cast;
   using ::pdalboost::const_pointer_cast;
   using ::pdalboost::get_deleter;
   using ::pdalboost::weak_ptr;
   using ::pdalboost::enable_shared_from_this;

} }
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#else

#  if defined(BOOST_HAS_INCLUDE_NEXT) && !defined(BOOST_TR1_DISABLE_INCLUDE_NEXT)
#     include_next BOOST_TR1_HEADER(memory)
#  else
#     include <boost/tr1/detail/config_all.hpp>
#     include BOOST_TR1_STD_HEADER(BOOST_TR1_PATH(memory))
#  endif

#endif

#endif


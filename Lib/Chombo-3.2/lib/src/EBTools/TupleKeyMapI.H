#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TUPLEKEYMAPI_H_
#define _TUPLEKEYMAPI_H_

#include <map>
#include <iostream>
#include "NamespaceHeader.H"

#define Tuple2 std::pair<T1,     T2>
#define Tuple3 std::pair<Tuple2, T3>
#define Tuple4 std::pair<Tuple3, T4>
#define TupleTypenames typename T1, typename T2, typename T3, typename T4
#define TupleArgTypes T1,T2,T3,T4
#define TupleArgValues arg1,arg2,arg3,arg4
#define TupleArgDecls T1 arg1, T2 arg2, T3 arg3, T4 arg4


/**
  Map whose key is a 4-tuple.  The 4 is hardcoded but can be pretty easily changed
  by making obvious modifications to the #define's above and going up to Tuple<n>
  in the functions below.
  Tuple types (T1,...) must implement both operator==() and operator<().
*/
template<TupleTypenames, typename ValueT> class TupleKeyMap
{
  typedef std::map< Tuple4, ValueT >         RepType;

public:
  typedef typename RepType::const_iterator   ConstIteratorType;

  void    insert( TupleArgDecls, ValueT );

  bool    containsKey( TupleArgDecls, ConstIteratorType* i=0  ) const;

  ValueT  fetch( TupleArgDecls ) const;

  ValueT  fetch( ConstIteratorType ) const;

  void    report() const;

  void    clear();

  unsigned size()
  {
    return m_rep.size();
  }

private:
  RepType m_rep;
};


/**
  Returns the value associated with the tuple key.
*/
template<TupleTypenames, typename ValueT> ValueT
TupleKeyMap<TupleArgTypes,ValueT>::fetch( TupleArgDecls ) const
{
  CH_assert( this->containsKey( TupleArgValues ) );

  Tuple2 tuple2(arg1,arg2);
  Tuple3 tuple3(tuple2,arg3);
  Tuple4 tuple4(tuple3,arg4);

  Tuple4 const & item( tuple4 );
  typename RepType::const_iterator iter = m_rep.find( item );
  return iter->second;
}

template<TupleTypenames, typename ValueT> ValueT
TupleKeyMap<TupleArgTypes,ValueT>::fetch( ConstIteratorType iter ) const
{
  return iter->second;
}


template<TupleTypenames, typename ValueT> bool
TupleKeyMap<TupleArgTypes,ValueT>::containsKey( TupleArgDecls,
                                                ConstIteratorType* iter) const
{
  Tuple2 tuple2(arg1,arg2);
  Tuple3 tuple3(tuple2,arg3);
  Tuple4 tuple4(tuple3,arg4);

  Tuple4 const & item( tuple4 );
  if ( iter )
  {
    *iter = m_rep.find( item );
    return !((*iter) == m_rep.end());
  } else
  {
    ConstIteratorType it = m_rep.find( item );
    return !(it == m_rep.end());
  }
}

template<TupleTypenames, typename ValueT> void
TupleKeyMap<TupleArgTypes,ValueT>::insert( TupleArgDecls, ValueT value )
{
  // Unlike std::map, which lets you insert a duplicate key (and thereupon
  // silently doesn't insert it), we consider that an error.  This is because
  // checking that the key is already in the map is costing us.
  CH_assert( ! this->containsKey( TupleArgValues ) );

  Tuple2 tuple2(arg1,arg2);
  Tuple3 tuple3(tuple2,arg3);
  Tuple4 tuple4(tuple3,arg4);

  Tuple4 const & item( tuple4 );
  m_rep[item] = value;
}


template<TupleTypenames, typename ValueT>
void
TupleKeyMap<TupleArgTypes,ValueT>::report() const
{
#ifndef NDEBUG
  std::cout << "TupleKeyMap::report(): unique arg combinations: "
            << m_rep.size() << std::endl;
#endif
}


template<TupleTypenames, typename ValueT> void
TupleKeyMap<TupleArgTypes,ValueT>::clear()
{
    m_rep.clear();
}

#undef Tuple2
#undef Tuple3
#undef Tuple4
#undef TupleTypenames
#undef TupleArgTypes
#undef TupleArgDecls

#include "NamespaceFooter.H"
#endif // _TUPLEKEYMAPI_H_

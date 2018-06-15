#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TupleKeyMapI.H"
#include "NamespaceHeader.H"

////////////////////////////////////////////////////////////////////
#ifdef DO_DEMO

struct Arg1T
{
  int m_i;

  Arg1T( int a_i )
    :
    m_i(a_i)
  {
  }

  bool operator<( Arg1T const & that ) const
  {
    return this->i < that.m_i;
  }

  bool operator==( Arg1T const & that ) const
  {
    return this->i == that.m_i;
  }
};

struct Arg2T
{
  int m_i;

  Arg2T( int a_i )
    :
    m_i(a_i)
  {
  }

  bool operator<( Arg2T const & that ) const
  {
    return this->i < that.m_i;
  }

  bool operator==( Arg2T const & that ) const
  {
    return this->i == that.m_i;
  }
};


struct Arg3T
{
  int m_i;

  Arg3T( int a_i )
    :
    m_i(a_i)
  {
  }

  bool operator<( Arg3T const & that ) const
  {
    return this->i < that.m_i;
  }

  bool operator==( Arg3T const & that ) const
  {
    return this->i == that.m_i;
  }
};


struct Arg4T
{
  int m_i;

  Arg4T( int a_i )
    :
    m_i(a_i)
  {
  }

  bool operator<( Arg4T const & that ) const
  {
    return this->i < that.m_i;
  }

  bool operator==( Arg4T const & that ) const
  {
    return this->i == that.m_i;
  }
};


int main()
{
    TupleKeyMap<Arg1T,Arg2T,Arg3T,Arg4T, char> mymap;
    mymap.insert( 1, 11, 9, 0, 'a');
    mymap.insert( 2, 11, 9, 0, 'b');
    mymap.insert( 1, 12, 9, 0, 'a');
    mymap.insert( 2, 12, 9, 0, 'c');
    mymap.insert( 1, 12, 9, 1, 'c');
    mymap.insert( 1, 12, 0, 0, 'c');
    mymap.insert( 4, 12, 9, 0, 'c');
    mymap.insert( 5, 17, 9, 0, 'c');
    mymap.report();
    // Should print: "unique arg combinations: 8"
}

#endif // DO_DEMO
////////////////////////////////////////////////////////////////////
#include "NamespaceFooter.H"

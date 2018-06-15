#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Callbacks.H"
#include <iostream>
#include "BaseNamespaceHeader.H"

void
Callbacks::add( void (*f)() )
{
  m_funcPointers.push_back( f );
}

void
Callbacks::run() const
{
  for ( std::vector<PVF>::const_iterator i = m_funcPointers.begin();
       i != m_funcPointers.end();
       ++i )
  {
    std::cerr << "Callbacks::run()\n";
    (**i)();
  }
}

#ifdef DO_DEMO

struct SourceOfCallback
{
  static void foo()
  {
    std::cerr << "SourceOfCallback::foo()\n";
  }
};

struct TargetOfCallback
{
  static Callbacks* s_callbacks;
  void addCallback( void (*f)() )
  {
    if ( ! s_callbacks )
    {
      s_callbacks = new Callbacks;
    }
    s_callbacks->add( f );
  }

  void doSomething()
  {
    s_callbacks->run();
  }
};

Callbacks* TargetOfCallback::s_callbacks(0);

int main()
{
  TargetOfCallback target;
  target.addCallback( SourceOfCallback::foo );
  target.doSomething();
}

#endif // DO_DEMO
#include "BaseNamespaceFooter.H"

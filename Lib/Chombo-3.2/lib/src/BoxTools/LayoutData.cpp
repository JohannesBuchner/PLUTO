#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LayoutData.H"

#include "NamespaceHeader.H"

// template < > void LayoutData<Real>::allocate()
// {
//  m_callDelete = true;

//   for (unsigned int i = 0; i < m_vector.size(); ++i)
//   {
//     delete m_vector[i];
//     m_vector[i] = NULL;
//   }

//   m_vector.resize(m_boxLayout.size(), NULL);

//   for (DataIterator it(dataIterator()); it.ok(); ++it)
//   {
//     unsigned int index = m_boxLayout.index(it());
//     if (m_vector[index] == NULL)
//     {
//       m_vector[index] = new Real;
//       if (m_vector[index] == NULL)
//       {
//         MayDay::Error("OutOfMemory in LayoutData::allocate");
//       }
//     }
//   }

// }

// template < > void LayoutData<int>::allocate()
// {
//  m_callDelete = true;

//   for (unsigned int i = 0; i < m_vector.size(); ++i)
//   {
//     delete m_vector[i];
//     m_vector[i] = NULL;
//   }

//   m_vector.resize(m_boxLayout.size(), NULL);

//   for (DataIterator it(dataIterator()); it.ok(); ++it)
//   {
//     unsigned int index = m_boxLayout.index(it());
//     if (m_vector[index] == NULL)
//     {
//       m_vector[index] = new int;
//       if (m_vector[index] == NULL)
//       {
//         MayDay::Error("OutOfMemory in LayoutData::allocate");
//       }
//     }
//   }

// }

#include "NamespaceFooter.H"

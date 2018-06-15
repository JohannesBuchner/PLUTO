#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DATAITERATOR_H_
#define _DATAITERATOR_H_

#include "Vector.H"
#include "DataIndex.H"
#include "BoxLayout.H"
#include "SPMD.H"
#include "LayoutIterator.H"
#include "NamespaceHeader.H"

#ifdef CH_MPI

///An Iterator based on a BoxLayout object.
/**
   An Iterator based on a BoxLayout object.  It does not
   support a dereferencing operation, since it is intended
   to work with all of BoxLayouts, DisjointBoxLayouts, BoxLayoutDatas
   LevelDatas, and any object that is built on top
   of a BoxLayout object.

   DataIterator accesses the data in a BoxLayout-based
   object in a <b> data-parallel </b> manner.
   This means that it skips over entries for Boxes not
   assigned to this processor.  The order of access is not defined.

   In serial execution mode, there is no difference between
   DataIterator and LayoutIterator.
*/
class DataIterator
{
  //Note: parallel version doesn't inherit LayoutIterator, serial version does
public:

  /// a null constructed DataIterator will return false on ok()
  DataIterator();

  DataIterator(const BoxLayout& a_layout)
  {
    *this = a_layout.dataIterator();
  }

  virtual ~DataIterator()
  {}

  /// return the index that this iterator is at
  /** Aborts if the iterator is not ok() */
  inline const DataIndex& operator()() const ;
  /// return a copy of the index that this iterator is at
  /** Aborts if the iterator is not ok() */
  DataIndex i() const
  {
    return this->operator()();
  }

  // const DataIndex& operator[](int ivec) const
  // {
  //   return (*m_indices)[ivec];
  // }

  // DataIndex& operator[](int ivec)
  // {
  //   return (*m_indices)[ivec];
  // }

  DataIndex operator[](int ivec) const
  {
    return (DataIndex)((*m_indices)[ivec]);
  }


  /// move the iterator to the next index in the layout
  inline virtual void operator++();
  /// move the iterator to the next index in the layout
  void incr()
  {
    ++(*this);
  }

  /// return true if this iterator is still in the layout
  virtual bool ok() const;

  /// initialize this iterator to the first index in the layout
  void begin();

  /// same as begin()
  void reset();

  /// move this iterator to after the last index in the layout
  /** The iterator will be !ok() afterwards. */
  void end();

  int size() const
  {
    return m_indices->size();
  }

private:

  friend class BoxLayout;
  friend class DisjointBoxLayout;
  friend class TimedDataIterator;

  DataIterator(const BoxLayout& boxlayout, const int* layoutID);

  BoxLayout  m_layout;
  //  RefCountedPtr<Vector<DataIndex> > m_indices;
protected:
  const Vector<DataIndex>* m_indices;

  int m_current;
};

inline DataIterator::DataIterator()
  : m_current(-1)
{
}

inline const DataIndex& DataIterator::operator()() const
{
  CH_assert(ok());
  return m_indices->operator[](m_current);
}

inline void DataIterator::operator++()
{
  ++m_current;
}

inline bool DataIterator::ok() const
{
  return m_current < m_indices->size();
}

inline void DataIterator::reset()
{
  begin();
}

inline void DataIterator::begin()
{
  m_current=0;
}

#else

// serial version

class DataIterator : public LayoutIterator
{
public:
  virtual ~DataIterator()
  {}

  DataIterator()
  {}

  DataIterator(const BoxLayout& a_layout)
  {
    *this = a_layout.dataIterator();
  }

  /// return the index that this iterator is at
  /** Aborts if the iterator is not ok() */
  const DataIndex& operator()() const
  {
    return (const DataIndex&)(LayoutIterator::operator()());
  };

  /// return a copy of the index that this iterator is at
  /** Aborts if the iterator is not ok() */
  DataIndex i() const
  {
    return this->operator()();
  }

  int size() const
  {
    return this->m_layout.size();
  }

  DataIndex operator[](int ivec) const
  {
    return (DataIndex)((*m_indicies)[ivec]);
  }

private:
  friend class BoxLayout;
  friend class DisjointBoxLayout;
  friend class TimedDataIterator;

protected:
  DataIterator(const BoxLayout& boxlayout, const int* layoutID)
    :LayoutIterator(boxlayout, layoutID)
  {}
};

#endif  /*CH_MPI*/

#define DATAITERATOR(CLASS, BOXLAYOUT)                \
        DataIterator dit = BOXLAYOUT .dataIterator(); \
        for (dit.begin(); dit.ok(); ++dit)            \
          {                                           \
            DataIndex di = dit();                     \
            MT_BEGIN1(CLASS, DataIndex, di)

#define ENDITERATOR(CLASS)                            \
            MT_END1(CLASS, DataIndex, di)             \
          }

#define DATAITERATOR1(CLASS, BOXLAYOUT, TYPE1, VAL1)     \
        DataIterator dit = BOXLAYOUT .dataIterator();    \
        for (dit.begin(); dit.ok(); ++dit)               \
          {                                              \
            DataIndex di = dit();                        \
            MT_BEGIN2(CLASS, TYPE1, VAL1, DataIndex, di)

#define ENDITERATOR1(CLASS, TYPE1, VAL1)                 \
            MT_END2(CLASS, TYPE1, VAL1, DataIndex, di)   \
          }

#define DATAITERATOR2(CLASS, BOXLAYOUT, TYPE1, VAL1, TYPE2, VAL2)     \
        DataIterator dit = BOXLAYOUT .dataIterator();                 \
        for (dit.begin(); dit.ok(); ++dit)                            \
          {                                                           \
            DataIndex di = dit();                                     \
            MT_BEGIN3(CLASS, TYPE1, VAL1, TYPE2, VAL2, DataIndex, di)

#define ENDITERATOR2(CLASS, TYPE1, VAL1, TYPE2, VAL2)                 \
            MT_END3(CLASS, TYPE1, VAL1, TYPE2, VAL2, DataIndex, di)   \
          }

#include "NamespaceFooter.H"
#endif

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <stdio.h>
#include <cstring>
#include "BitSet.H"
#include <cstdlib>
#include <cstdio>
#include "MayDay.H"
#include "SPMD.H"

#include "BaseNamespaceHeader.H"

BitSet BitSetIterator::emptyBitSet = BitSet();

BITSETWORD BitSet::trueMasks[BITSETWORDSIZE];

int BitSet::initialize()
{
  trueMasks[BITSETWORDSIZE-1] = 1;
  for (int i=BITSETWORDSIZE-2; i>=0 ; --i)
    {
      trueMasks[i] = trueMasks[i+1] << 1;
    }
  return 0;
}

BitSet::BitSet()
  :
  m_bits(NULL),
  m_size(0),
  m_length(0)
{
}

long int BitSet::bytes = BitSet::initialize();
long int BitSet::peak = 0;

int BitSet::linearSize() const
{
  return 2* CH_XD::linearSize<int>(0) + m_length* sizeof(BITSETWORD);
}

void BitSet::linearIn(const void* const inBuf)
{
  if (m_bits != NULL)
    {
      free(m_bits);
      bytes -= m_length*sizeof(BITSETWORD);
    }

  char* buf = (char*)inBuf;
  CH_XD::linearIn<int>(m_size, buf);
  buf+=CH_XD::linearSize<int>(m_size);
  CH_XD::linearIn<int>(m_length, buf);
  buf+=CH_XD::linearSize<int>(m_length);
  m_bits = (BITSETWORD*)malloc(m_length*sizeof(BITSETWORD));
  memcpy(m_bits, buf, m_length*sizeof(BITSETWORD));
  bytes += m_length*sizeof(BITSETWORD);

}

void BitSet::linearOut(void* const a_outBuf) const
{
  char* buf = (char*)a_outBuf;
  CH_XD::linearOut<int>(buf, m_size);
  buf+=CH_XD::linearSize<int>(m_size);
  CH_XD::linearOut<int>(buf, m_length);
  buf+=CH_XD::linearSize<int>(m_length);
  memcpy(buf, m_bits, m_length*sizeof(BITSETWORD));
}

BitSet::BitSet(int bits, bool init)
  :
  m_bits(NULL),
  m_size(bits),
  m_length(0)
{
  define(bits, init);
}

bool BitSet::operator<( const BitSet& rhs ) const
{
  if ( m_length < rhs.m_length )
  {
    return true;
  } else
  if ( m_length > rhs.m_length )
  {
    return false;
  }

  for ( int i=0; i<m_length; ++i )
  {
    if ( m_bits[i] < rhs.m_bits[i] )
    {
      return true;
    } else
    if ( m_bits[i] > rhs.m_bits[i] )
    {
      return false;
    }
  }

  return false;
}

void BitSet::define(int bits, bool init)
{
  if (m_bits != NULL)
  {
    free(m_bits);
    bytes-=m_length*sizeof(BITSETWORD);
  }
  m_size = bits;
  CH_assert(bits >= 0);
  m_length = bits/BITSETWORDSIZE;
  if (bits == 0)
  {
    m_bits = NULL; return;
  }
  if (m_length*BITSETWORDSIZE != bits) m_length +=1;
  m_bits = (BITSETWORD*)malloc(m_length*sizeof(BITSETWORD));
  bytes+= m_length*sizeof(BITSETWORD);
  if (bytes > peak) peak = bytes;

  if (m_bits == NULL)
    {
      char mesg[1024];
      sprintf(mesg, "Memory Error in BitSet::BitSet(int bits=%d, bool init)", bits);
      MayDay::Error(mesg);
    }
  if (init)
    {
      setAllTrue();
    }
  else
    {
      setAllFalse();
    }
}

BitSet::BitSet(const BitSet& rhs)
  :
  m_bits(NULL),
  m_size(rhs.m_size),
  m_length(rhs.m_length)
{
  if (m_length == 0) return;
  m_bits = (BITSETWORD*)malloc(m_length*sizeof(BITSETWORD));
  bytes+= m_length*sizeof(BITSETWORD);
  if (bytes > peak) peak = bytes;
  if (m_bits == NULL)
    {
      MayDay::Error("Memory Error in BitSet::BitSet(const BitSet& rhs)");
    }

  memcpy(m_bits, rhs.m_bits, m_length*sizeof(BITSETWORD));
}

bool BitSet::isEmpty() const
{
  if (m_length == 0) return true;
  for (int i=0; i<m_length; ++i)
    {
      if (m_bits[i] > 0) return false;
    }
  /*
  for (int i=(m_length-1)*BITSETWORDSIZE; i<m_size; ++i)
    {
      if (this->operator[](i)) return false;
    }
  */
  return true;
}

bool BitSet::isFull() const
{
  BITSETWORD g=0;
  g=~g;
  for (int i=0; i<m_length-1; ++i)
    {
      if (m_bits[i] != g) return false;
    }
  for (int i=(m_length-1)*BITSETWORDSIZE; i<m_size; ++i)
    {
      if (!this->operator[](i)) return false;
    }
  return true;
}

BitSet& BitSet::operator=(const BitSet& rhs)
{
  if (m_bits != NULL)
    {
      free(m_bits);
      bytes -= m_length*sizeof(BITSETWORD);
      m_bits = NULL;
    }
  m_length=rhs.m_length;
  m_size =rhs.m_size;

  if (rhs.m_bits != NULL)
    {
      m_bits = (BITSETWORD*)malloc(m_length*sizeof(BITSETWORD));
      bytes+= m_length*sizeof(BITSETWORD);
      if (bytes > peak) peak = bytes;
    }
  else
    {
      m_bits = NULL;
      return *this;
    }
  if (m_bits == NULL)
    {
      MayDay::Error("Memory Error in BitSet::operator=(const BitSet& rhs)");
    }

  memcpy(m_bits, rhs.m_bits, m_length*sizeof(BITSETWORD));
  return *this;
}

BitSet::~BitSet()
{
  if (m_bits != NULL)
    {
      free(m_bits);
      bytes-=m_length*sizeof(BITSETWORD);
    }
}

void BitSet::setAllTrue()
{
  memset(m_bits, ~0, (m_length-1)*sizeof(BITSETWORD));
  BITSETWORD* word = m_bits+(m_length-1);
  *word = 0;
  int count = m_size - (m_length-1)*BITSETWORDSIZE;
  for (int i=0; i<count; ++i)
    {
      *word = *word | trueMasks[i];
    }
}

void BitSet::setAllFalse()
{
  memset(m_bits, 0, m_length*sizeof(BITSETWORD));
}

void BitSet::setFalse(int i)
{
  // BITSETWORD* tm = trueMasks;
  CH_assert(i>=0);
  CH_assert(i<m_size);
  int index = i/BITSETWORDSIZE;
  CH_assert(index < m_length);
  int remainder = i-BITSETWORDSIZE*index;
  BITSETWORD* word = m_bits+index;

  *word = *word & ~(trueMasks[remainder]);
}

void BitSet::setTrue(int i)
{
  CH_assert(i>=0);
  CH_assert(i<m_size);
  int index = i/BITSETWORDSIZE;
  CH_assert(index < m_length);
  int remainder = i-BITSETWORDSIZE*index;
  BITSETWORD* word = m_bits+index;

  *word = *word | trueMasks[remainder];
}

void BitSetIterator::end()
{
  m_index = m_length;
  m_remainder = m_partialBits;
}
void BitSetIterator::begin()
{
  m_index = 0;
  m_remainder = 0;
}
#include "BaseNamespaceFooter.H"

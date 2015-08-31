/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file
  \brief  Internal routine to sort an integer array

  Internal routine to sort an integer array
 
  \author A. Malagoli (University of Chicago)
  \date Jul 17, 1999  
*/
/* ///////////////////////////////////////////////////////////////////// */
#define swapi_(a,b)  ( temp=(a), (a)=(b), (b)=temp)

/* ********************************************************************* */
int AL_Sort_(int n, int *in, int *ind)
/*!
 * Sort and array of integers
 *           
 * This is a really simple implementation, since we do not really
 *  use this for large arrays.
 *
 * \param [in] n    size of input array (integer)
 * \param [in] in   input array
 * \param [in] ind  array of the sorted index arrays (max to min)
 *********************************************************************** */
{
  register int i, temp;
  int m, r, s;
  
  /* Start with the normal index ordering */
  for(i=0;i<n;i++) ind[i]=i;

  m=0;
  while(m<n) {
    r = m;
    s = in[m];
    for(i=m+1;i<n;i++) if( in[ind[i]] > s ) { r=i; s=in[ind[i]]; }
    swapi_(ind[m],ind[r]);
    m++;
  }
  return 0;
}

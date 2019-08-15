/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Documentation template for C source files.

  Detailed description of the file goes here.         

  \authors A. Mignone (mignone@ph.unito.it)\n
           C. Zanni   (zanni@oato.inaf.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
           T. Matsakos

 \b References
    - "PLUTO: A Numerical Code for Computational Astrophysics." \n
      Mignone et al, ApJS (2007) 170, 228
    - "The PLUTO Code for Adaptive Mesh Computations in Astrophysical Fluid
       Dynamics" \n
      Mignone et al, ApJS (2012) 198, 7M
    - "A conservative orbital advection scheme for simulations of magnetized
       shear flows with the PLUTO code"\n
      Mignone et al., A\&A (2012) 545A, 152M

  \date   Aug 20, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/*! This is a comment for the following macro */
#define NEW_MACRO(a)  (a+1)

int global_var; /**< This is the description of global_var1 */

/* ********************************************************************* */
void MyFunctionTemplate (int var)
/*!
 * Start commenting the function here.
 *
 * \param [in]      sweep pointer to Sweep structure
 * \param [in]      beg   initial index of computation 
 * \param [in]      end   final   index of computation
 * \param [in,out]  d     pointer to PLUTO Data structure
 * \param [out]     src   2D array of source terms
 * \param [in]      grid  pointer to an array of Grid structures
 * 
 * The following is an itemized list of HOWTO:\n
 * - Display an extended math formula:
 *   \f[
 *         A\cos\phi = B\exp\left(\frac{y}{q}\right)
 *   \f]
 *   while in the text simply use \f$ x^2-y^2 = 0\f$.
 * - To display something in bold use \b bold symbol.
 * - To display something in typewrite use \c font symbol
 *   (for multiple words use <tt>more than one word</tt>).
 * - To emphasize text use *text* or **text** or, alternatively,
 *   \e text (for multiple words <em> more than one word </em>).
 * - To reference to a global variable use ::global_var.
 * - To reference to a function use just, e.g., Boundary().
 * 
 * Indentation stops here. \n
 *
 * To produce an enumerated list:
 *
 * -#  first item
 * -#  second item
 * -# etc...
 *
 * 
 * To produce code block just leave one blank line before and after and
 * then indent like
 *
 *      FL -> swL*FL    + swR*I(FR)
 *      FR -> swL*I(FL) + swR*FR
 *
 * and come back.
 *
 * To display a piece of code use
 * \code
 *   restart;
 *   a := b;
 *   c := d-e;
 * \endcode
 *
 *
 *  
 * \b References
 *    - ".." \n
 *      Mignone et al, JCP (2010) xx, xx
 *
 * \return  This function has no return value.
 * \see
 * \attention
 * \note  --
 * \todo  --
 *
 *********************************************************************** */
{

/* ---------------------------------------------------------------- */
/*! This is a buil-int comment                                      */
/* ---------------------------------------------------------------- */

}


/* ********************************************************************* */
/*! Structure comment here
   ********************************************************************* */

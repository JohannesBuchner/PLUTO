/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief File parser utilities.

  This file provides a set of useful functions to open / read and
  parse the content of a parameter file (typically pluto.ini).
  The parameter file can contain a number of lines
  with the general structure

  <tt> label    value1   value2   ...  </tt>
 
  where the number of values following the label can be different
  for each line.

  ParamFileRead() reads the file and store its content into an array of
  lines (\c **fline).
  ParamFileGet() can be used to retrieve the n-th parameter value
  following a given label, while ParQuery() check whether a
  parameter actually exists.
  
  As an example consider the following file "myparam.txt":
  
  \verbatim
   ---- File myparam.txt ----- 
   nx        100                
   xdomain   15.0  30.0         
   ---------------------------
  \endverbatim

  The parameter can be read with the following code snippet:
  \code
   int nlines,nx;
   double xbeg, xend;
   nlines = ParamFileRead("myparam.txt");
   nx     = atoi(ParamFileGet("nx", 1));
   xbeg   = atof(ParamFileGet("xdomain", 1));
   xend   = atof(ParamFileGet("xdomain", 2));
  \endcode
  
  \authors A. Mignone (mignone@ph.unito.it)
  \date    June 18, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static int nlines;   /**< The total number of lines (including empty ones)
                          contained in the file */
static char **fline; /**< All of the lines (including empty ones) in the file */
static int ParamFileGetWords(char *line, char **);

/* ********************************************************************* */
int ParamFileRead (char *fname)
/*!
 *  Parse file *fname and store its content line by line in *fline.
 *  Blank lines are excluded.
 *
 * \param [in]  fname  the name of the file to be read
 * \return the number of liens successfully read.
 *********************************************************************** */
{
  char  sline[512];
  FILE *fp;

  if (fline == NULL) fline = ARRAY_2D(128,128,char);

  fp = fopen(fname,"r");
  if (fp == NULL) {
    printf ("! ParamFileRead: file %s not found\n", fname);
    QUIT_PLUTO(1);
  }
  nlines = 0;

  while ( fgets(sline, 512, fp) != NULL ) {
    if (strlen(sline) > 0) {
      strcpy (fline[nlines],sline); 
      nlines++;
    }
  }
  fclose(fp);

  return(nlines);
}

/* ********************************************************************* */
char *ParamFileGet (const char *label, int pos)
/*!
 * Search for \c *label in all the lines pointed to by \c **fline. 
 * If label exists, return a pointer to the string located pos
 * words after label. Issue an error if this string cannot
 * be located.
 *
 * \param [in]  label  the first word of a line to be searched
 * \param [in]  pos    an integer giving the position of the word
 * \return the n-th word contained in the same line
 *********************************************************************** */
{
  int         k, nwords, nw;
  static char **words;

/* -------------------------------------------------------------
    Allocate memory and initialize words[][] string array with
    to avoid unpleasant situations occurring when there're
    empty words.
   ------------------------------------------------------------- */

  if (words == NULL) words = ARRAY_2D(128,128,char);
  for (k = 0; k < 127; k++) sprintf (words[k],"\0");

  for (k = 0; k < nlines; k++) {        /* Loop over lines */
    nwords = ParamFileGetWords(fline[k],words);  /* get words for k-th line */

    if (nwords > 0 && strcmp(words[0],label) == 0){  /* check if 1st word  */
      if (pos <= nwords) return words[pos];          /* matches label      */
      else {
        printf ("! ParamFileGet: field # %d does not exist\n",pos);
        QUIT_PLUTO(1);
      }
    }
  }

  printf ("! ParamFileGet: label '%s' was not found\n",label);
  QUIT_PLUTO(1);

  return NULL;    
}

/* ********************************************************************* */
int ParamFileHasBoth (const char *label1, const char *label2)
/*!
 * Locate the line beginning with label1 and return 1 if label2 
 * can be found on the same line. Return 0 otherwise.
 *
 * \param [in] label1  The first word of the line to be searched
 * \param [in] label2  A word containined in the same line beginning with 
 *                     label1
 * \return 1 If label2 and label1 are in the same line. 0 otherwise.
 *********************************************************************** */
{ 
  int         k, nwords, nw;
  static char **words;

  if (words == NULL) words = ARRAY_2D(128,128,char);

  for (k = 0; k < nlines; k++) {                /* Loop over lines */
    nwords = ParamFileGetWords(fline[k],words); /* get words for k-th line */
    if (strcmp(words[0],label1) == 0){          
      for (nw = 1; nw < nwords; nw++){
        if (words[nw][0] == '#') break;  /* comment detected:
                                            no need to read any further */
        if (strcmp(words[nw], label2) == 0) return 1;
      }
    }
  }
  return 0;
}

/* ********************************************************************* */
int ParamExist (const char *label)
/*!
 * Check whether *label exists in any of the lines (**fline).
 *
 * \param [in]  label  
 * \return 0 on success, 1 if label cannot be found.
 *********************************************************************** */
{
  int    k;
  char   sline[512], *str;
  const char  delimiters[] = " \t\r\f";

/* ---------------------------------------
     search if label exists
   --------------------------------------- */

  for (k = 0; k < nlines; k++) {
    sprintf (sline,"%s",fline[k]);
    str = strtok(sline,delimiters);
    if (strcmp(str,label) == 0) return(1);
  }

  return (0);
}

/* ********************************************************************* */
int ParamFileGetWords(char *line, char **words)
/*!
 *  Return the words and their number contained in a single line.
 *
 *********************************************************************** */
{
  int    nw, nwords=0;
  char        *str, cstr[512];
  const char  delimiters[] = " \t\r\f";

/* -- make a local copy of 'line' since strtok will modify it -- */

  strcpy (cstr,line);   
  str = strtok(cstr,delimiters);  

/* -- loop until end of line, copy and count words -- */

  while (str != NULL && str[0] != '\n') {
    strcpy(words[nwords],str);
    str = strtok(NULL, delimiters);
    nwords++;
  }

/* -- get rid of newline character at the end of each word (if any) -- */

  for (nw = 0; nw < nwords; nw++){
    strcpy (cstr, words[nw]);
    sprintf (words[nw],"%s",strtok(cstr,"\n")); 
  }

  return nw;
}

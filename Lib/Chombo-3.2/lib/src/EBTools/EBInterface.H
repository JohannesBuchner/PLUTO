#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*  Serial programming interface to Chombo EB HDF5 file format*/

#ifndef _EBINTERFACE_H_
#define _EBINTERFACE_H_

#include <hdf5.h>
#include <math.h>
#include "NamespaceHeader.H"

struct
{
  int i;
  int j;
} typedef intvect2d ;

struct
{
  int i;
  int j;
  int k;
} typedef intvect3d ;

struct
{
  intvect2d lo;
  intvect2d hi;
} typedef box2d ;

struct
{
  intvect3d lo;
  intvect3d hi;
} typedef box3d;

union
{
  box2d b2;
  box3d b3;
} typedef box;

struct
{
  intvect3d index;
  double    volfrac;
  double    areafrac[6];
} typedef regvof ; /* regular volume-of-fluid */

struct
{
  double    areafrac;
  intvect3d neighbour;
  int       nindex;
} typedef irregface;

struct
{
  intvect3d index;
  double volfrac;
  int numfaces;
  irregface* faces;
} typedef irregvof;  /* irregular volume-of-fluid */

/* gives values to *_id objects */
/* must be called, and called only once, before other functions in this
   header will function properly */
void initializeHDF5datatypes();

/*  domain is broken up into many non-overlapping boxes that completely cover
        the domain.

        for each box, 'i', there is a regoffset, and an irregoffset:

        if regoffset[i] = -1  the box is completely REGULAR
        if regoffset[i] = -2  the box is completely COVERED

        otherwise regoffset contains an index into regvof* array where reading
        is to begin of regular volume-of-fluids,
        regvof[ regoffset [i]] to regvof[ regoffset [i+1]]

        irregular volumes are then written for each box,
        irregvof[ irregoffset [i]] to irregvof[ irregoffset [i+1]]

        both indexings do not include the endpoint [...) , so that emtpy
        sets can be given with both indexes being the same.

        return codes;

        0 : success
        -1: error opening file
        -2: error writing attributes to header
        -4: error writing boxes to data file
        -20: error closing file on exit
*/

int  writeEBChomboFile(const char* filename,
                                           box3d domain,
                                           int length,
                                           box3d* boxes,
                                           long*  regoffset,
                                           long*  irregoffset,
                                           long numreg,
                                           long numirreg,
                                           regvof* regularVofs,
                                           irregvof* irregularVofs);

/*================================*/
/*  END OF USER INTERFACE SECTION */
/*================================*/

/*  these functions prototypes are used by EBInterface.c */

int isEmpty(const  box2d*);
int numPnts2(const box2d*);
int numPnts3(const box3d*);

void grow2(box2d*, intvect2d*);
void grow3(box3d*, intvect3d*);

void refine2(box2d*, int refinement);
void refine3(box3d*, int refinement);
int cut(const box3d* box, int dir, int index);
box2d* crossSection(const box3d* box, int dir, box2d* box2);

static hid_t intvect2d_id, intvect3d_id, box2d_id, box3d_id ;

enum
{
  Float,
  Double
} typedef ChPrecision;

struct
{
  hid_t file_ID;
  hid_t group_ID;
  int   dim;
  ChPrecision precision;
}typedef HDF5Handle;

/* three next functions return 0 on success */

/* returns negative number if file open failed
           1  if 'Chombo_global' couldn't be found and ChDim and ChPrecision not set
   accessMode is one of H5F_ACC_RDWR, H5F_ACC_RDONLY, H5F_ACC_CREATE
           */
int Handleopen(HDF5Handle* handle, const char* filename, hid_t accessMode);
int HandlesetGroup(HDF5Handle* handle, const char* group);/* ie. "/level_0" */
int Handleclose(HDF5Handle handle);
int HandleCreateGroup(HDF5Handle* handle, const char* group);/* ie. "/level_0" */

#define ChTYPES 8
enum
{
  INTEGER,
  FLOAT,
  DOUBLE,
  CHAR,
  INTVECT2D,
  INTVECT3D,
  BOX2D,
  BOX3D
} typedef ChattributeType;
static const char* TYPE_NAMES[ChTYPES];

struct
{
  void* data;
  char* name;
  int attributeType;
  void* next;
} typedef Attribute;

/*  struct that holds all the attributes of an HDF5 group.

  each entry in the accessByType array is the head of a linked
  list of members terminated with a null.

  to access directly by type you can use code like this:

   HDF5attributes a;
   readHDF5attributes(&a, handle);
   int numChar = a.numByType(CHAR);
   Attribute* at = a.accessByType[CHAR];
   while (at != NULL)
   {
      printf("char attribute %s = %s",
      at->name,
      (char*)(at->data));
      at = at->next;
   }
   */
struct
{
  int numByType[ChTYPES];
  Attribute*  accessByType[ChTYPES];
  int   dim;
  ChPrecision precision;
} typedef HDF5attributes;

void freeHDF5attributes(HDF5attributes* attributes);
void printAttributes(HDF5attributes* attributes);

/* utility functions for appending attributes to an HDF5attributes struct
   attrib:  HDF5atributes struct */
extern Attribute* ADD_tmp;
extern Attribute* ADD_at;
#define ADD_ATTRIBUTE(attrib, Ttype, ChType, Sname, value)         \
ADD_tmp = attrib.accessByType[ChType] ;                            \
ADD_at = (Attribute*) malloc(sizeof(Attribute));                   \
attrib.accessByType[ChType] = ADD_at;                              \
ADD_at->name = malloc(strlen(Sname)+1);                              \
sprintf(ADD_at->name, Sname);                                      \
ADD_at->data = malloc(sizeof(Ttype));                              \
ADD_at->attributeType = ChType;                                    \
*((Ttype *)ADD_at->data) = value;                                  \
  ADD_at->next = ADD_tmp;                                          \
attrib.numByType[ChType]++;                                        \

#define ADD_CHAR_ATTRIBUTE(attrib,  Sname, value)                  \
ADD_tmp = attrib.accessByType[CHAR] ;                              \
ADD_at = (Attribute*) malloc(sizeof(Attribute));                   \
attrib.accessByType[CHAR] = ADD_at;                                \
ADD_at->name = malloc(strlen(Sname));                              \
sprintf(ADD_at->name, Sname);                                      \
ADD_at->data = malloc(strlen(value));                              \
ADD_at->attributeType = CHAR;                                      \
sprintf((char*)ADD_at->data, value);                               \
  ADD_at->next = ADD_tmp;                                          \
attrib.numByType[CHAR]++;                                          \

/* next eight functions return 0 on success */

/* - users responsibility to make sure HDF5Handle is pointing at a valid group
   in an open file.
   - User must use freeHDF5attributes when finished.
   - User must pass in a pointer to a valid HDF5attributes struct.*/
int readHDF5attributes(HDF5attributes* attr, HDF5Handle handle);

/* input: HDF5Handle
   users responsibility to cast *(boxes) to either box2d*
   or box3d* and free it */
int readBoxes(box** boxes, int* length, HDF5Handle handle);

/* - users responsibility to make sure HDF5Handle is pointing at a valid group
   in an open file.
   - User must use freeHDF5attributes when finished.
   - User must pass in a pointer to a valid HDF5attributes struct.

   returns 0 on success

   -1 attribute creation problems  */
int writeHDF5attributes(HDF5attributes* attr, HDF5Handle handle);

int writeBoxes(box* boxes,    int length, HDF5Handle handle);

/** create our own platform independent form that spans
 a lot of IEEE goofiness. */
#define CFINITE(X) !isnan(X) && (X) < HUGE_VAL && (X) > -HUGE_VAL

#include "NamespaceFooter.H"
#endif

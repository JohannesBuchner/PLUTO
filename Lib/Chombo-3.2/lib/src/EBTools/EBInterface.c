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

#include "EBInterface.H"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

Attribute* ADD_tmp;
Attribute* ADD_at;

static const char* TYPE_NAMES[ChTYPES] =
{
  "INTEGER",
  "FLOAT",
  "DOUBLE",
  "CHAR",
  "INTVECT2D",
  "INTVECT3D",
  "BOX2D",
  "BOX3D"
};

void initializeHDF5datatypes()
{
  static int initialized = 1;
  if (initialized == 1)
  {
  initialized = 0;
  intvect2d_id = H5Tcreate (H5T_COMPOUND, sizeof(intvect2d));
  H5Tinsert (intvect2d_id, "intvecti", HOFFSET(intvect2d, i), H5T_NATIVE_INT);
  H5Tinsert (intvect2d_id, "intvectj", HOFFSET(intvect2d, j), H5T_NATIVE_INT);

  intvect3d_id = H5Tcreate (H5T_COMPOUND, sizeof(intvect3d));
  H5Tinsert (intvect3d_id, "intvecti", HOFFSET(intvect3d, i), H5T_NATIVE_INT);
  H5Tinsert (intvect3d_id, "intvectj", HOFFSET(intvect3d, j), H5T_NATIVE_INT);
  H5Tinsert (intvect3d_id, "intvectk", HOFFSET(intvect3d, k), H5T_NATIVE_INT);

  /* old composite-of-composite style boxes
  box2d_id = H5Tcreate (H5T_COMPOUND, sizeof(box));
  H5Tinsert (box2d_id, "smallend", HOFFSET(box2d, lo), intvect2d_id);
  H5Tinsert (box2d_id, "bigend",   HOFFSET(box2d, hi), intvect2d_id);

  box3d_id = H5Tcreate (H5T_COMPOUND, sizeof(box));
  H5Tinsert (box3d_id, "smallend", HOFFSET(box3d, lo), intvect3d_id);
  H5Tinsert (box3d_id, "bigend",   HOFFSET(box3d, hi), intvect3d_id);
  */
  box2d_id = H5Tcreate (H5T_COMPOUND, sizeof(box));
  H5Tinsert (box2d_id, "lo_i", HOFFSET(box2d, lo.i), H5T_NATIVE_INT);
  H5Tinsert (box2d_id, "lo_j", HOFFSET(box2d, lo.j), H5T_NATIVE_INT);
  H5Tinsert (box2d_id, "hi_i", HOFFSET(box2d, hi.i), H5T_NATIVE_INT);
  H5Tinsert (box2d_id, "hi_j", HOFFSET(box2d, hi.j), H5T_NATIVE_INT);

  box3d_id = H5Tcreate (H5T_COMPOUND, sizeof(box));
  H5Tinsert (box3d_id, "lo_i", HOFFSET(box3d, lo.i), H5T_NATIVE_INT);
  H5Tinsert (box3d_id, "lo_j", HOFFSET(box3d, lo.j), H5T_NATIVE_INT);
  H5Tinsert (box3d_id, "lo_k", HOFFSET(box3d, lo.k), H5T_NATIVE_INT);
  H5Tinsert (box3d_id, "hi_i", HOFFSET(box3d, hi.i), H5T_NATIVE_INT);
  H5Tinsert (box3d_id, "hi_j", HOFFSET(box3d, hi.j), H5T_NATIVE_INT);
  H5Tinsert (box3d_id, "hi_k", HOFFSET(box3d, hi.k), H5T_NATIVE_INT);
  }
}

void printAttributes(HDF5attributes* attr)
{
  Attribute* at;
  int i;
  for (i=0; i<ChTYPES; ++i)
    {
      at = attr->accessByType[i];
      while (at != NULL)
    {
      printf("%s attribute %s \n",
        TYPE_NAMES[i],
         at->name);
      at = at->next;
    }
    }
}

void initHDF5attributes(HDF5attributes* attr)
{
  int i;
  for (i=0; i<ChTYPES; ++i)
    {
       attr->accessByType[i] = NULL;
       attr->numByType[i] = 0;
    }
}

int numPnts2(const box2d* box)
{
  return (box->hi.i-box->lo.i+1)*(box->hi.j-box->lo.j+1);
}

int numPnts3(const box3d* box)
{
  return (box->hi.i-box->lo.i+1)*(box->hi.j-box->lo.j+1)*(box->hi.k-box->lo.k+1);
}

int isEmpty(const box2d* box)
{
  return (box->hi.i < box->lo.i)||(box->hi.j < box->lo.j);
}

void refine2(box2d* b2, int refinement)
{
  b2->lo.i*=refinement;
  b2->lo.j*=refinement;
  b2->hi.i = (b2->hi.i + 1)*refinement - 1;
  b2->hi.j = (b2->hi.j + 1)*refinement - 1;
}

void refine3(box3d* b3, int refinement)
{
  b3->lo.i*=refinement;
  b3->lo.j*=refinement;
  b3->lo.k*=refinement;
  b3->hi.i = (b3->hi.i + 1)*refinement - 1;
  b3->hi.j = (b3->hi.j + 1)*refinement - 1;
  b3->hi.k = (b3->hi.k + 1)*refinement - 1;
}

void grow2(box2d* b2, intvect2d* iv)
{
  if (iv==NULL) return;
  b2->lo.i-=iv->i;
  b2->lo.j-=iv->j;
  b2->hi.i+=iv->i;
  b2->hi.j+=iv->j;
  return;
}

void grow3(box3d* b3, intvect3d* iv)
{
  if (iv==NULL) return;
  b3->lo.i-=iv->i;
  b3->lo.j-=iv->j;
  b3->lo.k-=iv->k;
  b3->hi.i+=iv->i;
  b3->hi.j+=iv->j;
  b3->hi.k+=iv->k;
  return;
}

int Handleopen(HDF5Handle* handle, const char* filename, hid_t access)
{
  herr_t ret = 0;
  hid_t attr, datatype, aid;
  float fl = 0;
  double dl = 0;

  if (access == H5F_ACC_RDONLY)
    {
      handle->file_ID = H5Fopen(filename, access, H5P_DEFAULT);
      if (handle->file_ID < 0) return -1;
      handle->group_ID  =  H5Gopen(handle->file_ID, "ChomboGlobal");
      if (handle->group_ID < 0) return 1;
      /*      printf("handle->file_ID=%i, handle->group_ID=%i",handle->file_ID, handle->group_ID); */
      attr = H5Aopen_name(handle->group_ID, "SpaceDim");
      ret  = H5Aread(attr, H5T_NATIVE_INT, &(handle->dim));
      ret  = H5Aclose(attr);

      attr = H5Aopen_name(handle->group_ID, "testReal");
      datatype = H5Aget_type(attr);
      if (H5Tget_precision(datatype) == H5Tget_precision(H5T_NATIVE_FLOAT))
    handle->precision = Float;
      else
    handle->precision = Double;
    }
  else if (access == H5F_ACC_TRUNC)
    {
      handle->file_ID = H5Fcreate(filename,  H5F_ACC_TRUNC,  H5P_DEFAULT, H5P_DEFAULT);
      handle->group_ID  =  H5Gcreate(handle->file_ID,  "ChomboGlobal", 0);
      aid  = H5Screate(H5S_SCALAR);
      attr = H5Acreate( handle->group_ID, "SpaceDim", H5T_NATIVE_INT, aid, H5P_DEFAULT);
      ret = H5Awrite(attr, H5T_NATIVE_INT, &(handle->dim));
      H5Aclose(attr);
      if (handle->precision == Float)
    {

      attr =  H5Acreate(handle->group_ID, "testReal",
                H5T_NATIVE_FLOAT, aid, H5P_DEFAULT);
      ret = H5Awrite(attr, H5T_NATIVE_FLOAT, &fl);
    }
      else
    {
      attr =  H5Acreate(handle->group_ID, "testReal",
                H5T_NATIVE_DOUBLE, aid, H5P_DEFAULT);
      ret = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dl);
    }
    }
  if (handle->group_ID > 0) H5Gclose(handle->group_ID);
  handle->group_ID = H5Gopen(handle->file_ID, "/");
  return 0;
}

int HandlesetGroup(HDF5Handle* handle, const char* group)
{
  H5Gclose(handle->group_ID);
  handle->group_ID = H5Gopen(handle->file_ID, group);
  if (handle->group_ID > 0) return 0;
  return -1;
}

int HandleCreateGroup(HDF5Handle* handle, const char* group)
{
  H5Gclose(handle->group_ID);
  handle->group_ID = H5Gcreate(handle->file_ID, group, 0);
  if (handle->group_ID > 0) return 0;
  return -1;
}

int Handleclose(HDF5Handle handle)
{
  H5Gclose(handle.group_ID);
  return H5Fclose(handle.file_ID);
}

void freeHDF5attributes(HDF5attributes* attributes)
{
  int i;  Attribute *next, *tmp;
  for (i=0; i<ChTYPES; ++i)
    {
      next = attributes->accessByType[i];
      while (next != NULL)
    {
      tmp = next;
      next=tmp->next;
      free(tmp->data);
      free(tmp->name);
      free(tmp);
    }
      attributes->accessByType[i] = NULL;
      attributes->numByType[i]    = 0;
    }

}

herr_t attributeScan(hid_t loc_id, const char *name, void *opdata)
{
  herr_t ret = 0;
  HDF5attributes* data = opdata;
  hid_t attr, atype, aclass;
  size_t size = 0;
  Attribute* attribute, *head;

  size = strlen(name) + 1;
  attribute = malloc(sizeof(Attribute));
  attribute->name = malloc(size);
  memcpy(attribute->name, name, size);

  attr   = H5Aopen_name(loc_id, name);
  atype  = H5Aget_type(attr);
  aclass = H5Tget_class(atype);

  switch(aclass)
  {
  case H5T_INTEGER :
    data->numByType[INTEGER]++;
    attribute->attributeType = INTEGER;
    attribute->data = malloc(sizeof(int));
    ret = H5Aread(attr, H5T_NATIVE_INT, attribute->data);
    break;
  case H5T_FLOAT:
    if (data->precision == Float)
      {
    data->numByType[FLOAT]++;
    attribute->attributeType = FLOAT;
    attribute->data = malloc(sizeof(float));
    ret = H5Aread(attr, H5T_NATIVE_FLOAT, attribute->data);
      }
    else
      {
    data->numByType[DOUBLE]++;
    attribute->attributeType = DOUBLE;
    attribute->data = malloc(sizeof(double));
    ret = H5Aread(attr, H5T_NATIVE_DOUBLE, attribute->data);
      }
    break;
  case H5T_STRING:
    data->numByType[CHAR]++;
    size = H5Tget_size(atype);
    attribute->attributeType = CHAR;
    attribute->data = malloc(size+1);
    ret = H5Aread(attr, atype, attribute->data );
    if (ret < 0) break;
    ((char *)(attribute->data))[size] = '\0';
    break;
  case H5T_COMPOUND:
    if (strcmp(H5Tget_member_name(atype, 0), "lo_i") == 0)
      {
    if (data->dim == 2)
      {
        data->numByType[BOX2D]++;
        attribute->attributeType = BOX2D;
        attribute->data = malloc(sizeof(box));
        ret = H5Aread(attr, box2d_id, attribute->data);
      }
    else if (data->dim == 3)
      {
        data->numByType[BOX3D]++;
        attribute->attributeType = BOX3D;
        attribute->data = malloc(sizeof(box));
        ret = H5Aread(attr, box3d_id, attribute->data);
      }
    break;
      }
    else if (strcmp(H5Tget_member_name(atype, 0), "intvecti") == 0)
      {
    if (data->dim == 2)
      {
        data->numByType[INTVECT2D]++;
        attribute->attributeType = INTVECT2D;
        attribute->data = malloc(sizeof(intvect2d));
        ret = H5Aread(attr, intvect2d_id, attribute->data);
      }
    else if (data->dim == 3)
      {
        data->numByType[INTVECT3D]++;
        attribute->attributeType = INTVECT3D;
        attribute->data = malloc(sizeof(intvect3d));
        ret = H5Aread(attr, intvect3d_id, attribute->data);
      }
    break;
      }
  default:
    /* don't know what the hell this thing is */
    free(attribute->name);
    free(attribute);
    return -1;
  }

  /* OK, lets tack this attribute to the right linked-list */
  head = data->accessByType[attribute->attributeType];
  data->accessByType[attribute->attributeType] = attribute;
  attribute->next = head;
  return ret;
}

int readHDF5attributes(HDF5attributes* attr, HDF5Handle handle)
{
  int i;

  attr->dim = handle.dim;
  attr->precision = handle.precision;
  for (i=0; i<ChTYPES; ++i)
    {
      attr->numByType[i] = 0;
      attr->accessByType[i] = NULL;
    }
  return H5Aiterate(handle.group_ID, NULL, attributeScan, attr);
}

int readBoxes(box** boxes, int* length, HDF5Handle handle)
{
  hid_t boxdataset, boxdataspace, memdataspace;
  hsize_t dims[1], maxdims[1];
  herr_t error;

  boxdataset = H5Dopen(handle.group_ID, "boxes");
  if (boxdataset < 0) return boxdataset;
  boxdataspace =  H5Dget_space(boxdataset);
  if (boxdataspace < 0) return boxdataspace;

  H5Sget_simple_extent_dims(boxdataspace, dims, maxdims);

  memdataspace = H5Screate_simple(1, dims, NULL);

  *length = dims[0];
  *boxes = malloc(dims[0]*sizeof(box));
  if (handle.dim == 2)
    {
      error = H5Dread(boxdataset, box2d_id, memdataspace, boxdataspace,
              H5P_DEFAULT, *boxes);
    }
  else if (handle.dim == 3)
    {
      error = H5Dread(boxdataset, box3d_id, memdataspace, boxdataspace,
              H5P_DEFAULT, *boxes);
    }

  H5Dclose(boxdataset);
  H5Sclose(boxdataspace);
  H5Sclose(memdataspace);
  return 0;
}

int writeBoxes(box* boxes, int length, HDF5Handle handle)
{
  herr_t ret;
  hssize_t offset[1];
  hsize_t  flatdims[1], count[1];
  hid_t boxdataspace, boxdataset, memdataspace;
  int i;

  count[0] = 1;
  flatdims[0] = length;

  boxdataspace = H5Screate_simple(1, flatdims, NULL);
  memdataspace = H5Screate_simple(1, count, NULL);
  if (handle.dim == 2)
    {
      boxdataset   = H5Dcreate(handle.group_ID, "boxes",  box2d_id,
                     boxdataspace, H5P_DEFAULT);
      if (boxdataset < 0) return boxdataset;
      for (i=0; i<length; ++i)
    {
       offset[0] = i;
       ret = H5Sselect_hyperslab (boxdataspace, H5S_SELECT_SET, offset, NULL,
                      count, NULL);
       if (ret < 0) return ret;
       ret = H5Dwrite(boxdataset, box2d_id, memdataspace, boxdataspace,
              H5P_DEFAULT, boxes + i);
       if (ret < 0) return ret;
    }
    }
  else
    {
      boxdataset   = H5Dcreate(handle.group_ID, "boxes",  box3d_id,
                     boxdataspace, H5P_DEFAULT);
      if (boxdataset < 0) return boxdataset;
      for (i=0; i<length; ++i)
    {
      offset[0] = i;
      ret = H5Sselect_hyperslab (boxdataspace, H5S_SELECT_SET, offset, NULL,
                     count, NULL);
      if (ret < 0) return ret;
      ret = H5Dwrite(boxdataset, box3d_id, memdataspace, boxdataspace,
             H5P_DEFAULT, boxes + i);
      if (ret < 0) return ret;
    }
    }

  H5Sclose(boxdataspace);
  H5Sclose(memdataspace);
  H5Dclose(boxdataset);
  return 0;
}

int writeHDF5attributes(HDF5attributes* attrib, HDF5Handle handle)
{
  H5E_auto_t efunc; void* edata;
  herr_t  ret;
  hid_t aid, attr, s_type;
  Attribute* atr;
  Attribute* at;

  H5Eget_auto(&efunc, &edata);

#define INSERT(Ttype, attributePtr, H5Ttype)                              \
  at = attributePtr;                                                      \
  while (at != NULL)                                                       \
    {                                                                     \
      aid  = H5Screate(H5S_SCALAR);                                       \
      H5Eset_auto(NULL, NULL);                                            \
      attr = H5Acreate(handle.group_ID, at->name, H5Ttype,                \
                 aid, H5P_DEFAULT);                           \
      if (attr < 0)                                                        \
    {                                                                 \
      H5Adelete(handle.group_ID, at->name);                           \
      attr = H5Acreate(handle.group_ID, at->name, H5Ttype,            \
               aid, H5P_DEFAULT);                             \
      if (attr < 0)                                                    \
        {                                                             \
          return -1;                                      \
        }                                                             \
    }                                                                 \
      H5Eset_auto(efunc, edata);                                          \
      ret = H5Awrite(attr, H5Ttype, at->data);                            \
      if (ret < 0) return ret;                                             \
      H5Sclose(aid);                                                      \
      H5Aclose(attr);                                                     \
      at = (Attribute *)at->next;                                         \
    }                                                                     \

    atr = attrib->accessByType[INTEGER];
    INSERT(int,  atr, H5T_NATIVE_INT);

    atr = attrib->accessByType[FLOAT];
    INSERT(float,   atr,   H5T_NATIVE_FLOAT);

    atr = attrib->accessByType[DOUBLE];
    INSERT(double,  atr,  H5T_NATIVE_DOUBLE);

    atr = attrib->accessByType[INTVECT2D];
    INSERT(intvect2d, atr,  intvect2d_id);

    atr = attrib->accessByType[INTVECT3D];
    INSERT(intvect3d, atr,  intvect3d_id);

    atr = attrib->accessByType[BOX2D];
    INSERT(box2d, atr,  box2d_id);

    atr = attrib->accessByType[BOX3D];
    INSERT(box3d, atr,  box3d_id);

    /* string is different, of course */
    at = attrib->accessByType[CHAR];
     while (at != NULL)
    {
      aid    = H5Screate(H5S_SCALAR);
      s_type = H5Tcopy(H5T_C_S1);
      H5Tset_size(s_type, strlen((char*)at->data)); /*extra requirement for strings*/
      H5Eset_auto(NULL, NULL);
      attr = H5Acreate(handle.group_ID, at->name, s_type,
                 aid, H5P_DEFAULT);
      if (attr < 0)
    {
      H5Adelete(handle.group_ID, at->name);
      attr = H5Acreate(handle.group_ID, at->name, s_type,
               aid, H5P_DEFAULT);
      if (attr < 0)
        {
          return -1;
        }
    }
      H5Eset_auto(efunc, edata);
      ret = H5Awrite(attr, s_type, at->data);
      if (ret < 0) return ret;
      H5Sclose(aid);
      H5Aclose(attr);
      H5Tclose(s_type);
      at = (Attribute *)at->next;
    }

    return 0;
}

int writeEBChomboFile(const char* filename,
                      box3d domain,
                      int length,
                      box3d* boxes,
                      long*  regoffset,
                      long*  irregoffset,
                      long numreg,
                      long numirreg,
                      regvof* regularVofs,
                      irregvof* irregularVofs)
{
  return 0;
}

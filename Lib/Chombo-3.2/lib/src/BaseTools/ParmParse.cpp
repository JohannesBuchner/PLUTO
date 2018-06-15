#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*
  DISCLAIMER:
  This is ugly because it consists of four files from boxlib
  concatenated together and munged hopelessly.  This was done
  1) to greatly reduce the size of the API of boxlib
  2) to avoid the godawful task of rewriting parmparse
  3) to avoid cluttering the global namespace
 */

#include <iostream>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
using std::cin;
using std::cout;
using std::cerr;
using std::istream;
using std::ostream;

#include "MayDay.H"
#include "ParmParse.H"
#include "BaseNamespaceHeader.H"

int             ParmParse::xargc   = 0;
int             ParmParse::num_obj = 0;
char**          ParmParse::xargv   = 0;
PP_List<PP_entry*> ParmParse::table;

static
const char* const
tok_name[] =
{
   "NAME",
   "OPTION",
   "INTEGER",
   "UNSIGNEDLONG",
   "FLOAT",
   "DOUBLE",
   "STRING",
   "BOOL",
   "EQ_SIGN",
   "EOF"
};

PP_entry::PP_entry (PP_String&          name,
                    ParmParse::PPType typ,
                    PP_List<PP_String>&    vals)
    :
    defname(name),
    deftype(typ),
    val(vals.length())
{
   PP_ListIterator<PP_String> li(vals);
   for (int i = 0; li; i++, ++li)
      val[i] = vals[li];
}

ParmParse::ParmParse (int         argc,
                      char**      argv,
                      const char* prefix,
                      const char* parfile)
{
    define( argc, argv, prefix, parfile );
}

void
ParmParse::define (int         argc,
                   char**      argv,
                   const char* prefix,
                   const char* parfile)
{
    if (table.length() > 0)
       MayDay::Abort("ParmParse::define(): table already built");
    num_obj++;
    xargc = argc;
    xargv = argv;
    if (prefix != 0)
       thePrefix = prefix;
    ppinit(parfile);
}

ParmParse::ParmParse (const char* prefix)
{
  if (prefix != 0)
  {
    thePrefix = prefix;
  }
  // don't count this object if there isn't another one already --
  // just assume define() will be called later
  if ( num_obj > 0 ) num_obj++;
}

void ParmParse::dumpTable (ostream& os) const
{
   for (PP_ListIterator<PP_entry*> li(table); li; ++li)
      li()->dump(os);
}

//
// Initialize ParmParse.
//

void
ParmParse::ppinit (const char* parfile)
{
    if (parfile != 0)
       read_file(parfile,table);

    if (xargc > 0)
    {
        PP_String argstr;
        const char SPACE = ' ';
        for (int i = 0; i < xargc; i++)
        {
            argstr += xargv[i];
            argstr += SPACE;
        }
        PP_List<PP_entry*> arg_table;
        bldTable(argstr.c_str(), argstr.length()+1, arg_table);
        //
        // Append arg_table to end of existing table.
        //
        table.catenate(arg_table);
    }
}

ParmParse::~ParmParse ()
{
   if (--num_obj == 0)
   {
      for (PP_ListIterator<PP_entry*> li(table); li; ++li)
         delete table[li];
      table.clear();
   }
}

//
// Keyword aware string comparison.
//

static
bool
ppfound (const char*    keyword,
         const PP_String& key,
         const PP_String& prefix)
{
    //
    // Return true if key==keyword || key == prefix.keyword.
    //
    if (!prefix.isNull())
    {
        PP_String tmp(prefix);
        tmp += '.';
        tmp += keyword;
        return (key == tmp);
    }
    else
        return (key == keyword);
}

//
// Return number of occurences of parameter name.
//

int
ParmParse::countname (const char* name) const
{
    int cnt = 0;
    for (PP_ListIterator<PP_entry*> li(table); li; ++li)
       if (ppfound(name,li()->defname,thePrefix))
           cnt++;
    return cnt;
}

//
// Return true if name in table.
//

bool
ParmParse::contains (const char* name) const
{
    for (PP_ListIterator<PP_entry*> li(table); li; ++li)
    {
       if (ppfound(name,li()->defname,thePrefix))
           return true;
    }
    return false;
}

bool
ParmParse::contains (const std::string& name) const
{
    for (PP_ListIterator<PP_entry*> li(table); li; ++li)
    {
       if (ppfound(name.c_str(), li()->defname, thePrefix))
           return true;
    }
    return false;
}

//
// Return the index of the n'th occurence of a parameter name,
// except if n==-1, return the index of the last occurence.
// Return 0 if the specified occurence does not exist.
//

const PP_entry*
ParmParse::ppindex (int         n,
                    const char* name) const
{
    const PP_entry* fnd = 0;
    if (n < 0)
    {
        //
        // Search from back of list.
        //
        for (PP_ListIterator<PP_entry*> li = table.last(); li; --li)
        {
            if (ppfound(name,li()->defname,thePrefix))
            {
                fnd = li();
                break;
            }
        }
    }
    else
    {
        for (PP_ListIterator<PP_entry*> li(table); li; ++li)
        {
            if (ppfound(name,li()->defname,thePrefix))
            {
                fnd = li();
                n--;
                if (n < 0)
                    break;
            }
        }
        if (n >= 0)
            fnd = 0;
    }
    return fnd;
}

void
ParmParse::getval (const char*  name,
                   const PPType type,
                   void*        ptr,
                   int          ival,
                   int          occurence) const
{
    if (queryval(name,type,ptr,ival,occurence) == 0)
    {
        cerr << "ParmParse::getval ";
        if (occurence >= 0)
            cerr << "occurence number "
                 << occurence
                 << " of ";

        if (!thePrefix.isNull())
            cerr << thePrefix << '.';

        cerr << "ParmParse::getval(): "
             << name
             << " not found in table"
             << '\n';
        dumpTable(cerr);
        MayDay::Abort();
    }
}

int
ParmParse::queryval (const char*  name,
                     const PPType type,
                     void*        ptr,
                     int          ival,
                     int          occurence) const
{
    //
    // Get last occurrance of name in table.
    //
    const PP_entry *def = ppindex(occurence,name);
    if (def == 0)
        return 0;
    //
    // Does it have ival values?
    //
    if (ival >= def->val.length())
    {
        cerr << "ParmParse::queryval no value number"
             << ival << " for ";
        if (occurence < 0)
            cerr << "last occurence of ";
        else
            cerr << " occurence " << occurence << " of ";
        cerr << def->defname << '\n';
        def->dump(cerr);
        MayDay::Abort();
    }

    const PP_String& valname = def->val[ival];

    int ok;
    double val_dbl;
    //
    // Retrieve value.
    //
    switch (type)
    {
    case ppInt:
        ok = isInteger(valname,*(int*)ptr);
        break;
    case ppFloat:
        ok = isDouble(valname,val_dbl);
        if (ok)
            *(float*)ptr = (float) val_dbl;
        break;
    case ppDouble:
        ok = isDouble(valname,val_dbl);
        *(double*)ptr = val_dbl;
        break;
    case ppString:
        ok = true;
        *(PP_String*)ptr = valname;
        break;
    case ppBool:
        ok = isBool(valname,*(bool*)ptr) ;
        break ;
    case ppUnsignedLong:
        ok = isUnsigned(valname,*(unsigned long*)ptr) ;
        break ;
    default:
        ok = false;
        MayDay::Abort("ParmParse::queryval invalid type");
    }
    if (!ok)
    {
        cerr << "ParmParse::queryval type mismatch on value number "
             << ival << " of " << '\n';
        if (occurence < 0)
            cerr << " last occurence of ";
        else
            cerr << " occurence number " << occurence << " of ";
        cerr << def->defname << '\n';
        cerr << " Expected "
             << tok_name[type]
             << " value = "
             << valname << '\n';
        def->dump(cerr);
        MayDay::Abort();
    }
    return 1;
}

void
ParmParse::getarr (const char*  name,
                   const PPType type,
                   void*        ptr,
                   int          start_ix,
                   int          num_val,
                   int          occurence) const
{
    if (queryarr(name,type,ptr,start_ix,num_val,occurence) == 0)
    {
        cerr << "ParmParse::getarr ";
        if (occurence >= 0)
            cerr << "occurence number " << occurence << " of ";
        if (!thePrefix.isNull())
            cerr << thePrefix << '.';
        cerr << "ParmParse::getarr(): "
             << name
             << " not found in table"
             << '\n';
        dumpTable(cerr);
        MayDay::Abort();
    }
}

int
ParmParse::queryarr (const char*  name,
                     const PPType type,
                     void*        ptr,
                     int          start_ix,
                     int          num_val,
                     int          occurence) const
{
    //
    // Get last occurrance of name in table.
    //
    const PP_entry *def = ppindex(occurence,name);
    if (def == 0)
        return 0;
    //
    // Does it have sufficient number of values and are they all
    // the same type?
    //
    int stop_ix = start_ix + num_val - 1;
    if (stop_ix >= def->val.length())
    {
        cerr << "ParmParse::queryarr too many values requested for";
        if (occurence < 0)
            cerr << " last occurence of ";
        else
            cerr << " occurence " << occurence << " of ";
        cerr << def->defname << '\n';
        def->dump(cerr);
        MayDay::Abort();
    }

    for (int n = start_ix; n <= stop_ix; n++)
    {
       const PP_String& valname = def->val[n];
       //
       // Retrieve value.
       //
       int ok = false;
       double val_dbl;
       switch (type)
       {
       case ppInt:
           ok = isInteger(valname,*(int*)ptr);
           ptr = (int*)ptr+1;
           break;
       case ppFloat:
           ok = isDouble(valname,val_dbl);
           if (ok)
               *(float*)ptr = (float) val_dbl;
           ptr = (float*)ptr+1;
           break;
       case ppDouble:
           ok = isDouble(valname,*(double*)ptr);
           ptr = (double*)ptr+1;
           break;
       case ppString:
           ok = true;
           *(PP_String*)ptr = valname;
           ptr = (PP_String*)ptr+1;
           break;
       case ppBool:
           ok = isBool(valname,*(bool*)ptr) ;
           ptr = (bool*)ptr+1;
           break ;
       case ppUnsignedLong:
           ok = isUnsigned(valname,*(unsigned long*)ptr);
           ptr = (unsigned long*)ptr+1;
           break;
       default:
           MayDay::Abort("ParmParse::get invalid type");
       }
       if (!ok)
       {
           cerr << "ParmParse::queryarr type mismatch on value number "
                <<  n << " of ";
           if (occurence < 0)
               cerr << " last occurence of ";
           else
               cerr << " occurence number " << occurence << " of ";
           cerr << def->defname << '\n';
           cerr << " Expected "
                << tok_name[type]
                << " value = "
                << valname << '\n';
           def->dump(cerr);
           MayDay::Abort();
       }
    }

    return 1;
}

void
ParmParse::read_file (const char*      fname,
                      PP_List<PP_entry*>& tab)
{
    //
    // Space for input file if it exists.
    // Note: on CRAY, access requires (char*) not (const char*).
    //
    if (fname != 0 && fname[0] != 0)
    {
        FILE* pffd = fopen(fname, "rb");
        if (pffd == 0)
        {
            cerr << "ParmParse::read_file(): couldn't open \""
                 << fname
                 << "\"";
            MayDay::Abort();
        }
        //
        // Get the length.
        //
        fseek(pffd, 0, 2);
        int pflen = (int)ftell(pffd);
        rewind(pffd);
        char* str = new char[pflen+1];
        memset(str,0,pflen+1);
        int nread = fread(str, 1, pflen, pffd);
        if (!(nread == pflen))
        {
            cerr << "ParmParse::read_file(): fread() only "
                 << nread
                 << " bytes out of "
                 << pflen
                 << " from "
                 << fname;
            MayDay::Abort();
        }
        fclose(pffd);
        bldTable(str,pflen+1,tab);
        delete [] str;
    }
}

static
void
eat_garbage (const char* str,
             int&        i,
             int         len)
{
    for (;;)
    {
        if (i < len && str[i] == '#')
        {
            while (i < len && str[i] != '\n')
                i++;
        }
        else if (i < len && isspace(str[i]))
            i++;
        else
            break;
    }
}

//
// Simple lexical analyser.
//

enum lexState
{
    START,
    MINUS,
    SIGN,
    OPTION,
    STRING,
    QUOTED_STRING,
    INTEGER,
    START_FRACTION,
    FRACTION,
    START_EXP,
    SIGNED_EXP,
    EXP,
    PREFIX,
    SUFFIX,
    STAR
};

static
const char* const
state_name[] =
{
   "START",
   "MINUS",
   "SIGN",
   "OPTION",
   "STRING",
   "QUOTED_STRING",
   "INTEGER",
   "START_FRACTION",
   "FRACTION",
   "START_EXP",
   "SIGNED_EXP",
   "EXP",
   "PREFIX",
   "SUFFIX",
   "STAR"
};

ParmParse::PPType
ParmParse::getToken (const char* str,
                     int&        i,
                     int         slen,
                     char*       ostr)
{
#define ERROR_MESSAGE \
   ostr[k++] = '\0'; \
   cerr << "ParmParse::getToken(): invalid string = " << ostr << '\n'; \
   cerr << "STATE = " << state_name[state] \
        << ", next char = " << ch << '\n'; \
   cerr << ", rest of input = \n" << (str+i) << '\n'; \
   MayDay::Abort()
   //
   // Eat white space and comments.
   //
   eat_garbage(str,i,slen);
   //
   // Check for end of file.
   //
   if (i >= slen || str[i] == '\0')
       return ppEOF;
   //
   // Start token scan.
   //
   lexState state = START;
   int k = 0;
   while (1)
   {
       if (i >= slen)
       {
         MayDay::Abort("ParmParse::getToken - walked off the end of the string");
       }
       char ch = str[i];
       switch (state)
       {
       case START:
           if (ch == '=')
           {
               ostr[k++] = ch; i++;
               ostr[k++] = 0;
               return ppEQ_sign;
           }
           else if (ch == '"')
           {
               i++;
               state = QUOTED_STRING;
           }
           else if (ch == '*')
           {
               ostr[k++] = ch; i++;
               state = STAR;
           }
           else if (isalpha(ch) || ch == '_')
           {
               ostr[k++] = ch; i++;
               state = PREFIX;
           }
           else if (ch == '-')
           {
               ostr[k++] = ch; i++;
               state = MINUS;
           }
           else
           {
               ostr[k++] = ch; i++;
               state = STRING;
           }
           break;
       case MINUS:
           if (isalpha(ch) || ch == '_')
           {
               k--;
               ostr[k++] = ch; i++;
               state = OPTION;
           }
           else
           {
               ostr[k++] = ch; i++;
               state = STRING;
           }
           break;
       case OPTION:
           if (isalnum(ch) || ch == '_' || ch == '.')
           {
               ostr[k++] = ch; i++;
           }
           else if (isspace(ch) || ch == '\0')
           {
               ostr[k++] = 0;
               return ppOption;
           }
           else
           {
               ERROR_MESSAGE;
           }
           break;
       case STAR:
           if (ch == '.')
           {
               ostr[k++] = ch; i++;
               state = SUFFIX;
           }
           else
           {
               ostr[k++] = ch; i++;
               state = STRING;
           }
           break;
       case PREFIX:
           if (isalnum(ch) || ch == '_')
           {
               ostr[k++] = ch; i++;
           }
           else if (ch == '.')
           {
               ostr[k++] = ch; i++;
               state = SUFFIX;
           }
           else if (isspace(ch) || ch == '\0' || ch == '=')
           {
               ostr[k++] = 0;
               return ppDefn;
           }
           else
           {
               ostr[k++] = ch; i++;
               state = STRING;
           }
           break;
       case SUFFIX:
           if (isalnum(ch) || ch == '_')
           {
               ostr[k++] = ch; i++;
           }
           else if (isspace(ch) || ch == '\0' || ch == '=')
           {
               ostr[k++] = 0;
               return ppDefn;
           }
           else
           {
               ostr[k++] = ch; i++;
               state = STRING;
           }
           break;
       case QUOTED_STRING:
           if (ch == '"')
           {
               i++;
               ostr[k++] = 0;
               return ppString;
           }
           else
           {
               ostr[k++] = ch; i++;
           }
           break;
       case STRING:
           if (isspace(ch) || ch == '\0' || ch == '=')
           {
               ostr[k++] = 0;
               return ppString;
           }
           else
           {
               ostr[k++] = ch; i++;
           }
           break;
       default:
           ERROR_MESSAGE;
       }
   }
#undef ERROR_MESSAGE
}

void
ParmParse::bldTable (const char*      str,
                     int              lenstr,
                     PP_List<PP_entry*>& tab)
{
   PP_String       cur_name;
   PP_List<PP_String> cur_list;
   PP_String       tmp_str;
   PP_entry      *pp;

   int       i = 0;
   PPType    token;
   const int SCRATCH_STR_LEN  = 200;
   char      tokname[SCRATCH_STR_LEN];

   while (true)
   {
      token = getToken(str,i,lenstr,tokname);

      switch (token)
      {
      case ppEOF:
          addDefn(cur_name,cur_list,tab);
          return;
      case ppOption:
          addDefn(cur_name,cur_list,tab);
          tmp_str = tokname;
          pp = new PP_entry(tmp_str,ppOption,cur_list);
          tab.append(pp);
          break;
      case ppEQ_sign:
          if (cur_name.length() == 0)
              MayDay::Abort("ParmParse::bldTable() EQ with no current defn");
          if (cur_list.isEmpty())
              //
              // First time we see equal sign, just ignore.
              //
              break;
          //
          // Read one too far, remove last name on list.
          //
          tmp_str = cur_list.lastElement();
          cur_list.removeLast();
          addDefn(cur_name,cur_list,tab);
          cur_name = tmp_str;
          break;
      case ppDefn:
          if (cur_name.length() == 0)
          {
              cur_name = tokname;
              break;
          }
          //
          // Otherwise, fall through, this may be a string.
          //
      case ppInt:
      case ppUnsignedLong:
      case ppFloat:
      case ppDouble:
      case ppString:
      case ppBool:
          if (cur_name.length() == 0)
              MayDay::Abort("ParmParse::bldTable(): value with no defn");
          cur_list.append(tokname);
          break;
      }
   }
}

void
ParmParse::addDefn (PP_String&         def,
                    PP_List<PP_String>&   val,
                    PP_List<PP_entry*>& tab)
{
    static const PP_String FileKeyword("FILE");
    //
    // Check that defn exists.
    //
    if (def.length() == 0)
    {
        val.clear();
        return;
    }
    //
    // Check that it has values.
    //
    if (val.isEmpty())
    {
        cerr << "ParmParse::addDefn(): no values for definition " << def;
        MayDay::Abort();
    }
    //
    // Check if this defn is a file include directive.
    //
    if (def == FileKeyword && val.length() == 1)
    {
        //
        // Read file and add to this table.
        //
        const char* fname = val.firstElement().c_str();
        read_file(fname, tab);
    }
    else
    {
      PP_entry* pp = new PP_entry(def,ppDefn,val);
      tab.append(pp);

    }
    val.clear();
    def = PP_String();
}

void
PP_entry::dump (ostream& os) const
{
    static const char TokenInitial[] =
    {
      'N','O','I','F','D','S','=','E'
    };

    char tmp[1024];
    long nval = val.length();
    sprintf(tmp,
            "(%c,%1d) %15s :: ",
            TokenInitial[deftype],
            int(nval),
            defname.c_str());
    os << tmp;
    for (int i = 0; i < nval; i++)
    {
       os << " ("
          << TokenInitial[ParmParse::ppString]
          << ','
          << val[i]
          << ')';
    }
    os << '\n';

    if (os.fail())
        MayDay::Abort("PP_entry::dump(ostream&) failed");
}

void
PP_StringRep::resize (int n)
{
    if (n > bufferlength)
    {
        char* ns = new char[n];
        ::memcpy(ns,s,bufferlength);
        bufferlength = n;
        delete [] s;
        s = ns;
    }
}

PP_String::PP_String ()
    : p(new PP_StringRep(1)),
      len(0)
{
    p->s[0] = 0;
}

PP_String::PP_String (char c)
    : p(new PP_StringRep(2)),
      len(1)
{
    p->s[0] = c;
    p->s[1] = 0;
    if (c == '\0')
        len = 0;
}

PP_String::PP_String (int size)
    : p(new PP_StringRep(size+1)),
      len(0)
{
    ::memset(p->s,'\0',p->bufferlength);
}

PP_String::PP_String (const char* initialtext)
{
    CH_assert(initialtext != 0);
    len = ::strlen(initialtext);
    p = new PP_StringRep(len + 1);
    ::memcpy(p->s,initialtext,len+1);
}

PP_String::PP_String (const PP_String& initialstring)
    : p(initialstring.p),
      len(initialstring.len)
{
}

PP_String&
PP_String::operator= (const PP_String& rhs)
{
    p   = rhs.p;
    len = rhs.len;
    return *this;
}

PP_String&
PP_String::operator+= (const PP_String& val)
{
    copyModify();
    int clen = length() + val.length();
    p->resize(clen+1);
    ::memcpy(&(p->s[len]),val.p->s, val.length()+1);
    len = clen;
    return *this;
}

PP_String&
PP_String::operator+= (const char* s)
{
    CH_assert(s != 0);
    copyModify();
    int slen = ::strlen(s);
    int clen = length() + slen;
    p->resize(clen+1);
    ::memcpy(&(p->s[len]),s, slen+1);
    len = clen;
    return *this;
}

PP_String&
PP_String::operator+= (char c)
{
    if (!(c == '\0'))
    {
        copyModify();
        p->resize(len+2);
        p->s[len++] = c;
        p->s[len]   = 0;
    }
    return *this;
}

char&
PP_String::operator[] (int index)
{
    CH_assert(index >= 0 && index < len);
    copyModify();
    return p->s[index];
}

void
PP_String::copyModify ()
{
    if (!p.unique())
    {
        PP_StringRep* np = new PP_StringRep(len+1);
        ::memcpy(np->s,p->s,len+1);
        p = np;
    }
}

PP_String&
PP_String::toUpper ()
{
    copyModify();
    for (char *pp = p->s; *pp != 0; pp++)
        *pp = toupper(*pp);
    return *this;
}

PP_String&
PP_String::toLower ()
{
    copyModify();
    for (char *pp = p->s; *pp != 0; pp++)
        *pp = tolower(*pp);
    return *this;
}

istream&
operator>> (istream& is,
            PP_String& str)
{
    const int BufferSize = 128;
    char buf[BufferSize + 1];
    int index = 0;
    //
    // Nullify str.
    //
    str = "";
    //
    // Eat leading whitespace.
    //
    char c;
    do
    {
      is.get(c);
    } while (is.good() && isspace(c));
    buf[index++] = c;
    //
    // Read until next whitespace.
    //
    while (is.get(c) && !isspace(c))
    {
        buf[index++] = c;
        if (index == BufferSize)
        {
            buf[BufferSize] = 0;
            str += buf;
            index = 0;
        }
    }
    is.putback(c);
    buf[index] = 0;
    str += buf;
    if (is.fail())
        MayDay::Abort("operator>>(istream&,PP_String&) failed");
    return is;
}

ostream&
operator<< (ostream&       out,
            const PP_String& str)
{
    out.write(str.c_str(), str.len);
    if (out.fail())
        MayDay::Abort("operator<<(ostream&,PP_String&) failed");
    return out;
}

istream&
PP_String::getline (istream& is)
{
    char      c;
    const int BufferSize = 100000;
    char      buf[BufferSize + 1];
    int       index = 0;

    *this = "";
    //
    // Get those characters.
    // We read the newline but don't add it to the string.
    //
    while (is.get(c))
    {
        if (c == '\n')
            break;

        buf[index++] = c;

        if (index == BufferSize)
        {
            buf[BufferSize] = 0;
            *this += buf;
            index = 0;
        }
    }
    buf[index] = 0;
    *this += buf;

    if (!(is || is.eof()))
        MayDay::Abort("PP_String::getline(istream&) failed");

    return is;
}
#include "BaseNamespaceFooter.H"

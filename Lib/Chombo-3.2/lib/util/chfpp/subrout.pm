#!/usr/local/gnu/bin/perl -w
####
####    _______              __
####   / ___/ /  ___  __ _  / /  ___
####  / /__/ _ \/ _ \/  V \/ _ \/ _ \
####  \___/_//_/\___/_/_/_/_.__/\___/
####  Please refer to Copyright.txt, in Chombo's root directory.
####

# -------------------------------------------------------------------------
# this
# -------------------------------------------------------------------------

use strict;

require "Balanced.pm";
Text::Balance->import();

use ChF ;

package SubroutProc;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
$VERSION = 1.60;
@ISA = qw(Exporter);
@EXPORT = qw(&procSubrout);
@EXPORT_OK = qw(&procSubrout);
%EXPORT_TAGS = ();

sub SubroutProc::procSubrout
{
    use strict;
### Constants
    my ($FINFile, $COUTFile, $FOUTFile, $SpaceDim, $debug, $basename) = @_;
    $SubroutProc::indentstr = "\n      ";
    $SubroutProc::continstr = "\n     &           ";
    $SubroutProc::dimension = $SpaceDim;
    %SubroutProc::UsesCHFID = () ;  # for handling CHF_ID() array

    if($debug)
    {
        print "SubroutProc: \n";
        print "input file  = $FINFile \n";
        print "fortran output file = $FOUTFile \n";
        print "c prototype file = $COUTFile \n";
        print "SpaceDim    = $SpaceDim \n";
        print "basename    = $basename \n";
    }


    open(SubroutProc::COUT,">" . $COUTFile)
        or die "Error: cannot open output file " . $COUTFile . "\n";
    open(SubroutProc::FOUT,">" . $FOUTFile)
        or die "Error: cannot open output file " . $FOUTFile . "\n";
    open(SubroutProc::FIN, "<" . $FINFile)
        or die "Error: cannot open input file " . $FINFile . "\n";

### initialize the output files
    &startFFile();
    &startCFile($basename);


### read the input file and determine which subroutines
### use the CHF_ID() array; set the %UsesCHFID hash and
### reset FIN to the start of the file

    &findCHFID( $debug ) ;

### loop through lines in the input file and stuff them
### into one buffer
    my $callbuf = "" ;
    my $insubst = 0;
    my $incall  = 0;
    my $infirstline = 0;
    my $numopenparens = 0;
    my $linebuf = <FIN> ;
    my $nextbuf;
    my $endstatement = 0 ;
    while (defined $linebuf)
    {
        # look ahead to the next line to find statement ends
        if( $nextbuf = <FIN> )
        {
            # any char except "0" in column6 is a continuation line
            # or any line with "&" as the first non-whitespace
            if( $nextbuf !~ m/^     [^0 ]/ and $nextbuf !~ m/^\s*&/ ){
                $endstatement = 1;
            }else{
                $endstatement = 0;
            }
        }else{
            $endstatement = 1;
        }
        ###if we find a subroutine or a subroutine call
        ###go into buffering mode
        if ($linebuf =~ s/\! *no_chf//i)
        {
            # unless ! no_chf is detected in which case just strip the comment
            if($debug)
            {
                print "subroutProc::no ChF processing of line = $linebuf\n";
            }
        }
        elsif ($linebuf =~ m/^\s*(subroutine|chf_call)\s*/i)
        {
            if( $1 =~ m/chf_call/i ) { $incall  = 1; }
            # print everything up to the name of the subr
            if ($incall == 1)
            {
                print SubroutProc::FOUT "      call ";
                $linebuf =~ s/chf_call//i;
            }
            else
            {
                print SubroutProc::FOUT $&;
                $linebuf = $';  #'
            }

            if($debug)
            {
                print "subroutProc::first line = $linebuf\n";
            }
            # and remove it from the buffer

            # and set flags to remember what kind of statement this is
            $insubst = 1;
            $infirstline = 1;

        }
        elsif ($linebuf =~ m/\)\s*chf_call\s*/i)
        {
            # This is the tricky part - CALL can appear
            # after IF() but it could also be a variable
            # named call<something>.  To be rigorous, you
            # have to parse to the end of the parenthesized
            # expression and see if there is an "=".
            # But that's too hard so we're going to cheat
            # and assume it isn't a misleading variable name.

            # (DFM -- 10/29/10) this is a bit silly -- change chf_call->call
            # and then redo "if" statment to make the $& variable be the correct 
            # value. Hopefully somebody who knows perl can come back later and fix this
            # in some sort of reasonable way.

            $linebuf =~ s/chf_call/call/i;
            if ($linebuf =~ m/\)\s*call\s*/i)
            {
                # print everything up to the name of the subr     
                print SubroutProc::FOUT $`.$&;
                
                # and remove it from linebuf
                $linebuf = $' ;  #'
                    # and tell the arg list processor this is a CALL
                    # but not at the start of line (so dont set $insubst)
                    $incall = 1;
                # but fake being at the start of the statement
                # so the first 6 chars aren't stripped
                $infirstline = 1;
                if($debug)
                {
                    print "subroutProc::first call line = $linebuf\n";
                }
            }
        }
        ###if we are in buffering mode,
        ###add  the line to the buffer
        if($insubst||$incall)
        {
            ####this strips out beginning spaces and
            ####continuation characters
            ### only do this on non-first lines
            my $nosixbuf = "";
            if($infirstline)
            {
                $nosixbuf = $linebuf;
                $infirstline = 0;
            }
            else
            {
                $nosixbuf = substr($linebuf, 6);
            }
            $callbuf .= $nosixbuf;
            # has the whole statement been read?
            if($endstatement == 1)
            {
                my $stripbuf = $callbuf;
                ### strip out line breaks and white space
                $stripbuf =~ s/[\n\s]//sg;
                if($debug)
                {
                    print "subroutProc::whole buffer = $callbuf\n";
                    my $striplen = length($stripbuf);
                    print "subroutProc::length stripped buffer = $striplen \n";
                    print "subroutProc::stripped buffer = $stripbuf\n";
                }
                &procSubStatement($stripbuf,  $debug, $incall);
                ###reset in subrout flag to false
                ###and calling buffer to null
                $insubst = 0;
                $incall  = 0;
                $callbuf = "";
            }
        }
        else
        {
            ### not inside subroutine declaration
            ###just print it back out
            print SubroutProc::FOUT $linebuf;
        }
        $linebuf = $nextbuf ;
    }
    &finishCFile() ;

    ###close files.
    close(SubroutProc::COUT);
    close(SubroutProc::FIN);
    close(SubroutProc::FOUT);
    return 1;
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
sub SubroutProc::procSubStatement
{
    my ($fullstring,  $debug, $incall ) = @_;
    ###coming into this routine the string is everthing
    ###between subroutine and ")" possibly including ")"
    ###without any continuation characters or white space
    ###first get the name of the subroutine name
    ###which is the string before "("
    my $tempbuf = $fullstring;
    $tempbuf =~ m/\s*\(/g;
    my $offset = pos $tempbuf ;
    if( ! $offset )
    {
        # for a subroutine call with no parens
        if( $incall )
        {
            if( $debug ){ print "procSubStatement: subr call with no arglist [$fullstring]"; }
            $tempbuf .= "()";
            $tempbuf =~ m/\s*\(/g;
            $offset = pos $tempbuf ;
        }
        else
        {
            die "no beginning \( in sub string $fullstring\n";
        }
    }
    my $length = $offset - 1;
    my $subname = substr($fullstring, 0, $length);
    my $argstring = substr($fullstring, $offset-1);
    $argstring = SubroutProc::getList( $argstring );
    $fullstring =~ quotemeta $argstring ;
    my $remainder = $'; #'
    ###remove the now unnessesary parentheses
    $argstring =~ s/^\s*\(//ig;
    $argstring =~ s/\)\s*$//ig;
    if($debug)
    {
        print "procSubStatement:fullstring = $fullstring\n";
        print "procSubStatement:offset = $offset\n";
        print "procSubStatement:subname = $subname\n";
        print "procSubStatement:argstring = $argstring\n";
        if( $remainder ){ print "procSubStatement:remainder = $remainder\n"; }
    }
    &SubroutProc::doFortranProc($subname, $argstring, $debug, $incall ,$remainder);
    
    # if calling a subroutine, dont make a prototype
    if(!$incall) { &SubroutProc::doTimedCPrototype($subname, $argstring, $debug); }
}


# -------------------------------------------------------------------
# -------------------------------------------------------------------

sub SubroutProc::doFortranProc
{
    my ($subname, $argstring,  $debug, $incall, $remainder ) = @_;
###    print "in subroutProc\n";
###subname is the name of the subroutine
###argstring is the string of arguments
### without carraige returns, continuation characters
### or whitespace or parentheses
###incall is 1 if this is a subroutine call,
### else in a subroutine declaration
###remainder is the rest of the line to output after the arglist
###the number of indentations shall be eight
###and eight shall be the number of indents
    my $indentstr = $SubroutProc::indentstr;
    my $continstr = $SubroutProc::continstr;
    my $dimension = $SubroutProc::dimension;

    my %singtype = ( R => "REAL_T" , I => "integer" , C => "COMPLEX_T" );

    my $subnamecap = uc($subname);
    my @subargs = split(",",$argstring);
    my $arguments = "";
    my $moduledecs = "";
    my $declarations = "";
    my ( $smdim,$boxdec,$boxarg,$compnm,$fabint,@lboundarg,@uboundarg );

###create argument list and declaration stuff
    for(my $iarg = 0; $iarg <= $#subargs; $iarg++)
    {
        my $singarg = $subargs[$iarg];
        if($debug)
        {
            print "doFortranProc: singarg = $singarg\n";
        }
        ###the first argument of any declaration will
        ###only have a comma if it is not the first argument
        ###overall
        my $comma = "";
        if($iarg > 0)
        {
            $comma = ",";
        }
        # handle special cases of CHF macros that can appear in an arglist
        # but shouldn't be parsed here
        if($singarg =~ /CHF\_(IX|[LU]BOUND)/ig)
        {
            $arguments .= $continstr.$comma.$singarg;
        }
        else
        {
            my $singname = &SubroutProc::getNamedThing($singarg);
            my $chf_decl = 1 ; #assume this is a CHF_* declaration
            if($debug)
            {
                print "doFortranProc: singname= $singname\n";
            }
            if($singarg =~ /CHF\_USE/ig)
            {
                $moduledecs .= $indentstr."USE $singname";
            }
            elsif($singarg =~ /CHF(\_CONST)?\_INTVECT/ig)
            {
                #[NOTE: *INTVECT commands must come before *INT to match properly.]
                $smdim = $dimension -1;
                $arguments .= $continstr.$comma.$singname;
                $declarations .= $indentstr."integer $singname(0:".$smdim.")";
            }
            elsif($singarg =~ /CHF(\_CONST)?\_REALVECT/ig)
            {
                #[NOTE: *REALVECT commands must come before *REAL to match properly.]
                $smdim = $dimension -1;
                $arguments .= $continstr.$comma.$singname;
                $declarations .= $indentstr."REAL_T $singname(0:".$smdim.")";
            }
            elsif($singarg =~ /CHF(\_CONST)?\_INT/ig)
            {
                $arguments .= $continstr.$comma.$singname;
                $declarations .= $indentstr."integer $singname";
            }
            elsif($singarg =~ /CHF(\_CONST)?\_REAL/ig)
            {
                $arguments .= $continstr.$comma.$singname;
                $declarations .= $indentstr."REAL_T $singname";
            }
            elsif($singarg =~ /CHF(\_CONST)?\_COMPLEX/ig)
            {
                $arguments .= $continstr.$comma.$singname;
                $declarations .= $indentstr."COMPLEX_T $singname";
            }
            elsif($singarg =~ /CHF\_BOX/ig)
            {
                ### box declaration has all the neccessary carraige returns
                ### and continuation characters.
                $boxdec = &SubroutProc::getBoxDecl($singname);
                $boxarg = &SubroutProc::getBoxArgs($singname);
                $declarations .= $boxdec;
                $arguments    .= $continstr.$comma.$boxarg;
            }
            # this handles {CHF,CHF_CONST}_F{I,R,C}A[1] because they're all
            # the same except for the type and the numcomps variable
            elsif($singarg =~ /CHF(\_CONST)?\_F([RIC])A(1?)/ig)
            {
                ### single or multiple component real|integer|complex fab.
                ### box declaration has all the neccessary carraige returns
                ### and continuation characters.
                ### the fab ints is just the ddecl part of the declaration
                ### because we need to do real, integer and complex fabs and
                ### both single and multiple-component fabs
                $compnm = "n".$singname."comp";
                $boxdec = &SubroutProc::getBoxDecl($singname);
                $fabint = &SubroutProc::getFabIntArgs($singname);
                $boxarg = &SubroutProc::getBoxArgs($singname);

                ### arguments are the name, the boxints and maybe ncomp
                $arguments .= $continstr.$comma.$singname;
                $arguments .= $continstr.",".$boxarg;
                # the num_comps variable is not used for the *A1 commands
                if( $3 eq "" ){
                    $arguments    .= $continstr.",".$compnm;
                    $declarations .= $indentstr."integer $compnm";
                }
                ###declare the boxints and the fab
                $declarations .= $boxdec;
                if($debug){ print "doFortranProc: CHF*F*A[1]: |".defined($1)?$1:""."|$2|".defined($3)?$3:""."|\n" ;}
                $declarations .= $indentstr.$singtype{uc($2)}." $singname\(";
                $declarations .= $continstr.$fabint;
                if( $3 eq "" ){$declarations .= ",".$continstr."0:$compnm-1";}
                $declarations .= "\)";
            }
            # this handles {CHF,CHF_CONST}_{VI,VR,I1D,R1D} because they're all the same
            # except for the type
            elsif($singarg =~ /CHF\_(?:CONST\_)?(V[IRC]|[IRC]1D)/ig)
            {
                #NOTE: ChF code can use CHF_UBOUND[] on this, but not CHF_LBOUND[]
                $1 =~ /([IRC])/i;
                $uboundarg[0]  = "i".$singname."hi0";
                $arguments    .= $continstr.$comma.$singname;
                $arguments    .= $continstr.",".$uboundarg[0];
                $declarations .= $indentstr."integer $uboundarg[0]";
                $declarations .= $indentstr.$singtype{uc($1)}." $singname\(";
                $declarations .= $continstr."0:$uboundarg[0]\)";
            }
            # this handles {CHF,CHF_CONST}_{I,R,C}[CH]ARRAY
            elsif($singarg =~ /CHF(?:\_CONST)?\_([IRC])(?:CH)?ARRAY/ig)
            {
                my $artype = $1;
                $singname =~ / *(RANK_SPACEDIM_PLUS_)?(\d) */i;
                my $arrank = $2+0;
                if (defined($1) && $1 eq "RANK_SPACEDIM_PLUS_")
                {
                    $arrank += $dimension;
                }
                ++$iarg;
                $subargs[$iarg] =~ / *(.*) *]/i;
                $singname = $1;
                $arguments .= $continstr.$comma.$singname;
                for (my $i = 0; $i != $arrank; ++$i)
                {
                    $lboundarg[$i] = "i".$singname."lo".$i;
                    $arguments .= $continstr.",".$lboundarg[$i];
                    $declarations .= $indentstr."integer $lboundarg[$i]";
                }
                for (my $i = 0; $i != $arrank; ++$i)
                {
                    $uboundarg[$i] = "i".$singname."hi".$i;
                    $arguments .= $continstr.",".$uboundarg[$i];
                    $declarations .= $indentstr."integer $uboundarg[$i]";
                }
                $declarations .= $indentstr.$singtype{uc($artype)}." $singname\(";
                for (my $i = 0; $i != $arrank; ++$i)
                {
                    my $endc = ",";
                    if ($i == $arrank-1)
                    {
                        $endc = ")";
                    }
                    $declarations .= $continstr."$lboundarg[$i]:$uboundarg[$i]".$endc;
                }
            }
            # this handles {CHF,CHF_CONST}_{VECTOR,MATRIX}
            elsif($singarg =~ /CHF(?:\_CONST)?\_(VECTOR|MATRIX)/ig)
            {
                my $arrank = ($1 eq "VECTOR") ? 1 : 2;
                $arguments .= $continstr.$comma.$singname;
                for (my $i = 0; $i != $arrank; ++$i)
                {
                    $lboundarg[$i] = "i".$singname."lo".$i;
                    $arguments .= $continstr.",".$lboundarg[$i];
                    $declarations .= $indentstr."integer $lboundarg[$i]";
                }
                for (my $i = 0; $i != $arrank; ++$i)
                {
                    $uboundarg[$i] = "i".$singname."hi".$i;
                    $arguments .= $continstr.",".$uboundarg[$i];
                    $declarations .= $indentstr."integer $uboundarg[$i]";
                }
                $declarations .= $indentstr.$singtype{"R"}." $singname\(";
                for (my $i = 0; $i != $arrank; ++$i)
                {
                    my $endc = ",";
                    if ($i == $arrank-1)
                    {
                        $endc = ")";
                    }
                    $declarations .= $continstr."$lboundarg[$i]:$uboundarg[$i]".$endc;
                }
            }
            ## handle ordinary arguments
            elsif($incall)
            {
                $arguments .= $continstr.$comma.$singarg ;
                # this is not a CHF_* declaration
                $chf_decl = 0;
            }
            else
            {
              print "$declarations";
              die "error: doFortranProc: unable to process argument $singarg\n";
            }
            if( $chf_decl and $singname eq "" )
            {
                # CHF_* commands must have an argument
                die "error: doFortranProc: the CHF_ command |$singarg| does not have an argument in square brackets []\n" ;
            }
        }
    }

###start printing all this crap out
    if( $ChF::multiDim == 0 ){
        print SubroutProc::FOUT $subnamecap."(" ;
    }else{
        print SubroutProc::FOUT $subnamecap."_${ChF::SpaceDim}D(" ;
    }

    print SubroutProc::FOUT $arguments;

    print SubroutProc::FOUT $continstr.")$remainder\n" ;

###print preliminary stuff that is in every procedure.
###recall indentstr includes a carraige return
    if(!$incall){
        print SubroutProc::FOUT $moduledecs;
        print SubroutProc::FOUT $indentstr."implicit none";
        ## only declare CHF_ID array if it is used
        if($debug)
        {
            print "doFortranProc: UsesCHFID $subnamecap is ",defined($SubroutProc::UsesCHFID{$subnamecap}),"\n" ;
        }
        if( defined $SubroutProc::UsesCHFID{$subnamecap} )
        {
            print SubroutProc::FOUT $indentstr."integer CHF_ID(0:5,0:5)" ;
            print SubroutProc::FOUT $indentstr."data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /\n";
            print SubroutProc::FOUT "\n" ;
        }
        ##print declarations
        print  SubroutProc::FOUT  $declarations;
        print  SubroutProc::FOUT  "\n";
    }
}

#--------------------------------------------------------------
###get name out of  string of the form blah[name]
#--------------------------------------------------------------
sub SubroutProc::getNamedThing
{
    my ($argstring) = @_;
###find where there is a [
    my $teststring = $argstring;
    ###mulitple-argument things should not be sent here.
    if($teststring =~ /;/)
    {
        die "getNamedThing got something with a ; [$argstring].\n";
    }
    my $restring ;
    if( $argstring =~ /\[/ig )
    {
#       or ) die "no \[ in argument to getNamedThing\n";
        my $beginoffset = pos $argstring;
        $restring = substr($argstring, $beginoffset);
        $restring =~ s/\[//g;
        $restring =~ s/\]//g;
    }
    else
    {
        $restring = "" ;
    }
    return $restring;
}

#--------------------------------------------------------------
#--------------------------------------------------------------
sub SubroutProc::getBoxDecl
{
    my ($singname) = @_;
    my $indentstr = $SubroutProc::indentstr;
    my $continstr = $SubroutProc::continstr;
    my $declarations = "";
    my $loarg = "CHF\_DDECL\[";
    my $hiarg = "CHF\_DDECL\[";
    my $semicolon = "";
    for(my $idir =0; $idir < 6; $idir++)
    {
        if($idir > 0)
        {
            $semicolon = ";";
        }
        $loarg .= $semicolon."i".$singname."lo".$idir;
        $hiarg .= $semicolon."i".$singname."hi".$idir;
    }
    $loarg .= "\]";
    $hiarg .= "\]";
    $declarations .= $indentstr."integer $loarg";
    $declarations .= $indentstr."integer $hiarg";
    return $declarations;
}
#--------------------------------------------------------------
#--------------------------------------------------------------
sub SubroutProc::getFabIntArgs
{
    my ($singname) = @_;
    my $indentstr = $SubroutProc::indentstr;
    my $continstr = $SubroutProc::continstr;
    my $intargs = "CHF\_DDECL\[";
    my $semicolon = "";
    my $colon = ":";
    for(my $idir =0; $idir < 6; $idir++)
    {
        if($idir > 0)
        {
            $semicolon = ";$continstr";
        }
        $intargs .= $semicolon."i".$singname."lo".$idir;
        $intargs .= $colon."i".$singname."hi".$idir;
    }
    $intargs .= "\]";
    return $intargs;
}
#--------------------------------------------------------------
#--------------------------------------------------------------
sub SubroutProc::getBoxArgs
{
    my ($singname,$comma) = @_;
    my $indentstr = $SubroutProc::indentstr;
    my $continstr = $SubroutProc::continstr;
    my $arguments = "";
    my $loarg = "CHF\_DDECL\[";
    my $hiarg = "CHF\_DDECL\[";
    my $semicolon = "";
    for(my $idir =0; $idir < 6; $idir++)
    {
        if($idir > 0)
        {
            $semicolon = ";";
        }
        $loarg .= $semicolon."i".$singname."lo".$idir;
        $hiarg .= $semicolon."i".$singname."hi".$idir;
    }
    $loarg .= "\]";
    $hiarg .= "\]";
    $arguments .= $loarg."$continstr,".$hiarg;
}
#--------------------------------------------------------------
#--------------------------------------------------------------
sub SubroutProc::doCPrototype
{
    my ($subname, $argstring,  $debug) = @_;
    my $subnamecap = uc($subname);
    my $subnamelc = lc($subname);
    my @subargs = split(",",$argstring);
    my $arguments = "";
    my $indentstr = $SubroutProc::indentstr;
###create argument list and declaration stuff
    for(my $iarg = 0; $iarg <= $#subargs; $iarg++)
    {
        my $singarg = $subargs[$iarg];
        my $singname = &SubroutProc::getNamedThing($singarg);

        if($debug)
        {
            print "doCProtoType: singarg = $singarg\n";
            print "doCProtoType: singname= $singname\n";
        }
        my $comma = "";
        if($iarg > 0)
        {
            $comma = ",";
        }
        if($singarg =~ /CHF\_REALVECT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_REALVECT\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_REALVECT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_REALVECT\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_INTVECT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_INTVECT\($singname\)";
        }
        elsif($singarg =~ /CHF\_INTVECT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_INTVECT\($singname\)";
        }
        elsif($singarg =~ /CHF\_INT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_INT\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_INT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_INT\($singname\)";
        }
        elsif($singarg =~ /CHF\_REAL/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_REAL\($singname\)";
        }
        elsif($singarg =~ /CHF\_USE/ig)
        {
            ##ignore this arg
        }
        elsif($singarg =~ /CHF\_CONST\_REAL/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_REAL\($singname\)";
        }
        elsif($singarg =~ /CHF\_COMPLEX/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_COMPLEX\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_COMPLEX/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_COMPLEX\($singname\)";
        }
        elsif($singarg =~ /CHF\_BOX/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_BOX\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FIA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FIA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_FIA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FIA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FRA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FRA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_FRA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FRA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FCA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FCA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_FCA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FCA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FRA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FRA\($singname\)";
        }
        elsif($singarg =~ /CHF\_FRA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FRA\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FCA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FCA\($singname\)";
        }
        elsif($singarg =~ /CHF\_FCA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FCA\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FIA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FIA\($singname\)";
        }
        elsif($singarg =~ /CHF\_FIA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FIA\($singname\)";
        }
        elsif($singarg =~ /CHF\_((CONST\_)?V[IRC])/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_".uc($1)."\($singname\)";
        }
        elsif($singarg =~ /CHF\_((CONST\_)?[IRC]1D)/ig)
        {
            # generate the 2nd argument for the C++ prototype that isn't in
            # the Fortran code
            $arguments .= $indentstr.$comma."CHFp\_".uc($1)."\($singname,i${singname}hi0\)" ;
        }
        elsif($singarg =~ /CHF\_((?:CONST\_)?[IRC]ARRAY)/ig)
        {
            # generate the remaining arguments for the C++ prototype that aren't
            # in the Fortran code
            my $arname = $1;
            $singname =~ / *(RANK_SPACEDIM_PLUS_)?(\d) */i;
            my $arrank = $2+0;
            if (defined($1) && $1 eq "RANK_SPACEDIM_PLUS_")
            {
                $arrank += $ChF::SpaceDim;
                my $arrankspec = $1 . $2;
                die "$arrankspec cannot be used with basic arrays."
            }
            ++$iarg;
            $subargs[$iarg] =~ / *(.*) *]/i;
            $singname = $1;
            $arguments .= $indentstr.$comma."CHFp\_".uc($arname)."\($arrank,$singname";
            for (my $i = 0; $i != $arrank; ++$i)
            {
                $arguments .= ",i${singname}lo${i}";
            }
            for (my $i = 0; $i != $arrank; ++$i)
            {
                $arguments .= ",i${singname}hi${i}";
            }
            $arguments .= "\)";
        }
        elsif($singarg =~ /CHF\_((?:CONST\_)?[IRC]CHARRAY)/ig)
        {
            # generate the remaining arguments for the C++ prototype that aren't
            # in the Fortran code
            my $arname = $1;
            $singname =~ / *(RANK_SPACEDIM_PLUS_)?(\d) */i;
            my $arrankspec = $2+0;
            if (defined($1) && $1 eq "RANK_SPACEDIM_PLUS_")
            {
                $arrankspec = $1 . $2;
            }
            ++$iarg;
            $subargs[$iarg] =~ / *(.*) *]/i;
            $singname = $1;
            $arguments .= $indentstr.$comma."CHFp\_".uc($arname)."\($arrankspec,$singname\)";
        }
        elsif($singarg =~ /CHF\_((?:CONST\_)?VECTOR)/ig)
        {
            # generate the remaining arguments for the C++ prototype that aren't
            # in the Fortran code
            (my $arname = $1) =~ s/VECTOR/RCHARRAY/i;
            $arguments .= $indentstr.$comma."CHFp\_".uc($arname)."\(1,$singname\)";
        }
        elsif($singarg =~ /CHF\_((?:CONST\_)?MATRIX)/ig)
        {
            # generate the remaining arguments for the C++ prototype that aren't
            # in the Fortran code
            (my $arname = $1) =~ s/MATRIX/RCHARRAY/i;
            $arguments .= $indentstr.$comma."CHFp\_".uc($arname)."\(2,$singname\)";
        }
        else
        {
            die
            "doCPrototype::unable to process argument $singarg\n";
        }
    }
    ####preliminary stuff
    print SubroutProc::COUT "\n" ;
    print SubroutProc::COUT "#include \"CH_Timer.H\"\n\n" ;
    print SubroutProc::COUT "// Prototype for Fortran procedure $subname ...\n//\n" ;
    print SubroutProc::COUT "#define FORT_$subnamecap FORTRAN_NAME( $subnamecap ,$subnamelc )\n" ;
    print SubroutProc::COUT "void \nFORT_$subnamecap(" ;
    print SubroutProc::COUT $arguments;
    print SubroutProc::COUT " \);\n" ;
}
#--------------------------------------------------------------
sub SubroutProc::doTimedCPrototype
{
    my ($subname, $argstring,  $debug) = @_;
    my $subnamecap = uc($subname);
    my $subnamelc = lc($subname);
    my @subargs = split(",",$argstring);
    my $arguments = "";
    my $targuments= "";
    my $indentstr = $SubroutProc::indentstr;
###create argument list and declaration stuff
    for(my $iarg = 0; $iarg <= $#subargs; $iarg++)
    {
        my $singarg = $subargs[$iarg];
        my $singname = &SubroutProc::getNamedThing($singarg);
        if($debug)
        {
            print "doCProtoType: singarg = $singarg\n";
            print "doCProtoType: singname= $singname\n";
        }
        my $comma = "";
        if($iarg > 0)
        {
            $comma = ",";
        }
        if($singarg =~ /CHF\_REALVECT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_REALVECT\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_REALVECT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_REALVECT\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_INTVECT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_INTVECT\($singname\)";
        }
        elsif($singarg =~ /CHF\_INTVECT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_INTVECT\($singname\)";
        }
        elsif($singarg =~ /CHF\_INT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_INT\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_INT/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_INT\($singname\)";
        }
        elsif($singarg =~ /CHF\_REAL/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_REAL\($singname\)";
        }
        elsif($singarg =~ /CHF\_USE/ig)
        {
            ##ignore this arg
        }
        elsif($singarg =~ /CHF\_CONST\_REAL/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_REAL\($singname\)";
        }
        elsif($singarg =~ /CHF\_COMPLEX/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_COMPLEX\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_COMPLEX/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_COMPLEX\($singname\)";
        }
        elsif($singarg =~ /CHF\_BOX/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_BOX\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FIA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FIA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_FIA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FIA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FRA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FRA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_FRA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FRA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FCA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FCA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_FCA1/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FCA1\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FRA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FRA\($singname\)";
        }
        elsif($singarg =~ /CHF\_FRA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FRA\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FCA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FCA\($singname\)";
        }
        elsif($singarg =~ /CHF\_FCA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FCA\($singname\)";
        }
        elsif($singarg =~ /CHF\_CONST\_FIA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_CONST\_FIA\($singname\)";
        }
        elsif($singarg =~ /CHF\_FIA/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_FIA\($singname\)";
        }
        elsif($singarg =~ /CHF\_((CONST\_)?V[IRC])/ig)
        {
            $arguments .= $indentstr.$comma."CHFp\_".uc($1)."\($singname\)";
        }
        elsif($singarg =~ /CHF\_((CONST\_)?[IRC]1D)/ig)
        {
            # generate the 2nd argument for the C++ prototype that isn't in
            # the Fortran code
            $arguments .= $indentstr.$comma."CHFp\_".uc($1)."\($singname,i${singname}hi0\)" ;
        }
        elsif($singarg =~ /CHF\_((?:CONST\_)?[IRC]ARRAY)/ig)
        {
            # generate the remaining arguments for the C++ prototype that aren't
            # in the Fortran code
            my $arname = $1;
            $singname =~ / *(RANK_SPACEDIM_PLUS_)?(\d) */i;
            my $arrank = $2+0;
            if (defined($1) && $1 eq "RANK_SPACEDIM_PLUS_")
            {
                $arrank += $ChF::SpaceDim;
                my $arrankspec = $1 . $2;
                die "$arrankspec cannot be used with basic arrays."
            }
            ++$iarg;
            $subargs[$iarg] =~ / *(.*) *]/i;
            $singname = $1;
            $arguments .= $indentstr.$comma."CHFp\_".uc($arname)."\($arrank,$singname";
            for (my $i = 0; $i != $arrank; ++$i)
            {
                $arguments .= ",i${singname}lo${i}";
            }
            for (my $i = 0; $i != $arrank; ++$i)
            {
                $arguments .= ",i${singname}hi${i}";
            }
            $arguments .= "\)";
        }
        elsif($singarg =~ /CHF\_((?:CONST\_)?[IRC]CHARRAY)/ig)
        {
            # generate the remaining arguments for the C++ prototype that aren't
            # in the Fortran code
            my $arname = $1;
            $singname =~ / *(RANK_SPACEDIM_PLUS_)?(\d) */i;
            my $arrankspec = $2+0;
            if (defined($1) && $1 eq "RANK_SPACEDIM_PLUS_")
            {
                $arrankspec = $1 . $2;
            }
            ++$iarg;
            $subargs[$iarg] =~ / *(.*) *]/i;
            $singname = $1;
            $arguments .= $indentstr.$comma."CHFp\_".uc($arname)."\($arrankspec,$singname\)";
        }
        elsif($singarg =~ /CHF\_((?:CONST\_)?VECTOR)/ig)
        {
            # generate the remaining arguments for the C++ prototype that aren't
            # in the Fortran code
            (my $arname = $1) =~ s/VECTOR/RCHARRAY/i;
            $arguments .= $indentstr.$comma."CHFp\_".uc($arname)."\(1,$singname\)";
        }
        elsif($singarg =~ /CHF\_((?:CONST\_)?MATRIX)/ig)
        {
            # generate the remaining arguments for the C++ prototype that aren't
            # in the Fortran code
            (my $arname = $1) =~ s/MATRIX/RCHARRAY/i;
            $arguments .= $indentstr.$comma."CHFp\_".uc($arname)."\(2,$singname\)";
        }
        else
        {
            die
            "doCPrototype::unable to process argument $singarg\n";
        }
    }
    $targuments = $arguments;
    $targuments =~ s/CHFp/CHFt/g;
    ####preliminary stuff
    print SubroutProc::COUT "\n" ;
    print SubroutProc::COUT "#ifndef GUARD$subnamecap \n";
    print SubroutProc::COUT "#define GUARD$subnamecap \n";
    print SubroutProc::COUT "// Prototype for Fortran procedure $subname ...\n//\n" ;
    print SubroutProc::COUT "void FORTRAN_NAME( $subnamecap ,$subnamelc )(" ;
    print SubroutProc::COUT $arguments ;
    print SubroutProc::COUT " );\n\n" ;

    print SubroutProc::COUT "#define FORT_$subnamecap FORTRAN_NAME( inline$subnamecap, inline$subnamecap)\n";
    print SubroutProc::COUT "#define FORTNT_$subnamecap FORTRAN_NAME( $subnamecap, $subnamelc)\n\n";

    print SubroutProc::COUT "inline void FORTRAN_NAME(inline$subnamecap, inline$subnamecap)(" ;
    print SubroutProc::COUT $arguments;
    print SubroutProc::COUT " )\n{\n" ;
    print SubroutProc::COUT " CH_TIMELEAF(\"FORT_$subnamecap\");\n" ;
    print SubroutProc::COUT " FORTRAN_NAME( $subnamecap ,$subnamelc )(" ;
    print SubroutProc::COUT $targuments ;
    print SubroutProc::COUT " );\n" ;
    print SubroutProc::COUT "}\n" ;
    print SubroutProc::COUT "#endif  // GUARD$subnamecap \n";





}
# -------------------------------------------------------------------
# Subroutine: startFFile
# -------------------------------------------------------------------

sub SubroutProc::startFFile
{
    print SubroutProc::FOUT "#include \"REAL.H\"\n" ;
    print SubroutProc::FOUT "#include \"SPACE.H\"\n" ;
    print SubroutProc::FOUT "#include \"CONSTANTS.H\"\n" ;
    print SubroutProc::FOUT "\n" ;
}

# -------------------------------------------------------------------
# Subroutine: startCFile
# -------------------------------------------------------------------

sub startCFile
{
    # Chombo convention for the cpp macro that protects H file is _NAME_H_
    # so convert the filename to uppercase
    tr/a-z/A-Z/ for @_ ;
    my ($basename) = @_;
    if( $ChF::debug){ print "startcfile basename = $basename \n"; }
    print SubroutProc::COUT "#ifndef \_".$basename."\_F\_H\_\n" ;
    print SubroutProc::COUT "#define \_".$basename."\_F\_H\_\n" ;
    print SubroutProc::COUT "\n" ;
    print SubroutProc::COUT "#include \"FORT_PROTO.H\"\n" ;
    print SubroutProc::COUT "#include \"CH_Timer.H\"\n" ;
    print SubroutProc::COUT "#include \"REAL.H\"\n" ;
    print SubroutProc::COUT "\n" ;
#XXX -- dont do this.  Assume the _multidim_F.H file does it. <dbs>
#XXX    if( $ChF::multiDim != 0 ){
#XXX        print SubroutProc::COUT "#ifdef CH_MULTIDIM\n" ;
#XXX        my $hash_if = "#if  " ;
#XXX        for( my $idim=1 ; $idim<=$ChF::MAXDIM ; ++$idim ){
#XXX            print SubroutProc::COUT "${hash_if} CH_SPACEDIM == ${idim}\n";
#XXX            print SubroutProc::COUT "using namespace CH_${idim}D ;\n" ;
#XXX            $hash_if = "#elif" ;
#XXX        }
#XXX        print SubroutProc::COUT "#else\n" ;
#XXX        print SubroutProc::COUT "#error unsupported CH_SPACEDIM value\n" ;
#XXX        print SubroutProc::COUT "#endif\n" ;
#XXX        print SubroutProc::COUT "#endif\n" ;
#XXX        print SubroutProc::COUT "\n" ;
#XXX    }
    print SubroutProc::COUT "extern \"C\"\n" ;
    print SubroutProc::COUT "{\n" ;
}


# -------------------------------------------------------------------
# Subroutine: finishCFile
# -------------------------------------------------------------------

sub SubroutProc::finishCFile
{
    print SubroutProc::COUT "\n}\n\n" ;
    print SubroutProc::COUT "#endif\n" ;
}

sub SubroutProc::findCHFID
{
    my ($debug) = @_;
    ## find "subroutine" statements, remember the subroutine name, read the
    ## rest of the routine (until EOF or the next "subroutine" statement)
    ## and if CHF_ID() is used, put the subroutine name in the UsesCHFID hash
    my $subname = "" ;
    while( <SubroutProc::FIN> ){
        ## NOTE: this will fail on a statement like "subroutinefoo=1"
        if (m/^\s*subroutine/i)
        {
            if($debug)
            {
                print "findCHFID: subroutine line = $_\n";
            }
            # remove the keyword from the line
            s/^\s*subroutine\s*//;
            # find the subroutine name and remember the uppercase version
            m/([^\s]+)\s*\(/ ;
            if($debug)
            {
                print "findCHFID: subr name is [$1]\n";
            }
            $subname = uc($1);
        }
        elsif (m/chf_id\s*\(/i)
        {
            $SubroutProc::UsesCHFID{$subname} = 1 ;
            if($debug)
            {
                print "findCHFID: subroutine [$subname] uses CHF_ID.\n";
            }
        }
        elsif (m/chf_autoid\s*\[/i)
        {
            $SubroutProc::UsesCHFID{$subname} = 1 ;
            if($debug)
            {
                print "findCHFID: subroutine [$subname] uses CHF_AUTOID.\n";
            }
        }
    }
    seek( SubroutProc::FIN,0,0 ); #set file pointer to start-of-file
}

sub SubroutProc::getList
{
    my ($buf) = @_;
    my $str ;
    # skip until the first paren or quote
start:
    if( $buf =~ m/[\(\[\{\'\"]/g )
    {
        $buf = $& . $' ; #'
        # if a quote, skip the string and try again
        # else extract the parenthesized list
        if( $& eq "\'" or $& eq "\"" )
        {
            ($str ,$buf) = Text::Balanced::extract_delimited( $buf, "\'\"" ) ;
            goto start ;
        }
        else
        {
            $buf = Text::Balanced::extract_bracketed( $buf,"()[]{}" );
        }
    }
    else
    {
        # nothing to extract
        $buf = "" ;
    }
    return $buf;
}

## This is here because the Perl module loading mechanism expects the
## module to return a true value.  Otherwise the program aborts.
1;

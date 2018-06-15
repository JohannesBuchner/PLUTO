#!/usr/local/gnu/bin/perl -w

####  _______              __
#### / ___/ /  ___  __ _  / /  ___
####/ /__/ _ \/ _ \/  ' \/ _ \/ _ \
####\___/_//_/\___/_/_/_/_.__/\___/ 
####
####
#### This software is copyright (C) by the Lawrence Berkeley
#### National Laboratory.  Permission is granted to reproduce
#### this software for non-commercial purposes provided that
#### this notice is left intact.
#### 
#### It is acknowledged that the U.S. Government has rights to
#### this software under Contract DE-AC03-765F00098 between
#### the U.S.  Department of Energy and the University of
#### California.
####
#### This software is provided as a professional and academic
#### contribution for joint exchange. Thus it is experimental,
#### is provided ``as is'', with no warranties of any kind
#### whatsoever, no support, no promise of updates, or printed
#### documentation. By using this software, you acknowledge
#### that the Lawrence Berkeley National Laboratory and
#### Regents of the University of California shall have no
#### liability with respect to the infringement of other
#### copyrights by any part of this software.
####
####
#########################################
#### This is an augmented version of ChomboFortran.
#### It also removes the thrice-accursed D_TERM
####  DTGraves July 2000
#########################################
package ChF;
use strict;
BEGIN
{
##global variables god help us all
##all have the package name ChF    
    $ChF::CHFPP_Version = "1.23" ;
    $ChF::FComment = "C";
    $ChF::debug  = 0; 
    $ChF::FINFile   = "";
    $ChF::FOUTFile = "";
    $ChF::baseName = "";
    $ChF::basePath = "";
    
    
###parse the command line
###and generate appropriate file names
    &ChF::parseInputs();
    
###need two temp files because final output can be /dev/null
###Strip comments because they are evil
    require "stripsharp.pm";
    StripSharpProc->import();
    unless(my $flag = &StripSharpProc::StripSharp($ChF::debug))
    {
        die "problem in sharp stripper\n";
    }
    
###blow town
    exit 0 ;
    
    
    
#########################################
# Subroutine: parseInputs
# Inputs: none
# Returns: none
# list goes like this:
# -f <inputfilename>  (optional) -p <output filename> 
# if the fortran file output name is not specified it goes to baseName.F
# if the c header  output name is not specified it goes to baseName_F.H
# Globals (all outputs): 
# FINFile input file name
# FOUTFile fortran file name for output
# COUTFile c header file name for output
# baseName base name of input file
# basePath
#########################################
# -------------------------------------------------------------------
    sub ChF::parseInputs
    {

        package ChF::parseInputs;

        my $l_debug = 0;
        use Getopt::Std;
        
#if -d then debugging is enabled
        my %l_option;
        
      Getopt::Std::getopts("dc:f:p:D:" , \%l_option);
        if($l_option{d})
        {
            $l_debug = 1;
        }

        if($l_debug)
        {
            print "starting chombo fortran preprocessing ";
            print "with debug enabled\n";
        }
        $ChF::debug = $l_debug;

    }
}

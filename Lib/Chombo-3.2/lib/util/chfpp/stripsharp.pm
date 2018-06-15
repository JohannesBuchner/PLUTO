#!/usr/local/gnu/bin/perl -w
####
####    _______              __
####   / ___/ /  ___  __ _  / /  ___
####  / /__/ _ \/ _ \/  V \/ _ \/ _ \
####  \___/_//_/\___/_/_/_/_.__/\___/
####  Please refer to Copyright.txt, in Chombo's root directory.
####

#################################################
###  This is the DTerm processer.
### Interface is
### sub DTermProc::procDTermMacros(inputfile, outputfile,
###                                SpaceDim, debug)
###
###  reads in input file 
#################################################

package StripSharpProc;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
$VERSION = 1.23;
@ISA = qw(Exporter);
@EXPORT = qw(&stripSharp);
@EXPORT_OK = qw(&stripSharp);
@EXPORT_TAGS= ();

sub StripSharpProc::StripSharp
{

    use strict;
    my ($debug) = @_;
    if($debug)
    {
        print "StripSharpProc: \n";
    }
    
    while (defined(my $ibuf = <STDIN> )) 
    {

###     skip lines that start with # 
###     (or any number of spaces before it) -JNJ
        if($ibuf =~ m/^. \#/i)
        {
            next;
        }
        else
        {
            print $ibuf;
        }
    }
    
   return 1;   
}
###i have no idea why this is here.
###the perl cookbook book told me to put it there.
###really.
1;

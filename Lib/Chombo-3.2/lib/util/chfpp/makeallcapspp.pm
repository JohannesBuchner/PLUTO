#!/usr/local/gnu/bin/perl -w
####
####    _______              __
####   / ___/ /  ___  __ _  / /  ___
####  / /__/ _ \/ _ \/  V \/ _ \/ _ \
####  \___/_//_/\___/_/_/_/_.__/\___/
####  Please refer to Copyright.txt, in Chombo's root directory.
####

#################################################
###  This is the MultiDo processer.
### Interface is
### sub MultiDoProc::procMultiDoMacros(inputfile, outputfile,
###                                SpaceDim, debug)
###
###  reads in input file 
#################################################

package MakeAllCapsProc;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
$VERSION = 1.23;
@ISA = qw(Exporter);
@EXPORT = qw(&procMakeAllCaps);
@EXPORT_OK = qw(&procMakeAllCaps);
@EXPORT_TAGS= ();

sub MakeAllCapsProc::makeAllCaps
{
    use strict;
    my ($FINFile, $FOUTFile, $SpaceDim, $debug) = @_;
    if($debug)
    {
        print "MakeAllCapsProc: \n";
        print "input file  = $FINFile \n";
        print "output file = $FOUTFile \n";
        print "SpaceDim    = $SpaceDim \n";
    }
    
    open(FOUT,">" . $FOUTFile) 
        or die "Error: cannot open output file " . $FOUTFile . "\n";
    open(FIN,"<" . $FINFile) 
        or die "Error: cannot open input file " . $FINFile . "\n";

    my $ibuf = "";
    my $obuf = "";
    ###put the entire input file into a string buffer
    while (defined( $ibuf = <FIN> )) 
    {
        ##only uppercase lines that do not start in # 
        if($ibuf =~ m/^\#/i)
        {
            $obuf .= $ibuf;
        }
        else
        {
            $obuf .= uc($ibuf);
        }
    }
    
    print FOUT $obuf;

    ###need to close files to make this modular.
    close(FIN);
    close(FOUT);
    return 1;   
}
###i have no idea why this is here.
###the perl cookbook book told me to put it there.
###really.
1;

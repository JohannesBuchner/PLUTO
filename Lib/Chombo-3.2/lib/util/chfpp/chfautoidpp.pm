#!/usr/local/gnu/bin/perl -w
####
####    _______              __
####   / ___/ /  ___  __ _  / /  ___
####  / /__/ _ \/ _ \/  V \/ _ \/ _ \
####  \___/_//_/\___/_/_/_/_.__/\___/
####  Please refer to Copyright.txt, in Chombo's root directory.
####

#################################################
###  This is the Chfautoid processer.
### Interface is
### sub ChfAutoidProc::procChfAutoidMacros(inputfile, outputfile,
###                                SpaceDim, debug)
###
### The specification for MultiDo is:
###    CHF_AUTOID[var;dir{;factor}]
### where {;factor} is optional, defaulting to 1.
###
### Method: reads in input file 
#################################################

package ChfAutoidProc;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
$VERSION = 1.23;
@ISA = qw(Exporter);
@EXPORT = qw(&procChfAutoidMacros);
@EXPORT_OK = qw(&procChfAutoidMacros);
@EXPORT_TAGS= ();

sub ChfAutoidProc::procChfAutoidMacros
{

    use strict;
    my ($FINFile, $FOUTFile, $SpaceDim, $debug) = @_;
    if($debug)
    {
        print "ChfAutoidProc: \n";
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
        $obuf .= $ibuf;
    }
    my $indentstr = "\n      ";
    my $beginstring = "CHF_AUTOID";
    my $beginlen = 10;
    my $endstring = "]";
    my $endlen = 1;
    my $beginoffset = 0;
    my $endoffset = 0;
    my $length = 0;
    my $firstlen = 0;
    my $firststring = " ";
    my $newstring = " ";
    while ($obuf  =~ m/$beginstring/ig )
    {
        $beginoffset = pos $obuf;

        $firstlen = $beginoffset-$endoffset-$beginlen;
        $firststring = substr($obuf, $endoffset, $firstlen);
        print FOUT $firststring;

        my $tempbuf = $obuf;
        pos $tempbuf = $beginoffset;
        $tempbuf =~ m/$endstring/ig;
        $endoffset = pos $tempbuf;
        $length = $endoffset-$beginoffset-$endlen;

        $newstring = substr($obuf, $beginoffset, $length);
        ###newstring now is of form [var;dir] or [var;dir;factor]
        $newstring =~ s/\[//g;
        $newstring =~ s/\]//g;
        ###newstring now is of form var;dir or var;dir;factor without whitespace
        my @dimarg = split(";", $newstring);
        my $printstring ="";
        ### variable name is first in the list, so remove it
        my $varname = shift @dimarg;
        ### index of direction is next in the list, so remove it
        my $dirind = shift @dimarg;
        ### optional factor is next in the list, so remove it
        my $factor = shift @dimarg;
        for( my $idir=0 ; $idir < $SpaceDim ; $idir++ ) 
        {
            $printstring .= $indentstr."$varname$idir=";
            if (defined($factor))
            {
                $printstring .= "$factor*";
            }
            $printstring .= "CHF_ID($idir,$dirind)\n";
        }
        print FOUT $printstring;
    }
    $newstring = substr($obuf, $endoffset);
    print FOUT $newstring;

    ###need to close files to make this modular.
    close(FIN);
    close(FOUT);
    return 1;   
}
###i have no idea why this is here.
###the perl cookbook book told me to put it there.
###really.
1;

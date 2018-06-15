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

package DTermProc;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
$VERSION = 1.23;
@ISA = qw(Exporter);
@EXPORT = qw(&procDTermMacros);
@EXPORT_OK = qw(&procDTermMacros);
@EXPORT_TAGS= ();

sub DTermProc::procDTermMacros
{

    use strict;
    my ($FINFile, $FOUTFile, $SpaceDim, $debug) = @_;
    if($debug)
    {
        print "DTermProc: \n";
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
    
    my $beginstring = "CHF_DTERM";
    my $beginlen = 9;
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
        ###newstring now is of form [arg0;arg1;arg2]
        $newstring =~ s/\[//g;
        $newstring =~ s/\]//g;
        ###newstring now is of form arg0;arg1;arg2
        my @dimarg = split(";", $newstring);
        my $printstring ="";
        for( my $idir=0 ; $idir <= $#dimarg ; $idir++ ) 
        {
            if($idir < $SpaceDim)
            {
                $printstring .= $dimarg[$idir];
            }
 ##          print "$idir  $dimarg[$idir] \n";
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

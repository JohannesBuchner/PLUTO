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
### The specification for MultiDo is:
###    CHF_MULTIDO[boxvar;indexvar1;indexvar2;indexvar3{;stride}]
###        ....
###    CHD_ENDDO
### where {;stride} is optional and must be an integer constant.
### It defaults to 1.
### It may be negative, in which case the loop runs backwards.
###
###
### Method: reads in input file 
#################################################

package MultiDoProc;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
$VERSION = 1.24;
@ISA = qw(Exporter);
@EXPORT = qw(&procMultiDoMacros);
@EXPORT_OK = qw(&procMultiDoMacros);
@EXPORT_TAGS= ();

sub MultiDoProc::procMultiDoMacros
{

    use strict;
    my ($FINFile, $FOUTFile, $SpaceDim, $debug) = @_;
    if($debug)
    {
        print "MultiDoProc: \n";
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
    
    ###replace CHF_ENDDO with the right number of enddos
    my $indentstr = "\n      ";
    my $enddostring = "";
    for(my $idir = 0; $idir < $SpaceDim; $idir++)
    {
        $enddostring .= $indentstr."enddo";
    }
    $obuf =~ s/CHF\_ENDDO/$enddostring/ig;

    my $beginstring = "CHF_MULTIDO";
    my $beginlen = 11;
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
        ###newstring now is of form [arg0;arg1;arg2] or [arg0;arg1;arg2;stride]
        $newstring =~ s/\[//g;
        $newstring =~ s/\]//g;
        $newstring =~ s/\s//g;
        ###newstring now is of form arg0;arg1;arg2 or arg0;arg1;arg2;stride without white space
        my @dimarg = split(";", $newstring);
        my $printstring ="";
        ###box variable is first in the list, so remove it
        my $boxname = shift @dimarg ;
        my $stride = 1;
        ###if the last argument is a number, then assume it is the optional stride
        if( $dimarg[$#dimarg] =~ m/^[\+\-]?[0-9]+$/ ){
            $stride = pop @dimarg ;
            if( $stride == 0 ){
                die "error: CHF_MULTIDO with stride 0 is not allowed\n";
            }
        }
        for( my $idir=$#dimarg ; $idir >= 0 ; $idir-- ) 
        {
            if($idir < $SpaceDim)
            {
                $printstring .= $indentstr."do $dimarg[$idir] = ";
                ## add the stride if not 1 and reverse the bounds if <0
                if( $stride > 0 ){
                    $printstring .= "CHF\_LBOUND\[$boxname; $idir\],";
                    $printstring .= "CHF\_UBOUND\[$boxname; $idir\]";
                    if( $stride != 1 ){
                        $printstring .= ",$stride" ;
                    }
                }else{
                    $printstring .= "CHF\_UBOUND\[$boxname; $idir\],";
                    $printstring .= "CHF\_LBOUND\[$boxname; $idir\]";
                    $printstring .= ",$stride" ;
                }
            }
        }
        $printstring .= "\n";
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

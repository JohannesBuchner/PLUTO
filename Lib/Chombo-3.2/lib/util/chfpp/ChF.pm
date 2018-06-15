####
####    _______              __
####   / ___/ /  ___  __ _  / /  ___
####  / /__/ _ \/ _ \/  V \/ _ \/ _ \
####  \___/_//_/\___/_/_/_/_.__/\___/
####  Please refer to Copyright.txt, in Chombo's root directory.
####

package ChF;

BEGIN
{
##global variables god help us all
##all have the package name ChF
    $ChF::CHFPP_Version = "1.60" ;
    $ChF::FComment = "C";
    $ChF::MAXDIM   = 6;
    $ChF::debug    = 0;
    $ChF::FINFile  = "";
    $ChF::FOUTFile = "";
    $ChF::COUTFile = "";
    $ChF::baseName = "";
    $ChF::basePath = "";
    $ChF::SpaceDim = -1;
    $ChF::multiDim = 0 ;
#   The $$ is the process ID and this will work better
#   for compiling with multiple threads.  (bvs, dbs, ndk)
#    $ChF::T1File = "chfpp03079.tmp";
#    $ChF::T2File = "chfpp94720.tmp";
    $ChF::T1File = "$$.tmp1";
    $ChF::T2File = "$$.tmp2";

}

#########################################
# Subroutine: parseInputs
# Inputs: none
# Returns: none
# list goes like this:
# -f <inputfilename>  (optional) -p < fort output filename> -c <coutfilename>
# -D <SpaceDim> [-m]
# if the fortran file output name is not specified it goes to baseName.F
# if the c header  output name is not specified it goes to baseName_F.H
# if -m is given then multidimension mode is enabled
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
    use File::Basename;
###        package ChF::parseInputs;
    my $l_FINFile = "";
    my $l_FOUTFile= "";
    my $l_COUTFile= "";
    my $l_baseName= "";
    my $l_basePath = "";
    my $l_debug = 0;
    my $l_SpaceDim = -1;
    use Getopt::Std;

#if -d then debugging is enabled
    my %l_option;

    Getopt::Std::getopts("dc:f:mp:D:" , \%l_option);
    if($l_option{d})
    {
        $l_debug = 1;
    }
    if ($l_option{m})
    {
        # enable Ted's multiple-dimensions feature
        $ChF::multiDim = 1 ;
    }
    if ($l_option{f})
    {
        $l_FINFile = $l_option{f};
    }
    else
    {
        die "no input file specified\n";
    }

    if ($l_option{D})
    {
        $l_SpaceDim = $l_option{D};
        if($l_SpaceDim > $ChF::MAXDIM)
        {
            die "bogus spacedim specified= ".$l_SpaceDim.
                "maxSpaceDim=".$ChF::MAXDIM;
        }
        if(($l_SpaceDim < 1))

        {
            die "negative spacedim specified= ".$l_SpaceDim;
        }

    }
    else
    {
        die "no SpaceDim specified (-D option)\n";
    }

    use File::Basename;
    ($l_baseName,$l_basePath) = fileparse( $l_FINFile ,'\..*' ) ;
    if ($l_option{c})
    {
        $l_COUTFile = $l_option{c};
    }
    else
    {
        $l_COUTFile = $l_baseName . "_F.H" ;
    }
    if ($l_option{p})
    {
        $l_FOUTFile = $l_option{p};
    }
    else
    {
        $l_FOUTFile = $l_baseName . ".F" ;
    }
    if($l_debug)
    {
        print "starting chombo fortran preprocessing ";
        print "with debug enabled\n";
        print "input file  = $l_FINFile\n";
        print "baseName= $l_baseName \n";
        print "basePath= $l_basePath \n";
        if( $ChF::multiDim != 0 ){ print "multidimensional functionality enabled\n"; }
        print "outputting fortran source to $l_FOUTFile \n";
        print "outputting header for c++ to $l_COUTFile \n";
    }
    $ChF::debug = $l_debug;
    $ChF::FINFile = $l_FINFile;
    $ChF::FOUTFile = $l_FOUTFile;
    $ChF::COUTFile = $l_COUTFile;
    $ChF::baseName = $l_baseName;
    $ChF::basePath = $l_basePath;
    $ChF::SpaceDim = $l_SpaceDim;
}

1;

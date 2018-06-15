#!/usr/bin/env perl

use strict;
use warnings;

my $compiler;
if ($#ARGV != 0) {
   #warning "usage: ./compilerVersion [compiler]  ignoring others\n";
}

$compiler = $ARGV[0];

my $found=0;
my $version=0;

#current matching criteria on the output of compiler version command:
# simply look for the _first_ " digit.digit " or " digit.digit.digit " or "digit.digit-digit "
# will not match " digit "  (eg a year)
# will not match "digit.digit "  or " digit.digit" -- some 'space' char needs to be on either side.
#     (ie v=10.1 would be missed   or a version number within a pathname )
# will not match a version number with a non-digit (ie version 10.3b would be missed)

# Test --version first
my $tv = `$compiler 2>&1 --version`;
if (length($tv) < 1024) { # try to filter command returning a lot of garbage (or the man page)
   if ($tv =~ /\s(\d+\.\d+(|(\.|-)\d+))\s/) {
      #print "version=($1)\n";
      $version=$1;
      $found=1;
   } else {
      #print "($v)\n";
   }
}

# if nothing found, try -V
if ($found==0) {
   # -V
   $tv = `$compiler 2>&1 -V`;
   if (length($tv) < 1024) { # try to filter command returning a lot of garbage (or the man page)
      if ($tv =~ /\s(\d+\.\d+(|(\.|-)\d+))\s/) {
         #print "version=($1)\n";
         $version=$1;
         $found=1;
      } else {
         #print "($v)\n";
      }
   }
}

# if nothing found, try -qversion (for AIX, xlc)  (add an extra ".digit" so version can be d.d.d.d)
if ($found==0) {
   # -qversion
   $tv = `$compiler 2>&1 -qversion`;
   if (length($tv) < 1024) { # try to filter command returning a lot of garbage (or the man page)
      if ($tv =~ /\s(\d+\.\d+(|\.\d+)(|\.\d+))\s/) {
         #print "version=($1)\n";
         $version=$1;
         $found=1;
      } else {
         #print "($v)\n";
      }
   }
}

if ($found==1) {
   print "$version\n";
} else {
   print "version not found.\n";
}

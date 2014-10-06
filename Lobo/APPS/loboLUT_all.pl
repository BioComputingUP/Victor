#!/usr/bin/perl -w
##  This file is part of Victor.

##    Victor is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    Victor is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

##    You should have received a copy of the GNU General Public License
##    along with Victor.  If not, see <http://www.gnu.org/licenses/>.
##
## -*- Perl -*-----------------------------------------------------------------
##  $Id: loboLUT_all,v 1.7 2007-12-17 09:23:14 biocomp Exp $
##
##  Author:             Silvio Tosatto
##
##  Project Name:       Lobo
##
##  Date:               08/03
##
##  Description:
##    This script serves to create an ensemble of LUTs for use with the "lobo"
##    set of programs. 
##
## ---------------------------------------------------------------------------
use Getopt::Std;

$| = 1;
my ($script) = ($0 =~ m|([^/]*)$|);

$Use = "$script
This script serves to to create an ensemble of LUTs for use with the \"lobo\" 
set of programs.
 Options: 
\t -a <len> \t Systematically create _all_ LUTs of size 2 to <len>
\t -c <len>  \t Create _only_ LUT of size <len> and its sub-LUTs
\t [-o]       \t Force _overwrite_ of existing LUTs
\t [-p <path>] \t Change LUT depository from default to <path>
\t [-m \"<text>\"]\t More options to pass on to loboLUT 
\n";


$root = $ENV{"VICTOR_ROOT"} . "/bin";

#
# START of main program
#

unless ($ARGV[0])
{
   die "Missing file specification. Aborting. (-h for help)\n";
}

getopts("ha:c:m:op:", \%args) or die "Aborting. (-h for help)\n"; 

$overwrite = $args{o};
$LUTpath = $ENV{"VICTOR_ROOT"} . "data/";

if ($args{p})
{
    $LUTpath = $args{p};
    print "Setting LUT depository path to $LUTpath ...\n";
}

$LUTopt = " --table $LUTpath";

if ($args{m})
{
    $LUTopt = $LUTopt . " $args{m}";
    print "Passing the following options to loboLUT: \t $args{m}\n";
}

if ($args{h}) 
{
   die "$Use";
}

if ($args{a})
{
   $max = $args{a};
   if ($max < 1)
   {
	die "Maximum value must be positive\n";
   }
   print "Creating LUTs of size up to $max ...\n"; 

   for ($i = $max; $i > 3; $i--)
   {
     push(@list, $i);
   }

   createLUT();

}
elsif ($args{c})
{
   $max = $args{c};
   if ($max < 1)
   {
	die "Maximum value must be positive\n";
   }
   print "Creating LUT of size $max, checking dependencies ...\n"; 

   push(@list,$max);
   $counter = 0;
   $maxlist = 1;
   while ($counter < $maxlist)
   {
     $temp = $list[$counter];
     $tmp1 = int ($temp / 2) + ($temp % 2);
     if (($tmp1 > 3) && ($tmp1 != $list[$#list]))
     {
       push(@list, $tmp1);
       $maxlist++;
       print "adding $tmp1 \t I  \t @list \n";
     }
     $tmp2 = int ($temp / 2);
     if (($tmp1 != $tmp2) && ($tmp2 > 3))
     {
       push(@list, $tmp2);
       $maxlist++;
	print "adding $tmp2 \t II \t @list \n";
     }
    $counter++;
   }
    
   createLUT();
   

}
else
{
   print "No valid option selected. (-h for help)\n"; 
}

print "Done.\n";

#
# END of main program
#


sub createLUT
{
    if (($max >= 2) && (askWriteLUT(2)))
    {
	print "Creating LUT aa2.lt\n";
	print `$root/loboLUT -A 1 -B 1 -O aa2.lt $LUTopt`;
    }
    
    if (($max >= 3) && (askWriteLUT(3))) 
    {
	print "Creating LUT aa3.lt\n";
	print `$root/loboLUT -A aa2.lt -B 1 -O aa3.lt $LUTopt`;
    }
    
    for ($j = $#list; $j >= 0; $j--)
    {
	if (askWriteLUT($list[$j]))
	{
	    $i = $list[$j];
	    $st = int ($i / 2) + ($i % 2);
	    $en = int ($i / 2);
	    print "Creating LUT aa$i.lt from aa$st.lt and aa$en.lt\n";
	    print `$root/loboLUT -A aa$st.lt -B aa$en.lt -O aa$i.lt $LUTopt`;
	}
    }
}

sub askWriteLUT($$)
{
    if (defined($overwrite) || !(-s "$LUTpath/aa@_.lt"))
    {
	return 1;
    }

    return 0;
}

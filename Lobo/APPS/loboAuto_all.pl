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
##  $Id: loboAuto_all,v 1.15 2007-12-17 09:23:14 biocomp Exp $
##
##  Author:             Silvio Tosatto
##
##  Project Name:       Nazgul
##
##  Date:               08/03
##
##  Description:
##    This script serves to automatically model (with loboAuto) all indels in
##    a given protein structure. 
##
## ---------------------------------------------------------------------------
use Getopt::Std;

#$MAX_SCWRL_WAIT = 300; # seconds to wait before SCWRL timeout

$| = 1;
my ($script) = ($0 =~ m|([^/]*)$|);

$Use = "$script
This script serves to automatically model (with loboAuto) all indels in
a given protein structure.
 Options: 
\t -I <file> \t Input INDEL file
\t -i <file> \t Input PDB file (single chain)
\t -o <file> \t Output PDB file
\t -s <file> \t Sequence file (in one-letter code)
\t [-b <path>] \t Change directory for temporary files to <path>
\t [-p <path>] \t Change LUT depository basename from default to <path>
\t [-m \"<text>\"]\t More options to pass on to loboAuto 
\t [-S]       \t Skip SCWRL side chain placement step
\t [-c]       \t Clean-up temporary files
\t [-v]       \t Verbose mode
\n";

$root = $ENV{"VICTOR_ROOT"} . "bin";

#
# START of main program
#

unless ($ARGV[0])
{
   die "Missing file specification. Aborting. (-h for help)\n";
}

getopts("hI:m:p:vi:o:s:cb:S", \%args) or die "Aborting. (-h for help)\n"; 

$LUTpath = $ENV{"VICTOR_ROOT"} . "data/aa";
if ($args{p})
{
    $LUTpath = $args{p};
    print "Setting LUT depository path to $LUTpath ...\n";
}

$LUTopt = " --table $LUTpath --maxWrite 10";

if ($args{m})
{
    $LUTopt = $LUTopt . " $args{m}";
    print "Passing the following options to loboAuto: \t $args{m}\n";
}

if ($args{h}) 
{
   die "$Use";
}

checkInput();

$infile = $args{I};
 print `echo "El archivo para loboAuto_all es: $infile"`;
$basepath = "./";
if ($args{b})
{
    $basepath = $args{b};
}

@idT = split(m/\//,$infile);       # get file basename: first remove path
@idT2 = split(m/\./,$idT[$#idT]);  # then remove extension
$id = $basepath . $idT2[0] . "_$$"; # add process id to name
 print `echo "El process id es: $id"`;

print '*-' x 30, "*\n";

readIndels();

if ($args{v})
{
    printContent();
}

doLoboAutoDel();
doLoboAutoIns();

print `cp $id\_insall.pdb $args{o}`;

if ($args{c})
{
    print "Cleaning up temporary files ...\n";
    print `rm $id\_*`;
}

print "Done.\n";
print '*-' x 30, "*\n";

#
# END of main program
#


# Helper methods follow:

sub readIndels
{
    open (IN,$infile) or die "Error: Indel file not found.\n";
    while (defined ($inline = <IN>))
    {
	chomp($inline);
	@inelem = split(/\W+/,$inline);
#	print "$inline\n";
#	print "$inelem[0]+++$inelem[1]+++$inelem[2]+++$inelem[3]+++\n";
	
	if ($inelem[3] eq "HEAD")
	{
	    push(@head, $inelem[1]);
	    push(@head, $inelem[2]);
	}
	elsif ($inelem[3] eq "DEL")
	{
	    push(@del, $inelem[1]);
	    push(@del, $inelem[2]);
	}
	elsif ($inelem[3] eq "INS")
	{
	    push(@ins, $inelem[1]);
	    push(@ins, $inelem[2]);
	}
	elsif ($inelem[3] eq "TAIL")
	{
	    push(@tail, $inelem[1]);
	    push(@tail, $inelem[2]);
	}
	else
	{
	    die "Error: Invalid keyword found in input file.\n";
	}
    }
    close IN;
}

sub doLoboAutoDel
{
    print "Modelling all deletions ...\n";

    print `cp $args{i} $id\_del0.pdb`;
    print `cp $args{i} $id\_del1.pdb`;

    for ($i = 0; $i <= $#del; $i += 2) 
    {
	print '*-' x 30, "*\n";
	$j = int ($i / 2);
	$k = $j + 1;

	print `cp $id\_del$j.pdb $id\_del$k.pdb`;

	if ($del[$i+1] > 1)
	{
	    $l = $del[$i+1] + 1;
	     print `echo "$root/loboLUT_all -c $l"`;

	    print `$root/loboLUT_all -c $l`; 
	    #creating LUTs, where needed
	}
  print `echo "$root/loboAuto -i $id\_del$j.pdb -s $del[$i] -l $del[$i+1] --del --seq $args{s} -o $id\_resD$i.pdb --scwrl $id\_resD$i.scwrl $LUTopt"`;

 	print `$root/loboAuto -i $id\_del$j.pdb -s $del[$i] -l $del[$i+1] --del --seq $args{s} -o $id\_resD$i.pdb --scwrl $id\_resD$i.scwrl $LUTopt`;

	if (-s "$id\_resD$i.pdb.000")
	{
	    print `cp $id\_resD$i.pdb.000 $id\_del$k.pdb`;
	    if ($args{v})
	    {
		print `/home/layla/data/SCWRL4/Scwrl4 -i $id\_resD$i.pdb.000 -o $id\_del$k.pdb -s $id\_resD$i.scwrl`;
	    }
	    else
	    {
		`/home/layla/data/SCWRL4/Scwrl4 -i $id\_resD$i.pdb.000 -o $id\_del$k.pdb -s $id\_resD$i.scwrl`;
	    }
	}
	else
	{
	    print "Warning: No model produced by loboAuto.\n";
	}
    }

    $numdel = int ($#del / 2) + 1;
    print `cp $id\_del$numdel.pdb $id\_delall.pdb`;
    print '*-' x 30, "*\n";
}


sub doLoboAutoIns
{
    print "Modelling all insertions ...\n";
    
    print `cp $id\_delall.pdb $id\_ins0.pdb`;
    print `cp $id\_delall.pdb $id\_ins1.pdb`;
    
    for ($i = 0; $i <= $#ins; $i += 2) 
    {
	print '*-' x 30, "*\n";
	$j = int ($i / 2);
	$k = $j + 1;
	print `cp $id\_ins$j.pdb $id\_ins$k.pdb`;
	
	if ($ins[$i+1] > 1)
	{   print `echo "El lobo all: $ins[$i+1] "`; 
	    $l = $ins[$i+1] + 1;
	    	

            print `echo "$root/loboLUT_all -c $l"`;
	    print `$root/loboLUT_all -c $l`; 
	    # creating LUTs, where needed
	}
 if($ins[$i+1] > 1 ){
  print `echo "$root/loboAuto -i $id\_ins$j.pdb -s $ins[$i] -l $ins[$i+1] --ins --seq $args{s} -o $id\_resI$i.pdb --scwrl $id\_resI$i.scwrl $LUTopt\n"`;
     
	print `$root/loboAuto -i $id\_ins$j.pdb -s $ins[$i] -l $ins[$i+1] --ins --seq $args{s} -o $id\_resI$i.pdb --scwrl $id\_resI$i.scwrl $LUTopt\n`;

	
 	if (-s "$id\_resI$i.pdb.000")
        {
            print `cp $id\_resI$i.pdb.000 $id\_ins$k.pdb`;
            if ($args{v})
            {
                print `/home/layla/data/SCWRL4/Scwrl4 -i $id\_resI$i.pdb.000 -o $id\_ins$k.pdb -s $id\_resI$i.scwrl`;
            }
            else
            {
                `/home/layla/data/SCWRL4/Scwrl4 -i $id\_resI$i.pdb.000 -o $id\_ins$k.pdb -s $id\_resI$i.scwrl`;
            }
        }


	else
	{
	    print "Warning: No model produced by loboAuto.\n";
	}
}
    }
    
    $numins = int ($#ins / 2) + 1;
    print `cp $id\_ins$numins.pdb $id\_insall.pdb`;
    print '*-' x 30, "*\n";
}


sub printContent
{
    print "INDEL file content:\n", '-' x 19, "\n";
    for ($i = 0; $i <= $#head; $i += 2)
    {
	print "HEAD:: \t$head[$i]\t$head[$i+1]\n";
    }
    print '-' x 19, "\n";
    for ($i = 0; $i <= $#del; $i += 2)
    {
	print "DEL :: \t$del[$i]\t$del[$i+1]\n";
    }
    print '-' x 19, "\n";
    for ($i = 0; $i <= $#ins; $i += 2)
    {
	print "INS :: \t$ins[$i]\t$ins[$i+1]\n";
    }
    print '-' x 19, "\n";
    for ($i = 0; $i <= $#tail; $i += 2)
    {
	print "TAIL:: \t$tail[$i]\t$tail[$i+1]\n";
    }
    print '*-' x 30, "*\n";

}

sub checkInput
{
    die "Error: No input INDEL file specified. (-h for help)\n"
	unless ($args{I});
    die "Error: No output file specified. (-h for help)\n"
	unless ($args{o});
    die "Error: No input PDB file specified. (-h for help)\n"
	unless ($args{i});
    die "Error: No input sequence file specified. (-h for help)\n"
	unless ($args{s});
    die "Error: Invalid input PDB file.\n" 
	unless open (A,"$args{i}");
    close A;
    die "Error: Invalid input sequence file.\n" 
	unless open (B,"$args{s}");
    close B;
}

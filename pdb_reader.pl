#!/usr/bin/perl

# This program reads a pdb file and saves the
# 3-D coordinates of the protein.

$ProteinInputFile = '2erl_pdb.txt';

open(PDBFILE, $ProteinInputFile);

$x = 0;
do {
	$atom[$x] = <PDBFILE>;
	chomp($atom[$x]);
	$atomname[$x] = substr($atom[$x], 0 , 5);
#	print $atomname[$x];
#	print "\n";
} until ($atomname[$x] =~ 'ATOM ');

do {
	if ($atomname[$x] =~ 'ATOM ') {
	
	if ((substr($atom[$x], 8, 1) != ' ')){
		$atomnum[$x] = (substr($atom[$x], 8, 3));
	}
	elsif ((substr($atom[$x], 9, 1) != ' ')){
		$atomnum[$x] = (substr($atom[$x], 9, 2));
	}
	else {
		$atomnum[$x] = (substr($atom[$x], 10, 1));
	}
	

#	Here we read in the residue number
	if ((substr($atom[$x], 24, 1) != ' ')) {
		$resnum[$x] = substr($atom[$x], 24, 2);
	}
	else {
		$resnum[$x] = substr($atom[$x], 25, 1);
	}

#	Here we read in the name of the atom (carbon, nitrogen, etc.)
	$atomtype[$x] = substr($atom[$x], 13, 3);
	chomp ($atomtype[$x]);


#	Here we read in what we call the x-coordinate
	if ((substr($atom[$x], 32, 1) != ' ')) {
		$xcoor[$x] = substr($atom[$x], 32, 6);
	}
	else {
		$xcoor[$x] = substr($atom[$x], 33, 5);
	}


#	Here we read in what we call the y-coordinate
	if (not(substr($atom[$x], 40, 1) =~ ' ')) {
		$ycoor[$x] = substr($atom[$x], 40, 6);
	}
	else {
		$ycoor[$x] = substr($atom[$x], 41, 5);
	}
	if ($x == 2 || $x == 7 || $x == 102) {
		$tmp = substr($atom[$x], 40, 1);
		print "temp is: $tmp\n";
	}


#	Here we read in what we call the z-coordinate
	if ((substr($atom[$x], 48, 1) != ' ')) {
		$zcoor[$x] = substr($atom[$x], 48, 6);
	}
	else {
		$zcoor[$x] = substr($atom[$x], 49, 5);
	}


	}

	$x = $x + 1;
	$atom[$x] = <PDBFILE>;
	chomp ($atom[$x]);
	$atomname[$x] = substr($atom[$x], 0 , 6);
} until (($atomname[$x] !~ 'ATOM  ') && ($atomname[$x] !~ 'ANISOU'));


close PDBFILE;

$tot = 0;
# We save the coordinate information to another array.
for ($count = 0; $count <= $x; $count++) {
	if ($atomtype[$count] =~ 'CA') {
		$xc[$tot] = $xcoor[$count];
		$yc[$tot] = $ycoor[$count];
		$zc[$tot] = $zcoor[$count];
		$tot++;
	}
}

# Atomname and number of the last residue.
print "atomname is: $atomname[$x]\n";
print "tot is:  $tot\n";


$DataOutputFile = "2erl_contact.txt";

print "Writing to file...\n";

open(CONFILE, ">$DataOutputFile");

print CONFILE "$tot\n";

$totx = 0;
$toty = 0;
$totz = 0;
$ave = 0;

# We add up the x, y, and z coordinates in order to find the average values.
for ($c = 0; $c < $tot; $c++) {
	print "x: $xc[$c], y: $yc[$c], z: $zc[$c]\n";
	$totx = $totx + $xc[$c];
	$toty = $toty + $yc[$c];
	$totz = $totz + $zc[$c];
	print "total x: $totx\ttotal y: $toty\ttotal z: $totz\n";
}

# Finding the mass center of the structure.
$avex = $totx / $tot, print "average x is: $avex\n";
$avey = $toty / $tot, print "average y is: $avey\n";
$avez = $totz / $tot, print "average z is: $avez\n";

$rsq = 0;

for ($c = 0; $c < $tot; $c++) {
	$rsq = $rsq + (($xc[$c] - $avex)**2 + ($yc[$c] - $avey)**2 + ($zc[$c] - $avez)**2);
}

	$rsq = $rsq / $tot;
	print "r_squared for this protein is: $rsq\n";

# We calculate a contact graph for the structure based on a distance measure between
# residues.
for ($c = 0; $c < $tot; $c++) {
	for ($d = $c+1; $d < $tot; $d++) {
		$diffx = $xc[$c] - $xc[$d];
		$diffy = $yc[$c] - $yc[$d];
		$diffz = $zc[$c] - $zc[$d];
		$dist = ($diffx * $diffx) + ($diffy * $diffy) + ($diffz * $diffz);
		$newc = $c+1;
		$newd = $d+1;
		if (($dist < 36) || ($newd - $newc == 1)) {
			print CONFILE "$newc $newd\n";
		}
	}
}
print CONFILE "-1 -1";

close CONFILE;

exit;

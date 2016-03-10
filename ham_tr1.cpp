#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <algorithm>
#include <stdarg.h>
#include <time.h>
#include <ctype.h>
#include <math.h>

#define SETDIM 30000


#include <time.h>                      // define time()
#include "randomc.h"                   // define classes for random number generators
#include "userintf.cpp"                // define system specific user interface

// include code for one of the random number generators:
#include "mersenne.cpp"                // members of class TRandomMersenne

// printout is used to print the conformations with a single sequence folding to them with lowest energy.
// The conformation along with the sequence of the model peptide is printed out.

void printout(int*** mmem, int iii1[2], int seq[21], double ggtot, double etot, double etot_old, int filenum)
{
	char ppath[199][199][2];
	int aaa, a1[2], aa1, aa2, aa3, hh1;
	char bbb[2], b1[2];
	FILE *out1;


	if (filenum == 1)
	{
		out1 = fopen("ham_tr1_new.txt", "a");
	}
	else if (filenum == 2)
	{
		out1 = fopen("ham_tr2_new.txt", "a");
	}

	fprintf(out1, "iii1[0] is: %d\n", iii1[0]);


//	This for loop prints out the conformation
	for (int x = 0; x < ggtot; x++)
	{
		for (int l2 = 1; l2 <= 6; l2++)
		{
			for (int l1 = 1; l1 <= 11; l1++)
			{
				if ((l1 < 11 || l2 < 6) && (l2 % 1 == 1))
				{
					ppath[2*l1][2*l2][0] = ' ';
					ppath[2*l1][2*l2][1] = 0;
				}
				else
				{
					ppath[2*l1][2*l2][0] = 0;
					ppath[2*l1][2*l2][1] = 0;
				}
				aaa = mmem[iii1[x]][l1][l2];
				a1[1] = aaa/10;
				a1[2] = aaa - a1[1]*10;
				for (int ia = 1; ia <= 2; ia++)
				{
					if(a1[ia] == 0) b1[ia] = '0';
					if(a1[ia] == 1) b1[ia] = '1';
					if(a1[ia] == 2) b1[ia] = '2';
					if(a1[ia] == 3) b1[ia] = '3';
					if(a1[ia] == 4) b1[ia] = '4';
					if(a1[ia] == 5) b1[ia] = '5';
					if(a1[ia] == 6) b1[ia] = '6';
					if(a1[ia] == 7) b1[ia] = '7';
					if(a1[ia] == 8) b1[ia] = '8';
					if(a1[ia] == 9) b1[ia] = '9';
				}


				ppath[2*l1-1][2*l2-1][0] = b1[1];
				ppath[2*l1-1][2*l2-1][1] = b1[2];

				aa1 = mmem[iii1[x]][l1+2][l2];
				if(abs(aa1-aaa) == 1 && aaa != 0 && aa1 != 0)
				{
					ppath[2*l1][2*l2-1][0] = '-';
					ppath[2*l1][2*l2-1][1] = '-';
				}
				else
				{
					ppath[2*l1][2*l2-1][0] = 0;
					ppath[2*l1][2*l2-1][1] = 0;
				}


				aa2 = mmem[iii1[x]][l1-1][l2+1];
				if (abs(aa2-aaa) == 1 && aaa != 0)
				{
					ppath[2*l1-2][2*l2][0] = '/';
					ppath[2*l1-1][2*l2][1] = ' ';
				}
				else if (ppath[2*l1-1][2*l2][0] != '\\')
				{
					ppath[2*l1-1][2*l2][0] = 0;
					ppath[2*l1-1][2*l2][1] = 0;
				}
				aa3 = mmem[iii1[x]][l1+1][l2+1];
				if (abs(aa3-aaa) == 1 && aaa != 0)
				{
					ppath[2*l1+1][2*l2][0] = '\\';
					ppath[2*l1+1][2*l2][1] = 0;
				}
				else
				{
					ppath[2*l1+1][2*l2][0] = 0;
					ppath[2*l1+1][2*l2][1] = 0;
				}

				printf("");
			}
		}

		for (int j1 = 1; j1 <= 11; j1++)
		{
			if (j1 % 2 == 1)
			{
				hh1 = ((11 - j1) /2);
				for(int h1 = 0; h1 < hh1; h1++)
				{
					fprintf(out1, " ");
				}
			}
			for (int i1 = 1; i1 <= 21; i1++)
			{
				if (j1 % 2 == 0 && (ppath[i1][j1][0] == '0' && ppath[i1][j1][1] == '0'))
				{
				//	fprintf(out1, "%c", ppath[i1][j1][0]);
				}
				else if(ppath[i1][j1][0] != '0' || ppath[i1][j1][1] != '0')
				{
					if (ppath[i1-2][j1][0] == '-')
					{
						fprintf(out1, "--");
					}
					else if (j1 % 2 == 1)
					{
						fprintf(out1, "%c%c", ppath[i1][j1][0], ppath[i1][j1][1]);
					}
					else
					{
						fprintf(out1, "%c", ppath[i1][j1][0]);
						if (i1 % 2 == 1)
							fprintf(out1, " ");
					}
				}
			}
			fprintf(out1, "\n");
		}
	}

//	Output the sequence itself.

	for (int v = 0; v < 21; v++)
	{
		if (seq[v] == 1)
			fprintf(out1, "H");
		else
			fprintf(out1, "P");
	}
	fprintf(out1, "\n");
	fprintf(out1, "There were %g conformations of %g energy for that sequence\n", ggtot, etot);
	fprintf(out1, "%g %g\n", etot, etot_old);

	fclose (out1);
	

}

int main()
{




//	program square
//
//	The program calculates all possible hamiltonian paths
//	on a triangular lattice.  The path need not be closed.
//	Head-tail symmetries are eliminated
//


	int x[100], y[100], nnab[100], nab[8][100], res[100];
	int iseq[100], nseq[100], kseq[100], map[101][101], map1[101][101];
	int elig[100], inuse[100], xsymm[100], ysymm[100];
	int xysymm[100], xyasymm[100], xjoint[100], yjoint[100];
	int tabmax, rot180[100], mx[100], my[100];
	int rot90[100], rotm90[100], mdiag[100], madiag[100];
	int ncnt, excl[1200][25], aa, a[2], a1, a2, inv[100], symm[100][25];
	char path[199][199][2], bb[2], b[2], p1, key;
	int nexcl, nsymm, m, n, msites, xmin, ymin, xmax, ymax;
	int mhalf, nhalf, ix, iy, jx, jy, istart, k, mm, ip;
	int n180, n90, nm90, nmx, nmy, ndiag, nadiag, in, s;
	int ndis, i12, nj, jj;
	int vv, qq;
	int flag, gntot, tmp;
	float tot, gtot, gtot_old, ctot, htot;
	float tot2, gtot2, gtot2_old, ctot2, htot2;
	int aaa[50];
	int mflag, ii2[2], ii3[2], ii1;
	int count[SETDIM], count2[SETDIM];
	FILE *out1;
	double r;

	int32 seed = time(0);                // random seed

	TRandomMersenne rg(seed);

	htot = 0;
	gtot = 0;
	ctot = 0;
	mflag = 0;

//	mem will hold all possible conformations.
	int*** mem = new int** [SETDIM];
	for (int i = 0; i < SETDIM; i++)
	{
		mem[i] =  new int* [40];
		for (int j = 0; j < 40; j++)
		{
			mem[i][j] = new int [40];
		}
		if (i % 5000 == 0)
			printf ("i is: %d\n", i);
	}

//	aaa will hold the binary aa sequence.
	for (int i = 0; i < 22; i++)
	{
		aaa[i] = (i % 2);
	}


	tabmax=100;

//	the output files.
	out1 = fopen("ham_tr1_new.txt", "w");
	fclose (out1);

	out1 = fopen("ham_tr2_new.txt", "w");
	fclose (out1);

	for (int i = 0; i < SETDIM; i++)
	{
		count[i] = 0;
		count2[i] = 0;
	}


//
//	x(i) and y(i) are the coordinates of the ith point.
//	nnab(i) is the number of neighbors of the point,
//	nab(j,i) is the number in the table of the j'th neighbor
//	for j <= nnab(i).
//	elig(i) = .true. if the i'th point can be an end point of the path
//	res(i)= is the residue mod 2
//	of the sum of the coordinates of the i'th point.
//	xsymm(i) = .true.  if point i  is on one of the symmetry 
//	axis perpendicular to the x-axis.
//	ysymm(i) = .true.  if point i  is on one of the symmetry 
//	axis perpendicular to the y-axis.
//	xysymm(i) and xyasymm(i) correspond to point i on the
//	diagonal and antidiagonal of the square
//

	nexcl=0;
	ncnt=0;
	nsymm=0;

//	printf ("do you want to print paths? enter y or n:");
//	cin >> p1;
//	cin.get(p1);


//	scanf("%c", &p1);


	msites=21;


	for (int j = 1; j<=tabmax; j++)
	{
		for (int i = 1; i <=tabmax; i++)
		{
			map1[i][j]=0;
			map[i][j]=0;
		}
	}

	xmin=tabmax+1;
	ymin=tabmax+1;
	xmax=0;
	ymax=0;

//	initiating the variables
	for (int i = 1; i<=msites; i++)
	{
		inuse[i] = 0;
		elig[i] = 0;
		nnab[i] = 0;
		for (int j = 1; j<=4; j++)
		{
			nab[j][i] = 0;
		}
	}


//	map is an array used to hold information about where the
//	residues are in space. 
	map[6][1] = 1;

	map[5][2] = 2;
	map[7][2] = 3;

	map[4][3] = 4;
	map[6][3] = 5;
	map[8][3] = 6;

	map[3][4] = 7;
	map[5][4] = 8;
	map[7][4] = 9;
	map[9][4] = 10;

	map[2][5] = 11;
	map[4][5] = 12;
	map[6][5] = 13;
	map[8][5] = 14;
	map[10][5] = 15;

	map[1][6] = 16;
	map[3][6] = 17;
	map[5][6] = 18;
	map[7][6] = 19;
	map[9][6] = 20;
	map[11][6] = 21;


//	nab holds the neighbors for each residue.

	nab[1][1] = 2;
	nab[2][1] = 3;

	nab[1][2] = 1;
	nab[2][2] = 3;
	nab[3][2] = 4;
	nab[4][2] = 5;

	nab[1][3] = 1;
	nab[2][3] = 2;
	nab[3][3] = 5;
	nab[4][3] = 6;

	nab[1][4] = 2;
	nab[2][4] = 5;
	nab[3][4] = 7;
	nab[4][4] = 8;

	nab[1][5] = 2;
	nab[2][5] = 3;
	nab[3][5] = 4;
	nab[4][5] = 6;
	nab[5][5] = 8;
	nab[6][5] = 9;

	nab[1][6] = 3;
	nab[2][6] = 5;
	nab[3][6] = 9;
	nab[4][6] = 10;

	nab[1][7] = 4;
	nab[2][7] = 8;
	nab[3][7] = 11;
	nab[4][7] = 12;

	nab[1][8] = 4;
	nab[2][8] = 5;
	nab[3][8] = 7;
	nab[4][8] = 9;
	nab[5][8] = 12;
	nab[6][8] = 13;

	nab[1][9] = 5;
	nab[2][9] = 6;
	nab[3][9] = 8;
	nab[4][9] = 10;
	nab[5][9] = 13;
	nab[6][9] = 14;

	nab[1][10] = 6;
	nab[2][10] = 9;
	nab[3][10] = 14;
	nab[4][10] = 15;

	nab[1][11] = 7;
	nab[2][11] = 12;
	nab[3][11] = 16;
	nab[4][11] = 17;

	nab[1][12] = 7;
	nab[2][12] = 8;
	nab[3][12] = 11;
	nab[4][12] = 13;
	nab[5][12] = 17;
	nab[6][12] = 18;

	nab[1][13] = 8;
	nab[2][13] = 9;
	nab[3][13] = 12;
	nab[4][13] = 14;
	nab[5][13] = 18;
	nab[6][13] = 19;


	nab[1][14] = 9;
	nab[2][14] = 10;
	nab[3][14] = 13;
	nab[4][14] = 15;
	nab[5][14] = 19;
	nab[6][14] = 20;

	nab[1][15] = 10;
	nab[2][15] = 14;
	nab[3][15] = 20;
	nab[4][15] = 21;

	nab[1][16] = 11;
	nab[2][16] = 17;

	nab[1][17] = 11;
	nab[2][17] = 12;
	nab[3][17] = 16;
	nab[4][17] = 18;

	nab[1][18] = 12;
	nab[2][18] = 13;
	nab[3][18] = 17;
	nab[4][18] = 19;

	nab[1][19] = 13;
	nab[2][19] = 14;
	nab[3][19] = 18;
	nab[4][19] = 20;

	nab[1][20] = 14;
	nab[2][20] = 15;
	nab[3][20] = 19;
	nab[4][20] = 21;

	nab[1][21] = 15;
	nab[2][21] = 20;

//	msites is the total number of residues for a model.

	msites = 21;

//	nnab holds the number of neighbors of a
//	residue on a lattice.

	nnab[1] = 2;
	nnab[2] = 4;
	nnab[3] = 4;
	nnab[4] = 4;
	nnab[5] = 6;
	nnab[6] = 4;
	nnab[7] = 4;
	nnab[8] = 6;
	nnab[9] = 6;
	nnab[10] = 4;
	nnab[11] = 4;
	nnab[12] = 6;
	nnab[13] = 6;
	nnab[14] = 6;
	nnab[15] = 4;
	nnab[16] = 2;
	nnab[17] = 4;
	nnab[18] = 4;
	nnab[19] = 4;
	nnab[20] = 4;
	nnab[21] = 2;


	elig[1] = 1;
	elig[2] = 1;
	elig[4] = 1;
	elig[5] = 1;
	elig[8] = 1;


//	x and y hold the coordinates of each position on the lattice. They
//	are used to calculate energy values.

	x[1] = 6;
	x[2] = 5;
	x[3] = 7;
	x[4] = 4;
	x[5] = 6;
	x[6] = 8;
	x[7] = 3;
	x[8] = 5;
	x[9] = 7;
	x[10] = 9;
	x[11] = 2;
	x[12] = 4;
	x[13] = 6;
	x[14] = 8;
	x[15] = 10;
	x[16] = 1;
	x[17] = 3;
	x[18] = 5;
	x[19] = 7;
	x[20] = 9;
	x[21] = 11;

	y[1] = 1;
	y[2] = 2;
	y[3] = 2;
	y[4] = 3;
	y[5] = 3;
	y[6] = 3;
	y[7] = 4;
	y[8] = 4;
	y[9] = 4;
	y[10] = 4;
	y[11] = 5;
	y[12] = 5;
	y[13] = 5;
	y[14] = 5;
	y[15] = 5;
	y[16] = 6;
	y[17] = 6;
	y[18] = 6;
	y[19] = 6;
	y[20] = 6;
	y[21] = 6;


//	initialization of inuse - right now no positions are used.		
	for (int ist = 1; ist <= msites; ist++)
	{
		inuse[ist] = 0;
	}

	flag = 0;


	for (int istart = 1; istart <= msites; istart++)
	{
//		only a few positions will be used to start the polypeptide chain.
		if (!elig[istart]) goto L888;
		iseq[1] = istart;

//
//	This whole block of code is needed to eliminate symmetrical conformations.
//
		if (iseq[1] == 1)
		{
			nab[1][1] = 2;
			nnab[1] = 1;
		}
		else
		{
			nab[1][1] = 2;
			nab[2][1] = 3;
			nnab[1] = 2;
		}
		if (iseq[1] == 5)
		{
			nab[1][5] = 2;
			nab[2][5] = 4;
			nab[3][5] = 8;
			nnab[5] = 3;
		}
		else
		{
			nab[1][5] = 2;
			nab[2][5] = 3;
			nab[3][5] = 4;
			nab[4][5] = 6;
			nab[5][5] = 8;
			nab[6][5] = 9;
			nnab[5] = 6;
		}

		if (iseq[1] == 8)
		{
			nab[1][8] = 4;
			nab[2][8] = 5;
			nab[3][8] = 9;
			nnab[8] = 3;
		}
		else
		{
			nab[1][8] = 4;
			nab[2][8] = 5;
			nab[3][8] = 7;
			nab[4][8] = 9;
			nab[5][8] = 12;
			nab[6][8] = 13;
			nnab[8] = 6;
		}

		inuse[istart] = 1;
		s = 1;
		ip = istart;





//
//	iseq(i) is the index number of the ith point in the path
//	nseq(i) is the number of neighbors of the ith point in the path.
//	kseq(i) <= nseq(i) indicates which neighbor we are currently considering
//	for (i+1)st point in the path.
//
//	we advance one point in our sequencing
//
L280:

//	first a test to see if we've reached the end of the peptide.
		if (s == msites) 
		{
			goto L400;
		}

//	we assign a position to our growing chain
		nseq[s] = nnab[iseq[s]];
		kseq[s] = 0;

//
//	we try the next neighbor
//
L360:


		kseq[s]++;
		k = kseq[s];


		if (k > nseq[s]) goto L380;

		mm = nab[k][ip];



L359:
//	A test to see if the position is already assigned.
		if(inuse[mm]) goto L360;
		inuse[mm] = 1;
		s++;
		iseq[s] = mm;
		ip = mm;
		goto L280;

//
//	we back up a step in our sequencing
//
L380:
		inuse[iseq[s]] = 0;
		s--;
		if (s == 0) goto L888;
		ip = iseq[s];
		goto L360;

//
//	we have completed a sequence, and print it out
//



L400:

//	ncnt holds the total number of conformations so far.
			
		ncnt++;


//	this is to keep track of how many conformations we are up to.
		if ((ncnt % 5000) == 0) printf ("ncnt is: %d\n", ncnt);

//	we assign the conformation to mem and keep it for later use.
		for (int i1 = 1; i1 <= msites; i1++)
		{
			mem[ncnt][x[iseq[i1]]][y[iseq[i1]]] = i1;
		}




		goto L360;
L888:
		printf("");
	}

//	The conformations are generated, we print out how many there are.

	printf ("ncnt is: %d\n", ncnt);

//
//	In the previous part of the program we generated all conformations
//	for our lattice. We now thread binary sequences onto all those conformations
//	and calculate an energy value for each threading. Sequences that fold to a
//	unique conformation with lowest energy are kept and printed out.
//

//	We can either go through each sequence exhaustively or we can generate sequences randomly.
	for (aaa[0] = 0; aaa[0] < 2; aaa[0]++)
	{
	for (aaa[1] = 0; aaa[1] < 2; aaa[1]++)
	{
	for (aaa[2] = 0; aaa[2] < 2; aaa[2]++)
	{
	for (aaa[3] = 0; aaa[3] < 2; aaa[3]++)
	{
	for (aaa[4] = 0; aaa[4] < 2; aaa[4]++)
	{
	for (aaa[5] = 0; aaa[5] < 2; aaa[5]++)
	{
	for (aaa[6] = 0; aaa[6] < 2; aaa[6]++)
	{
	for (aaa[7] = 0; aaa[7] < 2; aaa[7]++)
	{
	for (aaa[8] = 0; aaa[8] < 2; aaa[8]++)
	{
	for (aaa[9] = 0; aaa[9] < 2; aaa[9]++)
	{
	for (aaa[10] = 0; aaa[10] < 2; aaa[10]++)
	{
	for (aaa[11] = 0; aaa[11] < 2; aaa[11]++)
	{
	for (aaa[12] = 0; aaa[12] < 2; aaa[12]++)
	{
	for (aaa[13] = 0; aaa[13] < 2; aaa[13]++)
	{
	for (aaa[14] = 0; aaa[14] < 2; aaa[14]++)
	{
	for (aaa[15] = 0; aaa[15] < 2; aaa[15]++)
	{
	for (aaa[16] = 0; aaa[16] < 2; aaa[16]++)
	{
	for (aaa[17] = 0; aaa[17] < 2; aaa[17]++)
	{
	for (aaa[18] = 0; aaa[18] < 2; aaa[18]++)
	{
	for (aaa[19] = 0; aaa[19] < 2; aaa[19]++)
	{
	for (aaa[20] = 0; aaa[20] < 2; aaa[20]++)


//	For each sequence, we test all conformations. The following block is used
//	to generate a threading energy for a conformation and a sequence.
		for (ii1 = 1; ii1 <= ncnt; ii1++)
		{
			for (int tt = 2; tt < 7; tt++)
			{
				tmp = 6 - (tt - 1);
				for (int uu = 1; uu <= tt - 1; uu++)
				{
					if((abs(mem[ii1][tmp + 2 * (uu - 1)][tt] - mem[ii1][tmp + 2 * (uu)][tt]) != 1)
					&& (aaa[mem[ii1][tmp + 2 * (uu - 1)][tt] - 1] == 1 && aaa[mem[ii1][tmp + 2 * (uu)][tt] - 1] == 1))
					{
						tot++;
						tot2 = tot2 + 2.3;
					}
					if((abs(mem[ii1][tmp + 2 * (uu - 1)][tt] - mem[ii1][tmp + 2 * (uu)][tt]) != 1)
					&& (aaa[mem[ii1][tmp + 2 * (uu - 1)][tt] - 1] + aaa[mem[ii1][tmp + 2 * (uu)][tt] - 1] == 1))
					{
						tot2++;
					}
				}
			}

			for (int tt = 1; tt < 6; tt++)
			{
				tmp = 5 + tt;
				for (int uu = 1; uu < 7 - tt; uu++) 
				{
					if((abs(mem[ii1][tmp - (uu - 1)][uu + (tt - 1)] - mem[ii1][tmp - (uu - 2)][uu + 1 + (tt - 1)]) != 1)
					&& (aaa[mem[ii1][tmp - (uu - 1)][uu + (tt - 1)] - 1] == 1 && aaa[mem[ii1][tmp - (uu - 2)][uu + 1 + (tt - 1)] - 1] == 1))
					{
						tot++;
						tot2 = tot2 + 2.3;
					}
					if((abs(mem[ii1][tmp - (uu - 1)][uu + (tt - 1)] - mem[ii1][tmp - (uu - 2)][uu + 1 + (tt - 1)]) != 1)
					&& (aaa[mem[ii1][tmp - (uu - 1)][uu + (tt - 1)] - 1] + aaa[mem[ii1][tmp - (uu - 2)][uu + 1 + (tt - 1)] - 1] == 1))
					{
						tot2++;
					}
				}
			}

			for (int tt = 1; tt < 6; tt++)
			{
				tmp = 5 + tt;
				for (int uu = 1; uu < 7 - tt; uu++)
				{
					if((abs(mem[ii1][tmp - (uu - 1)][uu + (tt - 1)] - mem[ii1][tmp - uu][uu + 1 + (tt - 1)]) != 1)
					&& (aaa[mem[ii1][tmp - (uu - 1)][uu + (tt - 1)] - 1] == 1 && aaa[mem[ii1][tmp - uu][uu + 1 + (tt - 1)] - 1] == 1))
					{
						tot++;
						tot2 = tot2 + 2.3;
					}
					if((abs(mem[ii1][tmp - (uu - 1)][uu + (tt - 1)] - mem[ii1][tmp - uu][uu + 1 + (tt - 1)]) != 1)
					&& (aaa[mem[ii1][tmp - (uu - 1)][uu + (tt - 1)] - 1] + aaa[mem[ii1][tmp - uu][uu + 1 + (tt - 1)] - 1] == 1))
					{
						tot2++;
					}
				}
			}


//	After we finish finding the energy, we test to see if it's the lowest energy so far.

			if (tot == gtot && tot != 0)
			{
				ctot++;
			//	if (ctot == 2)
			//	{
			//		ii2[1] = ii1;
			//	}
			}

			if (tot2 == gtot2 && tot2 != 0)
			{
				ctot2++;
			//	if (ctot2 = 2)
			//	{
			//		ii3[1] = ii1;
			//	}
			}
		
//	We hold onto the value if it is lowest in energy.
			if (gtot < tot)
			{
				ii2[0] = ii1;
				gtot_old = gtot;
				gtot = tot;
				ctot = 1;
				flag++;
			}

			gtot2 = float(int(gtot2 * 10 + 0.5))/10;
			tot2 = float(int(tot2 * 10 + 0.5))/10;
			if (gtot2 < tot2)
			{
				gtot2_old = gtot2;
				gtot2 = tot2;
				ctot2 = 1;
				ii3[0] = ii1;
			}

			if (gtot_old < tot && tot < gtot)
			{
				gtot_old = tot;
			}

			if (gtot2_old < tot2 && tot2 < gtot2)
			{
				gtot2_old = tot2;
			}

L98:
			printf("");


			tot = 0;
			tot2 = 0;


L361:
			printf("");

L800:
			printf("");	
		}

L900:

//	We've finished threading one sequence onto all conformations.
//	If there is one conformation with lowest energy, we print out
// 	both the conformation and the sequence.
		if (ctot  == 1)
		{
			count[ii2[0]]++;
			printout(mem, ii2, aaa, ctot, gtot, gtot_old, 1);
		}

		if (ctot2  == 1)
		{
			count2[ii3[0]]++;
			printout(mem, ii3, aaa, ctot2, gtot2, gtot2_old, 2);
		}

		flag = 0;
		htot = 0;
		ctot = 0;
		gtot = 0;

		htot2 = 0;
		ctot2 = 0;
		gtot2 = 0;

	}}}}}}}}}}}}}}}}}}}}}

       goto L999;
//
// failure messages printed out
L920:
	printf("An error somewhere!\n");
L999:


	for (int i = 0; i < SETDIM; i++)
	{
		for (int j = 0; j < 40; j++)
		{
			delete[] mem[i][j];
		}
		delete[] mem[i];
	}
	delete[] mem;


	for (int i = 0; i < SETDIM; i++)
	{
//		delete[] set[i];
//		delete[] rot[i];
	}
//	delete[] set;
//	delete[] rot;

}



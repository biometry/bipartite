#include <R.h>
// #include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
using namespace std;

extern "C" {
//#include "IOfiles.h"



/*	Input and output files */


/* caracter LF:fin de linea, CE:cero y UN:uno */
#define LF 10
#define CE 48
#define UN 49

//	Para Numerical Recipes
#define NR_END 1
#define FREE_ARG char*

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 3.e-7
#define RNMX (1.0-EPS)

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M_IND 7
#define NSTACK 50
   
#define ITMAX 100
#define EPSB 3.0e-8
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//	El valor más alto de U econtrado es:
const double UMAX = 0.04145;

//	para el algoritmo genetico:
unsigned long NBRS;
int POPSIZE;	// Must be at least 15
int TOURSIZE;
int NBGENER;
int OVER;
int *toc;
int *tor;
int bmo;


string INFILE;
string OUTFILE;
FILE *out;

/*
 ***************************************************************** 
 */

void matrixSize(string inputFile, int& nrows, int& ncols, int& skip);

void endnote(FILE *f);

void readMatrix(string inputFile, int nrows, int ncols, int skip, int **m);

double matrixTemperature(bool &success,int sp, int **mat, int *c_ord, int *r_ord, int nr, int nc, long& idum);

void orderMatrix(int **mat, int *c_ord, int *r_ord, int nrows, int ncols, int& effInsect, int& effPlant);

void removeBlacks(int **mat, int *c_ord, int *r_ord, int effInsect, int effPlant, int& nbRows, 
				  int& nbCols, double& phi);

void calcZ(double phi, double& z);

void calcDistance(double z, double border[],double **mat,int nr,int nc);

double packMatrix(int sp, int **dataMat, int **intMat, double **d, int *c_ord, int *r_ord, int nr, int nc, 
				  int nrows, int ncols, long& idum);

double calcTemp(double **d, int **mat, int indr[], int indc[], int nr, int nc);

void calcIdiosyncTemp(double **d, int **mat, int indr[], int indc[], int nr, int nc);

void prePackMatrix(int **mat, int indr[], int indc[], int nr, int nc, double x);

void prePackrows(int **mat, int indr[], int indc[], int nr, int nc, double x);

void prePackcols(int **mat, int indr[], int indc[], int nr, int nc, double x);

void prePackNTC(int **mat, int indr[], int indc[], int nr, int nc);

void prePackNTCrows(int **mat, int indr[], int indc[], int nr, int nc);

void prePackNTCcols(int **mat, int indr[], int indc[], int nr, int nc);

void permute(long& idum, int n, int index[]);

void mutate(long& idum, int n, int index[]);

void choosePlayers(long& idum, int n,int m,int arr[]);

void crossOver(long& idum, int nr, int nc, int pr[],
              int pc[],
              int or_new[],
              int oc[])
              ;

void indexx(int n, int arr[], int indx[]);

void indexxD(int n, double arr[], int indx[]);

double zbrent(int nr, int nc, double u1, double v1, double z);

double func(double yy, int nr, int nc, double x1, double y1, double z);

double ran1(long& idum);

void avevar(double data[], unsigned long n, double& ave, double& var);

int *ivector(long nl, long nh);

void free_ivector(int *v, long nl, long nh);

double *vector(long nl, long nh);

void free_vector(double *v, long nl, long nh);

int **imatrix(long nrl, long nrh, long ncl, long nch);

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);

double **matrix(long nrl, long nrh, long ncl, long nch);

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);

//void Rf_error(char error_text[]);

/*
 ***************************************************************** 
 */

int bmn5(int *matrix, int *n_rows, int *n_cols, double *temperature,
        int *n_nullmodels, int *p_size, int *n_individuals, int *binmatout,
        int *n_generations, int *nullmodels,
        double *p_null1, double *avt_null1, double 	*avv_null1,
        double *p_null2, double *avt_null2, double 	*avv_null2,
        double *p_null3, double *avt_null3, double 	*avv_null3,
        int *poc, int *por
        )
{
	bool success;
	int ncols,nrows;//,skip;
	int i,j;
	int nPresI,nPresJ=0;
	int *colOrder,*rowOrder;
	int **interactMat;
	int **tempMat;
	unsigned long nr;
	long idum=-48367;
	//idum = time(0);
	double fill;
	double p,p1,p2,p3,r;
	double tInt,tTemp;
	double ave,var;
	double *pCol,*pRow;
	double *nullTemp;
	
	toc=poc;     // move pointer of column and row order
  tor=por;     // to R objects in bmn
  bmo = *binmatout;

	

/*
	int stoch;
	fprintf(stdout, "This program can run in two modes: stochastic and deterministic.\n");
	fprintf(stdout, "In the deterministic mode, if you analyse the same matrix several times\n");
	fprintf(stdout, "it will always give the same results. In stochastic mode, \n");
	fprintf(stdout, "it will find slightly different results each time you analyse the data.\n");
	fprintf(stdout, "Type 1 if you want to use the stochastic mode, \n");
	fprintf(stdout, "and type any other number if you want to use the deterministic version.\n");
	cin >> stoch;
	fprintf(stdout, "\n\n");

	if (stoch==1)
	{
		idum = time(0);			// auto seeding using system time
		idum *= -1;
	}
*/
/*	fprintf(stdout, "Type the name of the input data file,\n");
	fprintf(stdout, "then press ENTER\n");
*/	//getline(cin,INFILE,'\n');
//	INFILE = "web.txt";

/*	fprintf(stdout, "\n\nType the name of the file\n");
	fprintf(stdout, "where the results should be saved\n");
	fprintf(stdout, "(do not forget the txt extension),\n");
*///	fprintf(stdout, "then press ENTER\n");
	//getline(cin,OUTFILE,'\n');
	


/*	do {
		fprintf(stdout, "\n\nType the number of null matrices\n");
		fprintf(stdout, "that you want to use to calculate p values,\n");
		fprintf(stdout, "then press ENTER\n");
		cin >> NBRS;
		if (NBRS<=0) fprintf(stdout, "must be a positive integer\n");
	} while (NBRS<=0);

	do {
		fprintf(stdout, "\n\nType the population size\n");
		fprintf(stdout, "that you want to use in the GA\n");
		fprintf(stdout, "to calculate matrix temperature,\n");
		fprintf(stdout, "then press ENTER.\n");
		fprintf(stdout, "We recommend 30\n");
		cin >> POPSIZE;
		if (POPSIZE<15) fprintf(stdout, "must be at least 15\n");
	} while (POPSIZE<15);

	do {
		fprintf(stdout, "\n\nType the number of individuals\n");
		fprintf(stdout, "that you want to select\n");
		fprintf(stdout, "at each generation of the GA,\n");
		fprintf(stdout, "then press ENTER\n");
		fprintf(stdout, "The number must be at least 4 and smaller than %2i\n",POPSIZE);
		fprintf(stdout, "We recommend 7\n");
		cin >> TOURSIZE;
		if (TOURSIZE<4) fprintf(stdout, "must be at least 4\n");
		if (TOURSIZE>=POPSIZE) fprintf(stdout, "must be smaller than %2i\n",POPSIZE);
	} while (TOURSIZE<4 || TOURSIZE>=POPSIZE);

	do {
		fprintf(stdout, "\n\nType the number of generations\n");
		fprintf(stdout, "for the GA,\n");
		fprintf(stdout, "then press ENTER\n");
		fprintf(stdout, "We recommend 2000\n");
		cin >> NBGENER;
		if (NBGENER<=0) fprintf(stdout, "must be a positive integer\n");
	} while (NBGENER<=0);
*/


nrows=*n_rows;
ncols=*n_cols;
NBRS =*n_nullmodels;
POPSIZE =*p_size;
TOURSIZE =*n_individuals;
NBGENER = *n_generations;

OUTFILE ="binmat.out";

/*	fprintf(stdout, "\n\n\n");
	fprintf(stdout, "BINMATNEST is calculating the temperature of your matrix\n");
	fprintf(stdout, "and the corresponding p-values. Please be patient.\n");
	fprintf(stdout, "\n\n\n");
*/
//	matrixSize(INFILE,nrows,ncols,skip);
	interactMat = imatrix(1,nrows,1,ncols);
	//readMatrix(INFILE,nrows,ncols,skip,interactMat);
//check matrix handing over....


//	out = fopen(OUTFILE.c_str(),"a");
//	fprintf(out,"matrix from R\n\n");



for (int i=1; i<=nrows;i++)
{
for (int j=1; j<=ncols;j++)
{
	interactMat[i][j] = matrix[(j-1)*nrows+(i-1)];
//	fprintf(out,"%1i",interactMat[i][j]);
	}
//	fprintf(out,"\n");
}

//fclose(out);

	//	Some variables are initialised
	fill = 0.0;
	for (i=1;i<=nrows;i++)
	{
		for (j=1;j<=ncols;j++) 
		{
			fill += interactMat[i][j];
		}
	}
	fill /= (1.0*nrows*ncols);

    pRow = ::vector(1,nrows);
    pCol = ::vector(1,ncols);
	for (i=1;i<=nrows;i++)
	{
		pRow[i] = 0;
		for (j=1;j<=ncols;j++) pRow[i] += interactMat[i][j];
		pRow[i] /= 1.0*ncols;
	}
	for (j=1;j<=ncols;j++)
	{
		pCol[j] = 0;
		for (i=1;i<=nrows;i++) pCol[j] += interactMat[i][j];
		pCol[j] /= 1.0*nrows;
	}

	colOrder = ivector(1,ncols);
	rowOrder = ivector(1,nrows);
	tInt = matrixTemperature(success,1,interactMat,colOrder,rowOrder,nrows,ncols,idum);
 if (bmo==1)   { //save binmatnest in output file
	// out = fopen(OUTFILE.c_str(),"a");
	Rprintf("Matrix temperature: T = %11.5f\n",tInt); //deleted: out, before "Matrix
	Rprintf("\n\n");
	Rprintf("Matrix temperature = %11.5f\n",tInt);
	//fclose(out);
}
 	*temperature = tInt;
  // return packing order by pointers (tor and toc via poc and por)
	free_imatrix(interactMat,1,nrows,1,ncols);

if (*nullmodels==1) {  //should nullmodels be calculated??

	tempMat = imatrix(1,nrows,1,ncols);
    nullTemp = ::vector(1,NBRS);
	//	First null model:
	p1 = 0.0;
	for (nr=1;nr<=NBRS;nr++)
	{
		do {
			for (i=1;i<=nrows;i++)
			{
				nPresI = 0;
				while (nPresI==0)
				{
					for (j=1;j<=ncols;j++) 
					{
						p = fill;
						r = ran1(idum);
						tempMat[i][j] = (r<p ? 1 : 0);
						nPresI += tempMat[i][j];
					}
				}
			}
			for (j=1;j<=ncols;j++)
			{
				nPresJ = 0;
				for (i=1;i<=nrows;i++) nPresJ += tempMat[i][j];
				if (nPresJ==0) break;
			}
		} while (nPresJ==0);
		tTemp = matrixTemperature(success,0,tempMat,colOrder,rowOrder,nrows,ncols,idum);
		if (success)
		{
			nullTemp[nr] = tTemp;
			if (tTemp < tInt) p1 += 1.0;
		}
		else nr--;
	}
	p1 /= (1.0*NBRS);
	avevar(nullTemp,NBRS,ave,var);
  //save results of nullmodel 1
	*p_null1 = p1;
	*avt_null1 = ave;
	*avv_null1 = var;
 if (bmo==1)  { //save binmatnest in output file
	//out = fopen(OUTFILE.c_str(),"a");
	Rprintf("  Null model     p-value   Average T    Variance\n"); // was: fprintf
	Rprintf("       First %11.5f %11.5f %11.5f\n",p1,ave,var);
	Rprintf("First null model:  p1 = %11.5f\n",p1);
	//fclose(out);
}
	//	Second null model:
	p2 = 0.0;
	for (nr=1;nr<=NBRS;nr++)
	{
		do {
			for (i=1;i<=nrows;i++)
			{
				nPresI = 0;
				while (nPresI==0)
				{
					for (j=1;j<=ncols;j++) 
					{
						p = pCol[j];
						r = ran1(idum);
						tempMat[i][j] = (r<p ? 1 : 0);
						nPresI += tempMat[i][j];
					}
				}
			}
			for (j=1;j<=ncols;j++)
			{
				nPresJ = 0;
				for (i=1;i<=nrows;i++) nPresJ += tempMat[i][j];
				if (nPresJ==0) break;
			}
		} while (nPresJ==0);
		tTemp = matrixTemperature(success,0,tempMat,colOrder,rowOrder,nrows,ncols,idum);
		if (success)
		{
			nullTemp[nr] = tTemp;
			if (tTemp < tInt) p2 += 1.0;
		}
		else nr--;
	}
	p2 /= (1.0*NBRS);
	avevar(nullTemp,NBRS,ave,var);
  //save results null model2
	*p_null2 = p2;
	*avt_null2 = ave;
	*avv_null2 = var;
if (bmo==1)  { //save binmatnest in output file
	//out = fopen(OUTFILE.c_str(),"a");
	Rprintf("      Second %11.5f %11.5f %11.5f\n",p2,ave,var); //CFD
	Rprintf("Second null model: p2 = %11.5f\n",p2);
	//fclose(out);
}
	//	Third null model:
	p3 = 0.0;
	for (nr=1;nr<=NBRS;nr++)
	{
		do {
			for (i=1;i<=nrows;i++)
			{
				nPresI = 0;
				while (nPresI==0)
				{
					for (j=1;j<=ncols;j++) 
					{
						p = 0.5*(pRow[i]+pCol[j]);
						r = ran1(idum);
						tempMat[i][j] = (r<p ? 1 : 0);
						nPresI += tempMat[i][j];
					}
				}
			}
			for (j=1;j<=ncols;j++)
			{
				nPresJ = 0;
				for (i=1;i<=nrows;i++) nPresJ += tempMat[i][j];
				if (nPresJ==0) break;
			}
		} while (nPresJ==0);
		tTemp = matrixTemperature(success,0,tempMat,colOrder,rowOrder,nrows,ncols,idum);
		if (success)
		{
			nullTemp[nr] = tTemp;
			if (tTemp < tInt) p3 += 1.0;
		}
		else nr--;
	}
	p3 /= (1.0*NBRS);
	avevar(nullTemp,NBRS,ave,var);
	
	*p_null3 = p3;
	*avt_null3 = ave;
	*avv_null3 = var;

if (bmo==1)  { //save binmatnest in output file
	//out = fopen(OUTFILE.c_str(),"a");
	Rprintf("       Third %11.5f %11.5f %11.5f\n",p3,ave,var);
	Rprintf("Third null model:  p3 = %11.5f\n",p3);
	//fclose(out);
}
	free_ivector(rowOrder,1,nrows);
	free_ivector(colOrder,1,ncols);
	free_imatrix(tempMat,1,nrows,1,ncols);
	free_vector(pRow,1,nrows);
	free_vector(pCol,1,ncols);
	free_vector(nullTemp,1,NBRS);
	 }
 // end of nullmodels.....
/*
	fprintf(stdout, "\n\n\n");
	fprintf(stdout, "BINMATNEST has saved the results.\n");
	fprintf(stdout, "Type a number and press ENTER to continue\n");
        cin >> TOURSIZE;
*/
	return EXIT_SUCCESS;
}
/***
 *****************************************
 ***/
void matrixSize(string inputFile, int& nrows, int& ncols, int& skip)
{
	int j,k;
	char b;	/*caracter leido*/
	FILE *f=0;

	/*abro*/
	if( fopen(inputFile.c_str(),"r")==NULL)
		Rf_error((char*)"Error trying to open input file\n\n");

	endnote(f);   /*me trago la nota */

	nrows=0;
	ncols=0;

	/*	Me situo al principio de la matriz, saltando posibles blancos */
	b=fgetc(f);
	while (b!=CE && b!=UN)
	{
		b=fgetc(f);
		if (b==EOF) Rf_error((char*)"no data found in input matrix\n\n");
	} 

	/*leo la primera fila y asi ya se cuantas columnas hay*/
	while (b==CE || b==UN)  /*hasta fin de linea*/
	{
		ncols++;
		b=fgetc(f);
	}
	skip = 1;
	nrows++;  

	//	Second row: first, we must see how many character separate consecutive rows
	while (true)
	{
		b=fgetc(f);
		if ( (b==CE) || (b==UN) || (b==EOF) ) break;
		skip++;
	}

	if (b == EOF) 
	{
		fclose(f);
		return;
	}
	for (j=2;j<=ncols;j++)
	{
		b=fgetc(f);
		if ((b != CE) && (b !=UN)) Rf_error((char*)"all rows must have the same number of columns\n\n");
	}
	nrows++;  

	/*ahora leo fila a fila hasta el final*/
	while (true)		/*bucle infinito*/
	{
		for (k=1;k<=skip;k++)
		{
			b=fgetc(f);
			if (b == EOF) 
			{
				fclose(f);
				return;
			}
		}

		b=fgetc(f);				/*leo el primer char de la segunda fila*/

		if ( (b != CE) && (b !=UN) )
		{
			fclose(f);
			return;
		}

		/*ahora leo el resto de la fila*/
		for (j=2;j<=ncols;j++)
		{
			b=fgetc(f);
			if ((b != CE) && (b !=UN)) Rf_error((char*)"all rows must have the same number of columns\n\n");
		}
		nrows++;
		/*aumento el contador de filas e intenta leer otra*/
	}
}
/***
 *****************************************
 ***/
void readMatrix(string inputFile, int nrows, int ncols, int skip, int **m)
{
	int i,j,k;
	char b;	/*caracter leido*/
	FILE *f=0;

	/*abro*/
	if( fopen(inputFile.c_str(),"r")==NULL)
		Rf_error((char*)"Error trying to open input file\n\n");

	endnote(f);   /*me trago la nota y ya estoy situado en el inicio de la matriz*/

	/*leo la primera fila */
	/*	Me situo al principio de la matriz, saltando posibles blancos */
	b=fgetc(f);
	while (b!=CE && b!=UN)
	{
		b=fgetc(f);
		if (b==EOF) Rf_error((char*)"no data found in input matrix\n\n");
	} 
	m[1][1] = (b==CE ? 0 : 1);
	for (j=2;j<=ncols;j++)
	{
		b=fgetc(f);
		m[1][j] = (b==CE ? 0 : 1);
	}

	for (i=2;i<=nrows;i++) 		
	{
		for (k=1;k<=skip;k++)
		{
			b=fgetc(f);
			if (b == EOF) Rf_error((char*)"Error reading data");
		}

		/*ahora leo la fila*/
		for (j=1;j<=ncols;j++)
		{
			b=fgetc(f);
			if ((b != CE) && (b !=UN)) Rf_error((char*)"all rows must have the same number of columns\n\n");
			m[i][j] = (b==CE ? 0 : 1);
		}
	}
}
/***
 *****************************************
 ***/
void endnote(FILE *f)
{
	int i;
	char c[8];   /*aqui meto lo que leo*/
	char fin[8]={101,110,100,110,111,116,101,115}; /* codificacion de "endnotes" */
	char fim[8]={69,78,68,78,79,84,69,83}; /* codificacion de "ENDNOTES" */

	/*lee los primeros 8 caracteres*/
	for(i=0;i<8;i++)
	{
		c[i]=fgetc(f);
	}

	while (c[7]!=EOF)
	{
		if ( (c[0]==fin[0] || c[0]==fim[0]) && (c[1]==fin[1] || c[1]==fim[1]) && (c[2]==fin[2] || c[2]==fim[2]) 
			&& (c[3]==fin[3] || c[3]==fim[3]) && (c[4]==fin[4] || c[4]==fim[4]) && (c[5]==fin[5] || c[5]==fim[5]) 
			&& (c[6]==fin[6] || c[6]==fim[6])  && (c[7]==fin[7] || c[7]==fim[7]) )
		{
			c[0]=fgetc(f); /*se supone que despues de endnotes hay un final de linea. lo leo y me lo como*/
			return;
		}
		else
		{
			for (i=0;i<7;i++)	/*mueve a la izquierda*/
				c[i]=c[i+1];
			c[7]=fgetc(f);		/* y lee otro caracter al final */
		}
	}
	Rf_error((char*)"The word endnotes must appear in the input file before the matrix\n\n");
}
/***
 *****************************************
 ***/
double matrixTemperature(bool &success,int sp, int **dataMat, int *c_ord, int *r_ord, int nr, int nc, long& idum)
{
	int i,j;
	int effPlant,effInsect;
	int nbRows,nbCols;
	static int count = 0;
	int **mat;
	int **packed;
	double phi,z;
	double tMat=0;
	double *border;
	double **distance;

	success = true;
	mat = imatrix(1,nr,1,nc);
	for (i=1;i<=nr;i++)
	{
		for (j=1;j<=nc;j++) 
		{
			mat[i][j] = dataMat[i][j];
		}
	}

	//	Permuting rows and columns for pre-packing
	orderMatrix(mat,c_ord,r_ord,nr,nc,effInsect,effPlant);
	removeBlacks(mat,c_ord,r_ord,effInsect,effPlant,nbRows,nbCols,phi);
	packed = imatrix(1,nbRows,1,nbCols);

//	out = fopen(OUTFILE.c_str(),"a");      // ******************
//	fprintf(out,"prepacked matrix bmn5\n\n"); // ******************

	for (i=1;i<=nbRows;i++)
	{
		for (j=1;j<=nbCols;j++) 
		{
			packed[i][j] = mat[i][j];
      //print out pack matrix
//     	fprintf(out,"%1i",packed[i][j]); // ******************

      
		}
//		fprintf(out,"\n");
	}


//fclose(out);  // ******************


	if (nbCols>2 && nbRows>2)
	{
        border = ::vector(1,nbCols);
		calcZ(phi,z);
//    out = fopen(OUTFILE.c_str(),"a");      // ******************
//	fprintf(out,"phi bmn5 %10.5f \n",phi); // ******************
//	fprintf(out,"z bmn5 %10.5f \n",z); // ******************


//  fprintf(out,"border bmn5 border \n"); // ******************
//  for (j=1;j<=nbCols;j++) fprintf(out,"%10.5f",border[j]); // ******************
// 	fprintf(out,"\n"); //************************
//fclose(out);  // ******************
		//	Calculate distance matrix
		distance = matrix(1,nbRows,1,nbCols);
		calcDistance(z,border,distance,nbRows,nbCols);
//****************************
//	out = fopen(OUTFILE.c_str(),"a");      // ******************
//	fprintf(out,"distance matrix bmn5\n\n"); // ******************
//	for (i=1;i<=nbRows;i++)
//	{
//		for (j=1;j<=nbCols;j++)
//		{
		//if (j==nbCols) distance[i][j] *= -1;
//     	fprintf(out,"%10.5f",distance[i][j]); // ******************
//		}
//		fprintf(out,"\n");
//	}
	
//fclose(out);  // ******************
//****************************

		//	Final packing using GA. Returns the temperature of the matrix
		tMat = packMatrix(sp,dataMat,packed,distance,c_ord,r_ord,nbRows,nbCols,nr,nc,idum);

		free_matrix(distance,1,nbRows,1,nbCols);
		free_vector(border,1,nbCols);
	}
	else 
	{
		if (sp) Rf_error((char*)"input matrix must have more than two rows and columns after removing blancks");
		else
		{
			tMat = 0;
			success = false;
			count++;
			if (count>1000) Rf_error((char*)"random matrix has less than two rows or columns too often");
		}
	}
	
	free_imatrix(mat,1,nr,1,nc);
	free_imatrix(packed,1,nbRows,1,nbCols);

	return tMat;
}
/***
 *****************************************
 ***/
void orderMatrix(int **mat, int *c_ord, int *r_ord, int nrows, int ncols, int& effInsect, int& effPlant)
{
	//	We permute rows and columns in such a way that rows with more interactions are closer to the top, and
	//	columns with more interactions are closer to the left.
	int i,j;
	int *lP,*lI;
	int *indP,*indI;
	int **temp;

	lP = ivector(1,ncols);
	indP = ivector(1,ncols);
	lI = ivector(1,nrows);
	indI = ivector(1,nrows);
	temp = imatrix(1,nrows,1,ncols);

	effInsect = effPlant = 0;

	//	We count the number of insects that interact with at least one plant
	for (i=1;i<=nrows;i++)
	{
		indI[i] = i;
		lI[i] = 0;
		for (j=1;j<=ncols;j++) lI[i] -= mat[i][j];
		if (lI[i]<0) effInsect++;
	}
	indexx(nrows,lI,indI);
	for (i=1;i<=nrows;i++) r_ord[i] = indI[i];

	//	We count the number of plants that interact with at least one insect
	for (j=1;j<=ncols;j++)
	{
		indP[j] = j;
		lP[j] = 0;
		for (i=1;i<=nrows;i++) lP[j] -= mat[i][j];
		if (lP[j]<0) effPlant++;
	}
	indexx(ncols,lP,indP);
	for (j=1;j<=ncols;j++) c_ord[j] = indP[j];

	for (i=1;i<=nrows;i++)
	{
		for (j=1;j<=ncols;j++) temp[i][j] = mat[i][j];
	}
	for (i=1;i<=nrows;i++)
	{
		for (j=1;j<=ncols;j++) mat[i][j]=temp[indI[i]][indP[j]];
	}

	free_ivector(lP,1,ncols);
	free_ivector(indP,1,ncols);
	free_ivector(lI,1,nrows);
	free_ivector(indI,1,nrows);
	free_imatrix(temp,1,nrows,1,ncols);
}
/***
 *****************************************
 ***/
void removeBlacks(int **mat, int *c_ord, int *r_ord, int effInsect, int effPlant, int& nbRows, 
				  int& nbCols, double& phi)
{
	int i,j,k;
	int blackInsect=0,blackPlant=0;
	int nbinter=0;
	int **temp;

	temp = imatrix(1,effInsect,1,effPlant);

	//	We count the number of black rows
	for (i=1;i<=effInsect;i++)
	{
		k = 0;
		for (j=1;j<=effPlant;j++) k += mat[i][j];
		if (k==effPlant) blackInsect++;
		else break;
	}
	nbRows = (blackInsect>1 ? effInsect - blackInsect + 1 : effInsect);

	//	We count the number of black columns
	for (j=1;j<=effPlant;j++)
	{
		k = 0;
		for (i=1;i<=effInsect;i++) k += mat[i][j];
		if (k==effInsect) blackPlant++;
		else break;
	}
	nbCols = (blackPlant>1 ? effPlant - blackPlant + 1 : effPlant);

	if (blackInsect>1)
	{
		for (i=1;i<=effInsect;i++)
		{
			for (j=1;j<=effPlant;j++) temp[i][j] = mat[i][j];
		}
		for (i=1;i<=effInsect-blackInsect+1;i++)
		{
			r_ord[i] = r_ord[i+blackInsect-1];
			for (j=1;j<=effPlant;j++) mat[i][j] = temp[i+blackInsect-1][j];
		}
		for (i=effInsect-blackInsect+2;i<=effInsect;i++)
		{
			r_ord[i] = 0;
			for (j=1;j<=effPlant;j++) mat[i][j] = 0;
		}
	}

	if (blackPlant>1)
	{
		for (j=1;j<=effPlant;j++)
		{
			for (i=1;i<=effInsect;i++) temp[i][j] = mat[i][j];
		}
		for (j=1;j<=effPlant-blackPlant+1;j++)
		{
			c_ord[j] = c_ord[j+blackPlant-1];
			for (i=1;i<=effInsect;i++) mat[i][j] = temp[i][j+blackPlant-1];
		}
		for (j=effPlant-blackPlant+2;j<=effPlant;j++)
		{
			c_ord[j] = 0;
			for (i=1;i<=effInsect;i++) mat[i][j] = 0;
		}
	}

	//	We count the number of interactions
	for (i=1;i<=nbRows;i++)
	{
		for (j=1;j<=nbCols;j++) nbinter += mat[i][j];
	}
	if (nbRows>1 && nbCols>1) phi = (2.0*nbinter-nbRows-nbCols+1)/(2.0*(nbRows-1)*(nbCols-1));
	else phi = -1;

	free_imatrix(temp,1,effInsect,1,effPlant);
}
/***
 *****************************************
 ***/
void calcZ(double phi, double& z)
	//	Calculate the parameter for the "line of perfect order"
{
	int i;
	static double znVal[41] = {0.2, 0.224, 0.2509, 0.281, 0.3147, 0.3525, 0.3948, 0.4421, 0.4952, 
		0.5546, 0.6212, 0.6957, 0.7792, 0.8727, 0.9774, 1.0947, 1.2261, 1.3732, 1.538, 1.7226, 
		1.9293, 2.1608, 2.4201, 2.7105, 3.0357, 3.4, 3.808, 4.265, 4.7768, 5.35, 5.992, 6.711, 
		7.5163, 8.4183, 9.4285, 10.5599, 11.8271, 13.2464, 14.8359, 16.6162, 18.6102};
	static double propOc[41] = {0.996, 0.9921, 0.9855, 0.9751, 0.9599, 0.9389, 0.9116, 0.8776, 0.8371, 
		0.7907, 0.7394, 0.6844, 0.6272, 0.5691, 0.5114, 0.4555, 0.402, 0.3519, 0.3057, 0.2636, 0.2258, 
		0.1922, 0.1626, 0.1368, 0.1146, 0.0956, 0.0793, 0.0656, 0.0541, 0.0445, 0.0365, 0.0297, 
		0.0242, 0.0197, 0.0161, 0.013, 0.0106, 0.0086, 0.007, 0.0057, 0.0046};

	if (phi>=1) z = 1000;
	else if (phi>0)
	{
		if (phi>=propOc[0])
		{
			z = znVal[0]*(1-phi)/(1-propOc[0]);
		}
		else if (phi<=propOc[40])
		{
			z = znVal[40];
		}
		else
		{
			for (i=1;i<=40;i++)
			{
				if (propOc[i]<=phi) break;
			}
			z = znVal[i-1]+(znVal[i]-znVal[i-1])*(propOc[i-1]-phi)/(propOc[i-1]-propOc[i]);
		}
	}
	else z = -1;
}
/***
 *****************************************
 ***/
void calcDistance(double z,double border[],double **mat,int nr,int nc)
{
	int i,j;
	double d2,dmax2;
	double s,x,x0,x1,xx,y,y0,y1,yy;

	if (z>0 && z<100)
	{

//    out = fopen(OUTFILE.c_str(),"a");
  	for (j=1;j<=nc;j++)
		{
			x = (j-0.5)/nc;
			xx = (x-0.5/nc)*nc/(nc-1);
      if (xx==1) s=0; else    s = pow(1-xx,z);
      if (s==1) yy=0; else 		yy = pow(1-s,1/z);
			y = (0.5+(nr-1)*yy)/nr;
			border[j] = y*nr;
			
//   fprintf(out,"%10.5f",x);
//   fprintf(out,"%10.5f",xx);
//   fprintf(out,"%10.5f",s);
//   fprintf(out,"%10.5f",yy);
//   fprintf(out,"%10.5f",y);
//   fprintf(out,"%10.5f",border[j]);
//   fprintf(out,"  \n"); // ******************

		}
//    fclose(out);  // ******************

//  out = fopen(OUTFILE.c_str(),"a");      // ******************

//  fprintf(out,"border bmn5 border \n"); // ******************
//  for (j=1;j<=nc;j++) fprintf(out,"%10.5f",border[j]); // ******************
// 	fprintf(out,"\n"); //************************
//  fclose(out);  // ******************





//out = fopen(OUTFILE.c_str(),"a");      // ******************

//  fprintf(out,"i,j bmn5  \n"); // ******************




		for (i=1;i<=nr;i++)
		{
			for (j=1;j<=nc;j++)
			{
				if ((i==1 && j==nc) || (i==nr && j==1)) mat[i][j] = 0;
				else
				{
					x1 = (j-0.5)/nc;
					y1 = (nr-i+0.5)/nr;
//					dmax2 = (x1+y1<1 ? 2*(x1+y1)*(x1+y1) : 2*(2-x1-y1)*(2-x1-y1));
//					but we calculate only half of the distance:
					dmax2 = (x1+y1<1 ? (x1+y1)*(x1+y1) : (2-x1-y1)*(2-x1-y1));
					yy = zbrent(nr,nc,x1,y1,z);
					

					y = (0.5+(nr-1)*yy)/nr;
/*	
					d2 = (x1-x)*(x1-x) + (y1-y)*(y1-y);
					with:
					x = x1 + y1 - y;
					so that x1-x = y-y1
					and d2 = (y-y1)*(y-y1) + (y1-y)*(y1-y) = 2*(y1-y)*(y1-y);
					We calculate only half of the distance:
*/	
					d2 = (y1-y)*(y1-y);
					mat[i][j] = d2/dmax2;
					if ((nr-i+0.5)<border[j])  mat[i][j] *= -1;
 //                  fprintf(out,"%1i/%1i, ",i,j);} // ******************
				}
			}
		}
//		fclose(out);  // ******************

	}
	else if (z>100)
	{
		for (j=1;j<=nc;j++) border[j] = 0.5;

		for (i=1;i<=nr;i++)
		{
			for (j=1;j<=nc;j++)
			{
				if ((i==1 && j==nc) || (i==nr && j==1)) mat[i][j] = 0;
				else
				{
					x1 = j-0.5;
					y1 = nr-i+0.5;
					if (y1<=nr*(1-x1/nc))
					{
						x0 = x1 + (nc*y1)/nr;
						y0 = y1 + (nr*x1)/nc;
						dmax2 = x0*x0 + y0*y0;
						d2 = (1+(1.0*nc*nc)/(nr*nr))*(nr-i)*(nr-i);
						mat[i][j] = d2/dmax2;
						if ((nr-i+0.5)<border[j]) mat[i][j] *= -1;
					}
					else
					{
						x0 = 2*nr - x1 - (nc*y1)/nr;
						y0 = 2*nc - y1 - (nr*x1)/nc;
						dmax2 = x0*x0 + y0*y0;
						d2 = (1+(1.0*nr*nr)/(nc*nc))*(nc-j)*(nc-j);
						mat[i][j] = d2/dmax2;
						if ((nr-i+0.5)<border[j]) mat[i][j] *= -1;
					}
				}
			}
		}
	}
	else
	{
		border[1] = 0.5;
		for (j=2;j<=nc;j++) border[j] = nr - 0.5;

		for (i=1;i<=nr;i++)
		{
			for (j=1;j<=nc;j++)
			{
				if ((i==1 && j==nc) || (i==nr && j==1)) mat[i][j] = 0;
				else
				{
					x1 = j-0.5;
					y1 = nr-i+0.5;
					if (y1<=nr*(1-x1/nc))
					{
						x0 = x1 + (nc*y1)/nr;
						y0 = y1 + (nr*x1)/nc;
						dmax2 = x0*x0 + y0*y0;
						d2 = (1+(1.0*nr*nr)/(nc*nc))*(j-1)*(j-1);
						mat[i][j] = d2/dmax2;
						if ((nr-i+0.5)<border[j]) mat[i][j] *= -1;
					}
					else
					{
						x0 = 2*nr - x1 - (nc*y1)/nr;
						y0 = 2*nc - y1 - (nr*x1)/nc;
						dmax2 = x0*x0 + y0*y0;
						d2 = (1+(1.0*nc*nc)/(nr*nr))*(i-1)*(i-1);
						mat[i][j] = d2/dmax2;
						if ((nr-i+0.5)<border[j]) mat[i][j] *= -1;
					}
				}
			}
		}
	}
}
/***
 *****************************************
 ***/
double packMatrix(int sp, int **dataMat, int **intMat, double **d, int *c_ord, int *r_ord, int nr, int nc, 
				  int nrows, int ncols, long& idum)
	//	indr[i] indica la fila de la matriz que deben ocupar los datos del insecto i
	//	indc[j] indica la columna de la matriz que deben ocupar los datos de la planta j
{
	int i,j,l,k1,k2,ng;
	int *indr,*indc;
	int *pr,*or_new,*pc,*oc;
	int *players,*sorted;
	int *bestr,*bestc;
	int *unused;
	int **pRows,**pCols;
	int **mat;
	double tMin;
	double *temp,*tPlayers;
	FILE *out=NULL;

	indr = ivector(1,nr);
	indc = ivector(1,nc);
	pr = ivector(1,nr);
	pc = ivector(1,nc);
	or_new = ivector(1,nr);
	oc = ivector(1,nc);
	bestr = ivector(1,nr);
	bestc = ivector(1,nc);
	players = ivector(1,TOURSIZE);
	sorted = ivector(1,TOURSIZE);
    tPlayers = ::vector(1,TOURSIZE);
    temp = ::vector(1,POPSIZE);
	pRows = imatrix(1,POPSIZE,1,nr);
	pCols = imatrix(1,POPSIZE,1,nc);

	for (i=1;i<=nr;i++) pRows[1][i] = indr[i] = bestr[i] = i;
	for (j=1;j<=nc;j++) pCols[1][j] = indc[j] = bestc[j] = j;
	tMin = temp[1] = calcTemp(d,intMat,indr,indc,nr,nc);
	for (l=2;l<=12;l++)
	{
		for (i=1;i<=nr;i++) pRows[l][i] = i;
		for (j=1;j<=nc;j++) pCols[l][j] = j;
		prePackMatrix(intMat,pRows[l],pCols[l],nr,nc,0.1*(l-2));
	}
	//	this one is created with the algorithm from NTC
	for (i=1;i<=nr;i++) pRows[13][i] = i;
	for (j=1;j<=nc;j++) pCols[13][j] = j;
	prePackNTC(intMat,pRows[13],pCols[13],nr,nc);
	for (l=14;l<=POPSIZE;l++)
	{
		k1 = 1 + l%13;
		for (i=1;i<=nr;i++) pRows[l][i] = pRows[k1][i];
		for (j=1;j<=nc;j++) pCols[l][j] = pCols[k1][j];
		if (ran1(idum)<0.7) mutate(idum,nr,pRows[l]);
		if (ran1(idum)<0.7) mutate(idum,nc,pCols[l]);
	}

	for (l=2;l<=POPSIZE;l++)
	{
		temp[l] = calcTemp(d,intMat,pRows[l],pCols[l],nr,nc);
		if (temp[l]<tMin)
		{
			tMin = temp[l];
			for (i=1;i<=nr;i++) bestr[i] = pRows[l][i];
			for (j=1;j<=nc;j++) bestc[j] = pCols[l][j];
		}
	}

	for (ng=1;ng<=NBGENER;ng++)
	{
		choosePlayers(idum,TOURSIZE,POPSIZE,players);
		for (l=1;l<=TOURSIZE;l++) tPlayers[l] = temp[players[l]];
		indexxD(TOURSIZE,tPlayers,sorted);

		k1 = players[sorted[1]];
		for (i=1;i<=nr;i++) or_new[i] = pRows[k1][i];
		for (j=1;j<=nc;j++) oc[j] = pCols[k1][j];
		do
		{
			k2 = 1 + int(POPSIZE*ran1(idum));
		} while (k1==k2 || k2>POPSIZE);
		for (i=1;i<=nr;i++) pr[i] = pRows[k2][i];
		for (j=1;j<=nc;j++) pc[j] = pCols[k2][j];
		crossOver(idum,nr,nc,pr,pc,or_new,oc);
		l = players[sorted[TOURSIZE]];
		for (i=1;i<=nr;i++) pRows[l][i] = or_new[i];
		for (j=1;j<=nc;j++) pCols[l][j] = oc[j];
		temp[l] = calcTemp(d,intMat,or_new,oc,nr,nc);
		if (temp[l]<tMin)
		{
			tMin = temp[l];
			for (i=1;i<=nr;i++) bestr[i] = or_new[i];
			for (j=1;j<=nc;j++) bestc[j] = oc[j];
		}

		k1 = players[sorted[2]];
		for (i=1;i<=nr;i++) or_new[i] = pRows[k1][i];
		for (j=1;j<=nc;j++) oc[j] = pCols[k1][j];
		do
		{
			k2 = 1 + int(POPSIZE*ran1(idum));
		} while (k1==k2 || k2>POPSIZE);
		for (i=1;i<=nr;i++) pr[i] = pRows[k2][i];
		for (j=1;j<=nc;j++) pc[j] = pCols[k2][j];
		crossOver(idum,nr,nc,pr,pc,or_new,oc);
		l = players[sorted[TOURSIZE-1]];
		for (i=1;i<=nr;i++) pRows[l][i] = or_new[i];
		for (j=1;j<=nc;j++) pCols[l][j] = oc[j];
		temp[l] = calcTemp(d,intMat,or_new,oc,nr,nc);
		if (temp[l]<tMin)
		{
			tMin = temp[l];
			for (i=1;i<=nr;i++) bestr[i] = or_new[i];
			for (j=1;j<=nc;j++) bestc[j] = oc[j];
		}
	}
	
//if (bmo==1) 	out = fopen(OUTFILE.c_str(),"a");
	
if (bmo==1) {
	if (sp)
	{
	//	out = fopen(OUTFILE.c_str(),"a");

		fprintf(out,"Packed matrix:\n"); //CFD

		for (i=1;i<=nr;i++)
		{
			for (j=1;j<=nc;j++)
			{
				fprintf(out,"%1i",intMat[bestr[i]][bestc[j]]);
			}
			fprintf(out,"\n");
		}
		fprintf(out,"\n");   //CFD

		fclose(out);


		calcIdiosyncTemp(d,intMat,or_new,oc,nr,nc);
	}
}
	free_ivector(indr,1,nr);
	free_ivector(indc,1,nc);
	free_ivector(pr,1,nr);
	free_ivector(pc,1,nc);
	free_ivector(or_new,1,nr);
	free_ivector(oc,1,nc);
	free_ivector(players,1,TOURSIZE);
	free_ivector(sorted,1,TOURSIZE);
	free_vector(tPlayers,1,TOURSIZE);
	free_vector(temp,1,POPSIZE);
	free_imatrix(pRows,1,POPSIZE,1,nr);
	free_imatrix(pCols,1,POPSIZE,1,nc);


	if (sp)
	{

//    fopen(OUTFILE.c_str(),"a");
		sorted = ivector(1,nr);
		for (i=1;i<=nr;i++) sorted[i] = r_ord[i];
		for (i=1;i<=nr;i++) r_ord[i] = sorted[bestr[i]];
		free_ivector(sorted,1,nr);

		unused = ivector(1,nrows);
		for (i=1;i<=nrows;i++) unused[i] = 1;
		for (i=1;i<=nr;i++) unused[r_ord[i]] = 0;
		for (i=1;i<=nrows;i++) 
		{
			if (unused[i])
			{
				k1 = 0;
				for (j=1;j<=ncols;j++) k1 += dataMat[i][j];
				unused[i] = (k1 > 0? 1 : nr+1);
			}
		}

if (bmo==1)  fprintf(out,"Permutation of rows:\n");
if (bmo==1)  fprintf(out,"Row number            ");
		for (i=1;i<=nrows;i++)
		{
			if (unused[i]==1)  if (bmo==1) fprintf(out,"%5i",unused[i]);
		}
		for (i=1;i<=nr;i++)  if (bmo==1) fprintf(out,"%5i",i);
		for (i=1;i<=nrows;i++)
		{
			if (unused[i]>1)  if (bmo==1) fprintf(out,"%5i",unused[i]);
		}
		if (bmo==1) fprintf(out,"\n");
 		if (bmo==1) fprintf(out,"comes from position   ");
		for (i=1;i<=nrows;i++)
		{
			if (unused[i]==1) { if (bmo==1) fprintf(out,"%5i",i); tor[i-1]=i;}
		}
		for (i=1;i<=nr;i++) { if (bmo==1) fprintf(out,"%5i",r_ord[i]); tor[i-1] = r_ord[i];}
		for (i=1;i<=nrows;i++)
		{
			if (unused[i]>1) { if (bmo==1) fprintf(out,"%5i",i); tor[i-1] =i;}
		}
		if (bmo==1) fprintf(out,"\n\n");
		free_ivector(unused,1,nrows);

		sorted = ivector(1,nc);
		for (j=1;j<=nc;j++) sorted[j] = c_ord[j];
		for (j=1;j<=nc;j++) c_ord[j] = sorted[bestc[j]];
		free_ivector(sorted,1,nc);

		unused = ivector(1,ncols);
		for (j=1;j<=ncols;j++) unused[j] = 1;
		for (j=1;j<=nc;j++) unused[c_ord[j]] = 0;
		for (j=1;j<=ncols;j++) 
		{
			if (unused[j])
			{
				k1 = 0;
				for (i=1;i<=nrows;i++) k1 += dataMat[i][j];
				unused[j] = (k1 > 0? 1 : nc+1);
			}
		}

    if (bmo==1) fprintf(out,"Permutation of columns:\n");
 		if (bmo==1) fprintf(out,"Column number         ");
		for (j=1;j<=ncols;j++)
		{
			if (unused[j]==1)  if (bmo==1) fprintf(out,"%5i",unused[j]);
		}
		for (j=1;j<=nc;j++)  if (bmo==1) fprintf(out,"%5i",j);
		for (j=1;j<=ncols;j++)
		{
			if (unused[j]>1)  if (bmo==1) fprintf(out,"%5i",unused[j]);
		}
		if (bmo==1) fprintf(out,"\n");
 		if (bmo==1) fprintf(out,"comes from position   ");
		for (j=1;j<=ncols;j++)
		{
			if (unused[j]==1) { if (bmo==1) fprintf(out,"%5i",j); toc[j-1] = j;}
		}
		for (j=1;j<=nc;j++) { if (bmo==1) fprintf(out,"%5i",c_ord[j]);toc[j-1] = c_ord[j];}
		for (j=1;j<=ncols;j++)
		{
			if (unused[j]>1) { if (bmo==1) fprintf(out,"%5i",j); toc[j-1] = j;}
		}
 		if (bmo==1) fprintf(out,"\n\n");
		free_ivector(unused,1,ncols);

 		//fclose(out);

  }
  
if (bmo==1)  	fclose(out);
	mat = imatrix(1,nr,1,nc);
	for (i=1;i<=nr;i++)
	{
		for (j=1;j<=nc;j++) mat[i][j] = intMat[i][j];
	}
	for (i=1;i<=nr;i++)
	{
		for (j=1;j<=nc;j++) intMat[i][j] = mat[bestr[i]][bestc[j]];
	}
	free_imatrix(mat,1,nr,1,nc);
	free_ivector(bestr,1,nr);
	free_ivector(bestc,1,nc);
	return tMin;
}
/***
 *****************************************
 ***/
void prePackMatrix(int **mat, int indr[], int indc[], int nr, int nc, double x)
{
	int fold;

	if (nc > nr)
	{
		for (fold=1;fold<=4;fold++)
		{
			prePackcols(mat,indr,indc,nr,nc,x);
			prePackrows(mat,indr,indc,nr,nc,x);
		}
	}
	else
	{
		for (fold=1;fold<=4;fold++)
		{
			prePackrows(mat,indr,indc,nr,nc,x);
			prePackcols(mat,indr,indc,nr,nc,x);
		}
	}
}
/***
 *****************************************
 ***/
void prePackcols(int **mat, int indr[], int indc[], int nr, int nc, double x)
{
	int i,j;
	double *w;

    w = ::vector(1,nc);

	//	We calculate an index for each column (plant)
	for (j=1;j<=nc;j++)
	{
		w[j] = 0;
		for (i=1;i<=nr;i++) 
		{
			if (mat[indr[i]][j]==1) w[j] -= x*i*i;
			else w[j] += (1-x)*(nr-i+1)*(nr-i+1);
		}
	}

	indexxD(nc,w,indc);

	free_vector(w,1,nc);
}
/***
 *****************************************
 ***/
void prePackrows(int **mat, int indr[], int indc[], int nr, int nc, double x)
{
	int i,j;
	double *w;

    w = ::vector(1,nr);

	//	We calculate an index for each column (plant)
	for (i=1;i<=nr;i++)
	{
		w[i] = 0;
		for (j=1;j<=nc;j++)
		{
			if (mat[i][indc[j]]==1) w[i] -= x*j*j;
			else w[i] += (1-x)*(nc-j+1)*(nc-j+1);
		}
	}

	indexxD(nr,w,indr);

	free_vector(w,1,nr);
}
/***
 *****************************************
 ***/
void prePackNTC(int **mat, int indr[], int indc[], int nr, int nc)
{
	int fold;

	if (nc > nr)
	{
		for (fold=1;fold<=8;fold++)
		{
			prePackNTCcols(mat,indr,indc,nr,nc);
			prePackNTCrows(mat,indr,indc,nr,nc);
		}
	}
	else
	{
		for (fold=1;fold<=8;fold++)
		{
			prePackNTCrows(mat,indr,indc,nr,nc);
			prePackNTCcols(mat,indr,indc,nr,nc);
		}
	}
}
/***
 *****************************************
 ***/
void prePackNTCcols(int **mat, int indr[], int indc[], int nr, int nc)
	//	indr[i] indica la fila de la matriz que deben ocupar los datos del insecto i
	//	indc[j] indica la columna de la matriz que deben ocupar los datos de la planta j
{
	int i,j;
	int *ndum,*ns,*nt;
	double *s,*t;

	ndum = ivector(1,nc);
	ns = ivector(1,nc);
	nt = ivector(1,nc);
    s = ::vector(1,nc);
    t = ::vector(1,nc);

	//	We calculate an index for each column (plant)
	for (j=1;j<=nc;j++)
	{
		s[j] = t[j] = 0;
		for (i=1;i<=nr;i++) 
		{
			if (mat[indr[i]][j]==1) s[j] -= i*i;
			else t[j] += (nr-i+1)*(nr-i+1);
		}
	}
	indexxD(nc,s,ndum);
	indexx(nc,ndum,ns);
	indexxD(nc,t,ndum);
	indexx(nc,ndum,nt);
	for (j=1;j<=nc;j++) s[j] = ns[j] + nt[j];
	indexxD(nc,s,indc);

	free_ivector(ndum,1,nc);
	free_ivector(ns,1,nc);
	free_ivector(nt,1,nc);
	free_vector(s,1,nc);
	free_vector(t,1,nc);
}
/***
 *****************************************
 ***/
void prePackNTCrows(int **mat, int indr[], int indc[], int nr, int nc)
	//	indr[i] indica la fila de la matriz que deben ocupar los datos del insecto i
	//	indc[j] indica la columna de la matriz que deben ocupar los datos de la planta j
{
	int i,j;
	int *ndum,*ns,*nt;
	double *s,*t;

	ndum = ivector(1,nr);
	ns = ivector(1,nr);
	nt = ivector(1,nr);
    s = ::vector(1,nr);
    t = ::vector(1,nr);

	//	We calculate an index for each column (plant)
	for (i=1;i<=nr;i++)
	{
		s[i] = t[i] = 0;
		for (j=1;j<=nc;j++) 
		{
			if (mat[i][indc[j]]==1) s[i] -= j*j;
			else t[i] += (nc-j+1)*(nc-j+1);
		}
	}
	indexxD(nr,s,ndum);
	indexx(nr,ndum,ns);
	indexxD(nr,t,ndum);
	indexx(nr,ndum,nt);
	for (i=1;i<=nr;i++) s[i] = ns[i] + nt[i];
	indexxD(nr,s,indr);
}
/***
 *****************************************
 ***/
double calcTemp(double **d, int **mat, int indr[], int indc[], int nr, int nc)
{
	//	indr[i] indica la fila de la matriz que deben ocupar los datos del insecto i
	//	indc[j] indica la columna de la matriz que deben ocupar los datos de la planta j
	int i,j;
	int i1,j1;
	double unex;
	double tMat;
	unex = 0.0;
	for (i=1;i<=nr;i++)
	{
		i1 = indr[i];
		for (j=1;j<=nc;j++)
		{
			j1 = indc[j];
			if ((mat[i1][j1]==1 && d[i][j]<0) || (mat[i1][j1]==0 && d[i][j]>0))
			{
				unex += fabs(d[i][j]);
			}
		}
	}
	unex /= (nr*nc);
	tMat = 100*unex/UMAX;

	return tMat;
}
/***
 *****************************************
 ***/
void calcIdiosyncTemp(double **d, int **mat, int indr[], int indc[], int nr, int nc)
{
	//	indr[i] indica la fila de la matriz que deben ocupar los datos del insecto i
	//	indc[j] indica la columna de la matriz que deben ocupar los datos de la planta j
	int i,j;
	int i1,j1;
	double unex;
	double tIdiosync;
	FILE *out=0;
    
if (bmo==1) {
	//out = fopen(OUTFILE.c_str(),"a");
	//	Idiosyncratic temperature for rows:
	Rprintf("Idiosyncratic temperature for rows:\n");
	Rprintf("Row:                       ");
	for (i=1;i<=nr;i++) fprintf(out,"%10i",i);
	Rprintf("\n");
	Rprintf("Idiosyncratic temperature: ");
	for (i=1;i<=nr;i++)
	{
		unex = 0.0;
		i1 = indr[i];
		for (j=1;j<=nc;j++)
		{
			j1 = indc[j];
			if ((mat[i1][j1]==1 && d[i][j]<0) || (mat[i1][j1]==0 && d[i][j]>0))
			{
				unex += fabs(d[i][j]);
			}
		}
		unex /= nc;
		tIdiosync = 100*unex/UMAX;
		Rprintf("%10.5f",tIdiosync);
	}
	Rprintf("\n\n");

	//	Idiosyncratic temperature for columns:
	Rprintf("Idiosyncratic temperature for columns:\n");
	Rprintf("Column:                    ");
	for (j=1;j<=nc;j++) Rprintf("%10i",j);
	Rprintf("\n");
	Rprintf("Idiosyncratic temperature: ");
	for (j=1;j<=nc;j++)
	{
		unex = 0.0;
		j1 = indc[j];
		for (i=1;i<=nr;i++)
		{
			i1 = indr[i];
			if ((mat[i1][j1]==1 && d[i][j]<0) || (mat[i1][j1]==0 && d[i][j]>0))
			{
				unex += fabs(d[i][j]);
			}
		}
		unex /= nr;
		tIdiosync = 100*unex/UMAX;
		Rprintf("%10.5f",tIdiosync);
	}
	Rprintf("\n\n");

	//fclose(out);
}
}
/***
 *****************************************
 ***/
void permute(long& idum, int n, int index[])
//	Takes an array index[i] as imput, with 1=1,..., n and permutes it at random
{
	int i,j,k,m=n;
	int *perm;
	perm = ivector(1,n);
	for (i=1;i<=n;i++)
	{
		j = 1 + int(m*ran1(idum));
		if (j<=m)
		{
			perm[i] = index[j];
			for (k=j;k<m;k++) index[k] = index[k+1];
			m--;
		}
		else i--;
	}
	for (i=1;i<=n;i++) index[i] = perm[i];
	free_ivector(perm,1,n);
}
/***
 *****************************************
 ***/
void mutate(long& idum, int n, int index[])
//	Takes an array index[i] as imput, with 1=1,..., n and returns a random mutation
//	Mutations consist of a cyclic permutation of a part of the sequence
{
	int i,i1,i2,itoto;

	i1 = 1 + int(n*ran1(idum));
	i2 = 1 + int(n*ran1(idum));
	if (i1 < i2)
	{
		itoto = index[i2];
		for (i=i2;i>i1;i--) index[i] = index[i-1];
		index[i1] = itoto;
	}
	else if (i1 > i2)
	{
		itoto = index[i2];
		for (i=i2;i<i1;i++) index[i] = index[i+1];
		index[i1] = itoto;
	}
}
/***
 *****************************************
 ***/
void choosePlayers(long& idum, int n,int m,int arr[])
//	Chooses n numbers from integers 1, 2... m, and returns them as the elements of arr[]
{
	int i,j,k;
	int *index;

	index = ivector(1,m);
	if (n>m) Rf_error((char*)"n too large in choosePlayers");
	else if (n==m)
	{
		for (i=1;i<=n;i++) arr[i] = i;
	}
	else
	{
		for (i=1;i<=m;i++) index[i] = i;
		for (i=1;i<=n;i++)
		{
			j = 1 + int(m*ran1(idum));
			if (j<=m)
			{
				arr[i] = index[j];
				for (k=j;k<m;k++) index[k] = index[k+1];
				m--;
			}
			else i--;
		}
	}
	free_ivector(index,1,n);
}
/***
 *****************************************
 ***/
void crossOver(long& idum, int nr, int nc, int pr[], int pc[], int or_new[], int oc[])
//	At input, the arrays pr[] and or[] contain a permutation of row indexes from two individuals of the population,
//	while pc[] and oc[] contain permutations of the column indexes. At output, pr and pc are not changed, but or and oc
//	contain new permutations, obtained from "crossing" over the corresponding input strategies, and possibly by mutating 
//	the resulting strategy.
{
	int i,j,k,m;
	int flag=0;
	int *ind;

	while (flag==0)
	{
		if (ran1(idum)<0.5)
		{
			//	crossing over of rows
			ind = ivector(1,nr);
			for (i=1;i<=nr;i++) ind[i] = 0;
			k = 2 + int((nr-2)*ran1(idum));
			for (i=1;i<=k;i++) ind[or_new[i]] = 1;
			for (i=k+1;i<=nr;i++)
			{
				if (ind[pr[i]]==0)
				{
					ind[pr[i]] = 1;
					or_new[i] = pr[i];
				}
				else or_new[i] = 0;
			}
			m = 0;
			for (i=1;i<=nr;i++)
			{
				if (ind[i]==0)
				{
					m++;
					ind[m] = i;
				}
			}
			if (m>1) permute(idum,m,ind);
			if (m>0)
			{
				for (i=1;i<=nr;i++)
				{
					if (or_new[i]==0)
					{
						if (m<=0) Rf_error((char*)"problem in crossOver, rows");
						or_new[i] = ind[m];
						m--;
					}
				}
			}
			free_ivector(ind,1,nr);
			flag = 1;
		}
		if (ran1(idum)<0.5)
		{
			//	crossing over of columns
			ind = ivector(1,nc);
			for (j=1;j<=nc;j++) ind[j] = 0;
			k = 2 + int((nc-2)*ran1(idum));
			for (j=1;j<=k;j++) ind[oc[j]] = 1;
			for (j=k+1;j<=nc;j++)
			{
				if (ind[pc[j]]==0)
				{
					ind[pc[j]] = 1;
					oc[j] = pc[j];
				}
				else oc[j] = 0;
			}
			m = 0;
			for (j=1;j<=nc;j++)
			{
				if (ind[j]==0)
				{
					m++;
					ind[m] = j;
				}
			}
			if (m>1) permute(idum,m,ind);
			if (m>0)
			{
				for (j=1;j<=nc;j++)
				{
					if (oc[j]==0)
					{
						if (m<=0) Rf_error((char*)"problem in crossOver, columns");
						oc[j] = ind[m];
						m--;
					}
				}
			}
			free_ivector(ind,1,nc);
			flag = 1;
		}
	}
	//	now mutate:
	if (ran1(idum)<0.1) mutate(idum,nr,or_new);
	if (ran1(idum)<0.1) mutate(idum,nc,oc);
}
/***
 *****************************************
 ***/
void indexx(int n, int arr[], int indx[])
{
	int i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	int a;
   
	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M_IND) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1])
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) Rf_error((char*)"NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
/***
 *****************************************
 ***/
void indexxD(int n, double arr[], int indx[])
{
	int i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	double a;
   
	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M_IND) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1])
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) Rf_error((char*)"NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
/***
 *****************************************
 ***/
double zbrent(int nr, int nc, double u1, double v1, double z)
{
	int iter;
	double x1,x2=1,tol=1.0e-5;
	double a,b=x2,c=x2,d=0,e=0,min1,min2;
	double fa,fb,fc,p,q,r,s,tol1,xm;

	x1 = (u1+v1<1 ? 0 : u1+v1-1);
	a = x1;
	fa=func(a,nr,nc,u1,v1,z);
	fb=func(b,nr,nc,u1,v1,z);
   
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
	{
		Rf_error((char*)"Root must be bracketed in zbrent");
	}
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPSB*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=func(b,nr,nc,u1,v1,z);
	}
	Rf_error((char*)"Maximum number of iterations exceeded in zbrent");
	return 0.0;
}
/***
 *****************************************
 ***/
double func(double yy, int nr, int nc, double x1, double y1, double z)
{
	double f,s,x;
	x = nc*(x1+y1-(nr-1)*yy/nr-0.5/nc-0.5/nr)/(nc-1.0);
	if (fabs(x)<EPS) s = 1;
	else if (x < 1) s = pow(1-x,z);
	else s = 0;
	if (fabs(1-yy)<EPS) f = s;
	else if (yy>0) f = pow(yy,z) + s - 1;
	else f = s - 1;
	return f;
}
/***
 *****************************************
 ***/
double ran1(long& idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
   
	if (idum <= 0 || !iy) {
		if (-idum < 1) idum=1;
		else idum = -idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
/***
 *****************************************
 ***/
void avevar(double data[], unsigned long n, double& ave, double& var)
{
	unsigned long j;
	double s,ep;
   
	for (ave=0.0,j=1;j<=n;j++) ave += data[j];
	ave /= n;
	var=ep=0.0;
	for (j=1;j<=n;j++) {
		s=data[j]-ave;
		ep += s;
		var += s*s;
	}
	var=(var-ep*ep/n)/(n-1);
}
/***
 *****************************************
 ***/
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;
   
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) Rf_error((char*)"allocation failure in ivector()");
	return v-nl+NR_END;
}
/***
 *****************************************
 ***/
double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
   
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) Rf_error((char*)"allocation failure in dvector()");
	return v-nl+NR_END;
}
/***
 *****************************************
 ***/
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
/***
 *****************************************
 ***/
void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
/***
 *****************************************
 ***/
int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
   
	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) Rf_error((char*)"allocation failure 1 in imatrix()");
	m += NR_END;
	m -= nrl;
   
	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) Rf_error((char*)"allocation failure 2 in imatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
   
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
   
	/* return pointer to array of pointers to rows */
	return m;
}
/***
 *****************************************
 ***/
double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
   
	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) Rf_error((char*)"allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
   
	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) Rf_error((char*)"allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
   
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
   
	/* return pointer to array of pointers to rows */
	return m;
}
/***
 *****************************************
 ***/
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an integer matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
/***
 *****************************************
 ***/
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
/***
 *****************************************
 ***/
/* void Rf_error(char error_text[])
    // Numerical Recipes standard error handler 
{
	fprintf(stderr,"\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"Run-time error...\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"\n\n");
	fprintf(stdout, "Type a number and press ENTER to continue\n");
	fprintf(stderr,"\n");
        cin >> TOURSIZE;

	exit(1);
}   */

    /***
 *****************************************
 ***/



   
}

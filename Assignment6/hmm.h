#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// pi[N], A[N][N] and B[N][M]
#define N 2
#define M 27
//#define M 26

// character set
// Note: Size of character set must be M
#define ALPHABET {"abcdefghijklmnopqrstuvwxyz "}
//#define ALPHABET {"abcdefghijklmnopqrstuvwxyz"}

// maximum characters per line
#define MAX_CHARS 500

// other
#define EPSILON 0.00001
#define DABS(x) ((x) < (0.0) ? -(x) : (x))

// debugging and/or printing
//#define PRINT_OBS
//#define PRINT_GET_T
//#define CHECK_GAMMAS
#define PRINT_REESTIMATES

struct stepStruct
{
    int obs;
    double c;
    double alpha[N];
    double beta[N];
    double gamma[N];
    double diGamma[N][N];
};

void alphaPass(struct stepStruct *step,
               double pi[], 
               double A[][N],
               double B[][M],
               int T);

void betaPass(struct stepStruct *step,
              double pi[], 
              double A[][N],
              double B[][M],
              int T);

void computeGammas(struct stepStruct *step,
                   double pi[], 
                   double A[][N],
                   double B[][M],
                   int T);
                   
void reestimatePi(struct stepStruct *step, 
                  double piBar[]);
                  
void reestimateA(struct stepStruct *step, 
                 double Abar[][N], 
                 int T);

void reestimateB(struct stepStruct *step, 
                 double Bbar[][M], 
                 int T);

void initMatrices(double pi[], 
                  double A[][N], 
                  double B[][M],
                  int seed);

int GetT(char fname[],
         int startPos,
         int startChar,
         int maxChars);

int GetObservations(char fname[], 
                    struct stepStruct *step,
                    int T,
                    int startPos,
                    int startChar,
                    int maxChars);
                                  
void printPi(double pi[]);

void printA(double A[][N]);

void printBT(double B[][M]);

//
// Hidden Markov Model program for written English
//
// Consistent with "Revealing Introduction" version dated January 12, 2017
//
// The program reads in data (English text) ignoring
// all punctuation, leaving only 26 letters and
// word spaces. Consequently, the A matrix is N x N and 
// the B matrix is N x M, where M = 27. The program
// begins with random (approximately uniform) A, B and pi,
// then attempts to climb to the optimal A and B.
// This program implements the work described in 
// Cave and Neuwirth's paper "Hidden Markov Models
// for English"
//
// Note: The data file must not have any formatting!
//
// To compile: gcc -o hmm hmm.c -O3
//

#include "Q15.h"

int randSubs = 0;

int main(int argc, const char *argv[])
{
    int startPos,
        startChar,
        maxChars,
        maxIters,
        i,
        j,
        T,
        iter,
        No_Iter;
        
    int seed;

    double logProb,
           newLogProb;

    double pi[N],
           piBar[N],
           A[N][N],
           Abar[N][N],
           B[N][M],
           Bbar[N][M];
    int cnt = 0;
    int Zodiac_Cipher[24][17] = {
                                     { 1, 2, 3, 4, 5, 4, 6, 7, 2, 8, 9,10,11,12,13,11, 7},
        {14,15,16,17,18,19,20,21, 1,22, 3,23,24,25,26,19,17},
        {27,28,19,29, 6,30, 8,31,26,32,33,34,35,19,36,37,38},
        {39,40, 4, 1, 2, 7, 3, 9,10,41, 6, 2,42,10,43,26,44},
        { 8,29,45,27, 5,28,46,47,48,12,20,22,15,14,17,31,19},
        {23,16,26,18,36, 1,24,30,38,21,26,13,49,37,50,39,40},
        {10,34,33,30,19,44,43, 9, 1,26,18, 7,32,21,39, 2, 7},
        {45,46, 4, 3, 2, 7,23,13,26,44,22,27, 6,29,10,10, 8},
        {51, 5,24,26,12,30,38,14,26,25,49,37,45,27,47, 1,52},
        { 7, 3,36,10,16,28,11,21,48,34,40,17,44, 6,22, 8,20},
        { 5,51,12, 9,15,14,30,37,16,33,45,38,43,29,10,21,22},
        {30, 1,36,10,53,32,19,47,48,46,17, 4,23,13,28,35,41},
        { 3,37,27,49,10, 6,33, 2,45,38,34,15,44,24,22,11,18},
        {47,30,25,28, 8,37, 1,49,45,27,43,34,41,38, 5,40, 3},
        {50, 6,12, 8,41, 1,52, 7,15,14,48,16,15,32,33, 9, 3},
        {29,11,39,47,43,42, 6,17,21,31,36,50,18, 2, 2,25,27},
        {34, 8,38,39,51,44, 4, 1, 2, 2, 5,42,41, 3,52, 7,15},
        {12,17,13,26,14,26,53,20,52,49,51,16,23, 1,41, 1, 7},
        { 2, 9,32,37,10, 6,51,16,53,46,19,26,53,29,39,26,14},
        {15, 5,17,18,19,24,44,53,32,19,41, 1, 2,52,45,33,53},
        {22,25,20, 7,13, 1,50,13,41,36,46,48,31,45,25,11,26},
        {53,17,46,52,52,21,17,37, 3, 9,10,13,35,20, 2,18,51},
        { 5,23,28,32,33,26,53,49,28,30,16,47, 7, 3,35,14,21},
        {15,44,13,47, 1,14,30,21,26,44,22,27,38,11,19,30, 8}
    };

    char fname[80];

    struct stepStruct *step;
    // srand(seed);

    if(argc !=8)
    {
oops:   fprintf(stderr, "\nUsage: %s filename startPos startChar maxChars maxIters seed\n\n", argv[0]);
        fprintf(stderr, "where filename == input file\n");
        fprintf(stderr, "      startPos == starting position for each line (numbered from 0)\n");
        fprintf(stderr, "      startChar == starting character in file (numbered from 0)\n");
        fprintf(stderr, "      maxChars == max characters to read (<= 0 to read all)\n");
        fprintf(stderr, "      maxIters == max iterations of re-estimation algorithm\n");
        fprintf(stderr, "      seed == seed value for pseudo-random number generator (PRNG)\n\n");
        fprintf(stderr, "For example:\n\n      %s datafile 0 0 0 100 1241\n\n", argv[0]);
        fprintf(stderr, "will read all of `datafile' and perform a maximum of 100 iterations.\n\n");
        fprintf(stderr, "For the English text example, try:\n\n      %s BrownCorpus 15 1000 50000 200 22761\n\n", argv[0]);
        fprintf(stderr, "will read from `BrownCorpus' and seed the PRNG with 22761,\n");
        fprintf(stderr, "will not read characters 0 thru 14 of each new line in `BrownCorpus',\n");
        fprintf(stderr, "will not save the first 1000 characters read and\n");
        fprintf(stderr, "will save a maximum of 50k observations.\n\n");
        exit(0);
    }

//    sprintf(fname, argv[1]);
    strcpy(fname, argv[1]);
    startPos = atoi(argv[2]);
    startChar = atoi(argv[3]);
    maxChars = atoi(argv[4]);
    maxIters = atoi(argv[5]);
    seed = atoi(argv[6]);
    No_Iter = atoi(argv[7]);

    ////////////////////////
    // read the data file //
    ////////////////////////
    
    // determine number of observations
    printf("GetT... ");
    fflush(stdout);
    T = GetT(fname, 
             startPos,
             startChar,
             maxChars);
    
    printf("T = %d\n", T);

    // allocate memory
    printf("allocating %lu bytes of memory... ", (T + 1) * sizeof(struct stepStruct));
    fflush(stdout);
    if((step = calloc(T + 1, sizeof(struct stepStruct))) == NULL)
    {
        fprintf(stderr, "\nUnable to allocate alpha\n\n");
        exit(0);
    }
    // printf("done\n");

    // read in the observations from file
    // printf("GetObservations... ");
    fflush(stdout);
    T = GetObservations(fname, 
                        step, 
                        T,
                        startPos,
                        startChar,
                        maxChars,
                        0);

    printf("T = %d\n", T);


    digraph_freq_A(step, A, T);

    // printf("initial digraph A =\n");
    // printA(A);

    T = 408;

    printf("allocating %lu bytes of memory... ", (T + 1) * sizeof(struct stepStruct));
    fflush(stdout);
    if((step = calloc(T + 1, sizeof(struct stepStruct))) == NULL)
    {
        fprintf(stderr, "\nUnable to allocate alpha\n\n");
        exit(0);
    }
    printf("done\n");
   
    for(int a = 0; a < 24; a++)
    {
        for(int b = 0; b < 17; b++)
        {
            step[cnt++].obs = Zodiac_Cipher[a][b];
        }
    }
    
    // printf("The Zodia Cipher is:\n");

    // for(int a = 0; a < 24; a++)
    // {
    //     for(int b = 0; b < 17; b++)
    //     {
    //         printf("%d ",Zodiac_Cipher[a][b]);
    //     }
    //     printf("\n");
    // }
    
  
    srandom(seed);
    int max_match = 0;
    int max_hmm_model = 0;
    int model_no = 0;
    double score = 0.0,
            max_score = 0.0;
    
    for(int temp_iter = 0; temp_iter < No_Iter; temp_iter++)
    {
            

    // initialize pi[] and B[][]
    initMatrices(pi,B,seed);


    // printf("Value of T after cipherText read is T = %d\n", T);

   
    
    // print pi[], A[][] and B[][] transpose
    // printf("\nN = %d, M = %d, T = %d\n", N, M, T);
    // printf("initial pi =\n");
    // printPi(pi);
    // printf("initial A =\n");
    // printA(A);
    // printf("initial B^T =\n");
    // printBT(B);

    // initialization
    iter = 0;
    logProb = -1.0;
    newLogProb = 0.0;

    // main loop
    while(iter < maxIters)
    {

        // printf("\nbegin iteration = %d\n", iter);

        logProb = newLogProb;

        // alpha (or forward) pass
        //printf("alpha pass... ");
        fflush(stdout);
        alphaPass(step, pi, A, B, T);
        //printf("done\n");
        
        // beta (or backwards) pass
        //printf("beta pass... ");
        fflush(stdout);
        betaPass(step, pi, A, B, T);
        //printf("done\n");
        
        // compute gamma's and diGamma's
        //printf("compute gamma's and diGamma's... ");
        fflush(stdout);
        computeGammas(step, pi, A, B, T);
        //printf("done\n");
        
        // find piBar, reestimate of pi
        // printf("reestimate pi... ");
        fflush(stdout);
        reestimatePi(step, piBar);
       // printf("done\n");
        
        // find Abar, reestimate of A
        // printf("reestimate A... ");
        //fflush(stdout);
       // reestimateA(step, Abar, T);  We do not have to reestimate A matrix, hence commenting the call
        // printf("done\n");
        
        // find Bbar, reestimate of B
        // printf("reestimate B... ");
        fflush(stdout);
        reestimateB(step, Bbar, T);
        //printf("done\n");
        
/*#ifdef PRINT_REESTIMATES
        printf("piBar =\n");
        printPi(piBar);
        printf("Abar =\n");
        printA(Abar);
        printf("Bbar^T = \n");
        printBT(Bbar);
#endif // PRINT_REESTIMATES*/

        // assign pi, A and B corresponding "bar" values
        for(i = 0; i < N; ++i)
        {
            pi[i] = piBar[i];
        
            // for(j = 0; j < N; ++j)
            // {
            //     A[i][j] = Abar[i][j];
            // }

            for(j = 0; j < M; ++j)
            {
                B[i][j] = Bbar[i][j];
            }
            
        }// next i

        // compute log [P(observations | lambda)], where lambda = (A,B,pi)
        newLogProb = 0.0;
        for(i = 0; i < T; ++i)
        {
            newLogProb += log(step[i].c);
        }
        newLogProb = -newLogProb;

        // a little trick so that no initial logProb is required
        if(iter == 0)
        {
            logProb = newLogProb - 1.0;
        }

        // printf("completed iteration = %d, log [P(observation | lambda)] = %f\n", 
        //         iter, newLogProb);

        ++iter;

    }// end while
    
    // printf("\nT = %d, N = %d, M = %d, iterations = %d\n\n", T, N, M, iter);
    // printf("final pi =\n");
    // printPi(pi);
    // printf("\nfinal A =\n");
    // printA(A);
    // printf("\nfinal B^T =\n");
    // printBT(B);
    // printf("\nlog [P(observations | lambda)] = %f\n\n", newLogProb);

    int count = 0;
      for(int i = 0 ; i < N ; i++){
            double max = B[i][0];
            int elem = 0;
            for(int j = 1 ; j < M ; j++){
                if(B[i][j] > max){
                    max = B[i][j];
                    elem = j;
                }
            }
            // printf("\nThe %c element is mapped to %c", i+97, elem+97);

            if(elem % 54 == i)
                count++;

    }

    // printf("\n\n\nThe key is %d", randSubs);

    // printf("\nFinal Count value: %d", count);

    score = (double) count / N;

    // printf("\nThe percentage is %0.4f\n\n", score);

    if (count > max_match)
    {
        max_match = count;
        max_hmm_model = temp_iter;
    }

    printf("\nModel %d has score %8.5f", temp_iter + 1, score); 
    if(score > max_score)
        max_score = score;


    }

    printf("\nThe best score was for model %d with score %8.5f\n\n", max_hmm_model + 1, max_score);

}// end hmm

void digraph_freq_A(struct stepStruct *step, double A[][N], int T)
{
    double sum_row = 0.0;
    double row_sto = 0.0;
    for(int i = 0;i < 26; i++)
         for(int j = 0;j < 26;j++)
            A[i][j] = 0.0;
    
    for(int t =1;t < T;t++)
        A[step[t-1].obs][step[t].obs]++;

    for(int i = 0;i < 26; i++)
         for(int j = 0;j < 26;j++)
            A[i][j] += 5.0;

    for(int i = 0;i < 26; i++)
    {
        sum_row = 0.0;
         for(int j = 0;j < 26;j++)
            sum_row += A[i][j];
        
        for(int j = 0;j < 26;j++)
        {
            A[i][j] /= sum_row;
            A[i][j] += row_sto;
        }

        //  printf("\n The sum for row %d is %f", i+1, row_sto);
    }
    // printf("\n\n");
}
    

//
// alpha pass (or forward pass) including scaling
//
void alphaPass(struct stepStruct *step,
               double pi[], 
               double A[][N],
               double B[][M],
               int T)
{
    int i,
        j,
        t;
        
    double ftemp;
    
    // compute alpha[0]'s
    ftemp = 0.0;
    for(i = 0; i < N; ++i)
    {
        step[0].alpha[i] = pi[i] * B[i][step[0].obs];
        ftemp += step[0].alpha[i];
    }
    step[0].c = 1.0 / ftemp;

    // scale alpha[0]'s
    for(i = 0; i < N; ++i)
    {
        step[0].alpha[i] /= ftemp;
    }

    // alpha pass
    for(t = 1; t < T; ++t)
    {
        ftemp = 0.0;
        for(i = 0; i < N; ++i)
        {
            step[t].alpha[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                step[t].alpha[i] += step[t - 1].alpha[j] * A[j][i];
            }
            step[t].alpha[i] *= B[i][step[t].obs];
            ftemp += step[t].alpha[i];
        }
        step[t].c = 1.0 / ftemp;
        
        // scale alpha's
        for(i = 0; i < N; ++i)
        {
            step[t].alpha[i] /= ftemp;
        }
    
    }// next t
    
}// end alphaPass


//
// beta pass (or backwards pass) including scaling
//
void betaPass(struct stepStruct *step,
              double pi[], 
              double A[][N],
              double B[][M],
              int T)
{
    int i,
        j,
        t;

    // compute scaled beta[T - 1]'s
    for(i = 0; i < N; ++i)
    {
        step[T - 1].beta[i] = 1.0 * step[T - 1].c;
    }

    // beta pass
    for(t = T - 2; t >= 0; --t)
    {
        for(i = 0; i < N; ++i)
        {
            step[t].beta[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                step[t].beta[i] += A[i][j] * B[j][step[t + 1].obs] * step[t + 1].beta[j];
            }
            
            // scale beta's (same scale factor as alpha's)
            step[t].beta[i] *= step[t].c;
        }

    }// next t
        
}// end betaPass


//
// compute gamma's and diGamma's including optional error checking
//
void computeGammas(struct stepStruct *step,
                   double pi[], 
                   double A[][N],
                   double B[][M],
                   int T)
{
    int i,
        j,
        t;
        
    double denom;

#ifdef CHECK_GAMMAS
    double ftemp,
           ftemp2;
#endif // CHECK_GAMMAS

    // compute gamma's and diGamma's
    for(t = 0; t < T - 1; ++t)// t = 0,1,2,...,T-2
    {
        
#ifdef CHECK_GAMMAS
        ftemp2 = 0.0;
#endif // CHECK_GAMMAS

        for(i = 0; i < N; ++i)
        {
            step[t].gamma[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                step[t].diGamma[i][j] = (step[t].alpha[i] * A[i][j] * B[j][step[t + 1].obs] * step[t + 1].beta[j]);
                step[t].gamma[i] += step[t].diGamma[i][j];
            }

            #ifdef CHECK_GAMMAS
             // verify that gamma[i] == alpha[i]*beta[i] / sum(alpha[j]*beta[j])
            ftemp2 += step[t].gamma[i];
            ftemp = 0.0;
            for(j = 0; j < N; ++j)
            {
                ftemp += step[t].alpha[j] * step[t].beta[j];
            }
            ftemp = (step[t].alpha[i] * step[t].beta[i]) / ftemp;
            if(DABS(ftemp - step[t].gamma[i]) > EPSILON)
            {
                printf("gamma[%d] = %f (%f) ", i, step[t].gamma[i], ftemp);
                printf("********** Error !!!\n");
            }
#endif 
// CHECK_GAMMAS

        }// next i
            
#ifdef CHECK_GAMMAS
        if(DABS(1.0 - ftemp2) > EPSILON)
        {
            printf("sum of gamma's = %f (should sum to 1.0)\n", ftemp2);
        }
#endif // CHECK_GAMMAS
            
    }// next t
    
    // special case for t = T-1
    for(j = 0; j < N; ++j)  
    {
        step[T-1].gamma[j] = step[T-1].alpha[j];
    }

}// end computeGammas


//
// reestimate pi, the initial distribution
//
void reestimatePi(struct stepStruct *step, 
                  double piBar[])
{
    int i;
    
    // reestimate pi[]        
    for(i = 0; i < N; ++i)
    {
        piBar[i] = step[0].gamma[i];
    }
        
}// end reestimatePi


//
// reestimate the A matrix
//
void reestimateA(struct stepStruct *step, 
                 double Abar[][N], 
                 int T)
{
    int i,
        j,
        t;
    
    double numer,
           denom;
           
    // reestimate A[][]
    for(i = 0; i < N; ++i)
    {
        denom = 0.0;
        // t = 0,1,2,...,T-2
        for(t = 0; t < T - 1; ++t)
        {
            denom += step[t].gamma[i];
            
        }// next t

        for(j = 0; j < N; ++j)
        {
            numer = 0.0;

            // t = 0,1,2,...,T-2
            for(t = 0; t < T - 1; ++t)
            {
                numer += step[t].diGamma[i][j];
            }// next t

            Abar[i][j] = numer / denom;
        
        }// next j
        
    }// next i
        
} // end reestimateA    


//
// reestimate the B matrix
//
void reestimateB(struct stepStruct *step, 
                 double Bbar[][M], 
                 int T)
{
    int i,
        j,
        t;
    
    double numer,
           denom;
           
    // reestimate B[][]
    for(i = 0; i < N; ++i)
    {
        denom = 0.0;
        // t = 0,1,2,...,T-1
        for(t = 0; t < T; ++t)
        {
            denom += step[t].gamma[i];
            
        }// next t

        for(j = 0; j < M; ++j)
        {
            numer = 0.0;

            // t = 0,1,2,...,T-1
            for(t = 0; t < T; ++t)
            {
                if(step[t].obs == j)
                {
                    numer += step[t].gamma[i];
                }
                
            }// next t

            Bbar[i][j] = numer / denom;
        
        }// next j
        
    }// next i
        
}// end reestimateB


//
// initialize pi[], A[][] and B[][]
//
void initMatrices(double pi[],
                  double B[][M],
                  int seed)
{
    int i,
        j;
        
    double prob,
           ftemp,
           ftemp2;
    
    // initialize pseudo-random number generator
    // srandom(seed);

    // initialize pi
    prob = 1.0 / (double)N;
    ftemp = prob / 10.0;
    ftemp2 = 0.0;
    for(i = 0; i < N; ++i)
    {
        if((rand() & 0x1) == 0)
        {
            pi[i] = prob + (double)(random() & 0x7) / 8.0 * ftemp;
        }
        else
        {
            pi[i] = prob - (double)(random() & 0x7) / 8.0 * ftemp;
        }
        ftemp2 += pi[i];
        
    }// next i
    
    for(i = 0; i < N; ++i)
    {
        pi[i] /= ftemp2;
    }

    // initialize B[][]
    prob = 1.0 / (double)M;
    ftemp = prob / 10.0;
    for(i = 0; i < N; ++i)
    {
        ftemp2 = 0.0;
        for(j = 0; j < M; ++j)
        {
            if((random() & 0x1) == 0)
            {
                B[i][j] = prob + (double)(random() & 0x7) / 8.0 * ftemp;
            }
            else
            {
                B[i][j] = prob - (double)(random() & 0x7) / 8.0 * ftemp;
            }
            ftemp2 += B[i][j];
            
        }// next j
        
        for(j = 0; j < M; ++j)
        {
            B[i][j] /= ftemp2;
        }
        
    }// next i
    
}// end initMatrices


//
// read (but don't save) observations get T
//
int GetT(char fname[],
         int startPos,
         int startChar,
         int maxChars)
{
    FILE *in;

    int T,
        i,
        j,
        len,
        thisStartPos,
        totalNum,
        num;
        
    char temp[MAX_CHARS + 1];
    
    char space[1] = {" "};
    
    char alphabet[M] = ALPHABET;

    in = fopen(fname, "r");
    if(in == NULL)
    {
        fprintf(stderr, "\nError opening file %s\n\n", fname);
        exit(0);
    }
    
#ifdef PRINT_GET_T
    printf("\n");
#endif // PRINT_GET_T

    // count 'em
    totalNum = num = 0;
    while(fgets(temp, MAX_CHARS, in) != NULL)
    {
        len = strlen(temp);

        // each line should end with a single space
        while((strncmp(&temp[len - 1], space, 1) == 0) && (len > 0))
        {
            --len;
        }

        strncpy(&temp[len], space, 1);
        
        thisStartPos = startPos;
        
        // ignore leading spaces
        while((strncmp(&temp[thisStartPos], space, 1) == 0) && (thisStartPos < len))
        {
            ++thisStartPos;
        }
        
        for(i = thisStartPos; i <= len; ++i)
        {
            // find alphabetic characters, ignoring case
            // also drop all non-alphabet characters as well as space
            for(j = 0; j < M; ++j)
            {
                if(strncasecmp(&temp[i], &alphabet[j], 1) == 0)
                {
                    ++totalNum;
                    if(totalNum >= startChar)
                    {
#ifdef PRINT_GET_T
                        printf("%c %d\n", alphabet[j], num);
#endif // PRINT_GET_T
                        ++num;
                        if((maxChars > 0) && (num >= maxChars))
                        {
                            return(num);
                        }
                        
                    }// end if
                    
                    break;
                    
                }// end if
                
            }// next j

        }// next i
                
    }// end while
    
    fclose(in);

    return(num);

}// end GetT


//
// read and save observations

int GetObservations(char fname[], 
                    struct stepStruct *step,
                    int T,
                    int startPos,
                    int startChar,
                    int maxChars,int flag)
{
    FILE *in;

    int i,
        j,
        len,
        num,
        thisStartPos,
        totalNum;
   
        
    char temp[MAX_CHARS + 2];
    
    char space[1] = {" "};
    
    char alphabet[M] = ALPHABET;

    randSubs = 0;

    if(flag == 1)
    {
        while(randSubs == 0)
        {
            randSubs = rand() % 10;
        }
    }
    
    in = fopen(fname, "r");
    if(in == NULL)
    {
        fprintf(stderr, "\nError opening file %s\n\n", fname);
        exit(0);
    }
    
#ifdef PRINT_OBS
    printf("\n");
#endif // PRINT_OBS

    // read 'em in
    totalNum = num = 0;
    while(fgets(temp, MAX_CHARS, in) != NULL)
    {
        len = strlen(temp);
        
        // each line should end with a single space
        while((strncmp(&temp[len - 1], space, 1) == 0) && (len > 0))
        {
            --len;
        }
        strncpy(&temp[len], space, 1);
        
        thisStartPos = startPos;
        
        // ignore leading spaces
        while((strncmp(&temp[thisStartPos], space, 1) == 0) && (thisStartPos < len))
        {
            ++thisStartPos;
        }

        for(i = thisStartPos; i <= len; ++i)
        {
            // find alphabetic characters, ignoring case,
            // drop all non-alphabet characters including space as ALPHABET does not have space 
            for(j = 0; j < M; ++j)
            {
                if(strncasecmp(&temp[i], &alphabet[j], 1) == 0)
                {
                    ++totalNum;
                    if(totalNum >= startChar)
                    {
                        if(flag == 1)
                        {
                            step[num].obs = (j + randSubs) % 26;
                        }
                        else
                        {
                                step[num].obs = j;
                        }
                        
                        
#ifdef PRINT_OBS
                        printf("%c %d\n", alphabet[j], step[num].obs);
#endif // PRINT_OBS
                        ++num;
                        if(num > T)
                        {
                            printf("\nError --- T exceeded in GetObservations()\n\n");
                            exit(0);
                        }
                        if((maxChars > 0) && (num >= maxChars))
                        {
                            return(num);
                        }
                        
                    }// end if
                    
                    break;
                    
                }// end if
                
            }// next j

        }// next i

    }// end while
    
    fclose(in);

    return(num);

}// end GetObservations


//
// print pi[]
//
void printPi(double pi[])
{
    int i;
        
    double ftemp;

    ftemp = 0.0;
    for(i = 0; i < N; ++i)
    {
        printf("%8.5f ", pi[i]);
        ftemp += pi[i];
    }
    printf(",  sum = %f\n", ftemp);

}// end printPi


//
// print A[][]
//
void printA(double A[][N])
{
    int i,
        j;
        
    double ftemp;

    for(i = 0; i < N; ++i)
    {
        ftemp = 0.0;
        for(j = 0; j < N; ++j)
        {
            printf("%8.5f ", A[i][j]);
            ftemp += A[i][j];
        }
        printf(",  sum = %f\n", ftemp);
        
    }// next i

}// end printA


//
// print BT[][]
//
void printBT(double B[][M])
{
    int i,
        j;
        
    double ftemp;
    
    char alphabet[M] = ALPHABET;
    
    for(i = 0; i < M; ++i)
    {
        printf("%c ", alphabet[i]);
        for(j = 0; j < N; ++j)
        {
            printf("%8.5f ", B[j][i]);
        }
        printf("\n");
    }
    for(i = 0; i < N; ++i)
    {
        ftemp = 0.0;
        for(j = 0; j < M; ++j)
        {
            ftemp += B[i][j];
        }
        printf("sum[%d] = %f ", i, ftemp);
    }
    printf("\n");

}// end printB

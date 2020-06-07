//
//  main.c
//  HMM
//
//  Created by Yash Sahasrabuddhe on 2/3/20.
//  Copyright © 2020 Yash Sahasrabuddhe. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <string.h>
#include <stdbool.h>

int minIters, iters;
float old_log_prob, improve_threshold;

void generateMatrix(int N, int M, float *matrix[])
{
    int sum=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
            sum += matrix[i][j]=rand();
        for(int j=0;j<M;j++)
            matrix[i][j] /= sum;
        sum = 0;
    }
}

void generateObsSeq(int T, int ObsSeq[T])
{
    
    char ch, *path="corpus.dos/",file_path[30]="\0";
    int t=0;
    FILE *fp;
    DIR *d;
    struct dirent *file;
    d = opendir("corpus.dos");
    while((file = readdir(d)) && t<T)
    {
        if(strcmp(file->d_name,".") && strcmp(file->d_name,".."))
        {
            strcat(file_path,path);
            strcat(file_path,file->d_name);
            fp = fopen(file_path,"r");
            while(!feof(fp) && t<T)
            {
                ch = fgetc(fp);
                if(ch>=97 && ch<=122)
                    ObsSeq[t++] = (int)ch%97;
                else if(ch>=65 && ch<=90)
                    ObsSeq[t++] = ((int)ch+32)%97;
                else if(ch==32)
                    ObsSeq[t++] = 26;
                
            }
            memset(file_path,0,strlen(file_path));
            fclose(fp);
        }
    }
}

void alphaPass(float *A[],float *B[],float *Pi[],int T,int ObsSeq[T],int N,int M,float *c,float *alpha[])
{
    
    c[0] = 0;
    for(int i=0;i<N;i++)
    {
        alpha[i][0] = Pi[0][i]*B[i][ObsSeq[0]];  //computing alpha(0)
        c[0] = c[0] + alpha[i][0];
    }
    
    c[0] = 1/c[0];
    for(int i=0;i<N;i++)
        alpha[i][0] = c[0]*alpha[i][0];    //compute alpha0(i)
    
    for(int t=1;t<T;t++)
    {
        c[t] = 0;
        for(int i=0;i<N;i++)
        {
            alpha[i][t] = 0;
            for(int j=0;j<N;j++)
                alpha[i][t] += (alpha[j][t-1]*A[j][i]);
            alpha[i][t] *= B[i][ObsSeq[t]];
            c[t] += alpha[i][t];
        }
        c[T] = 1/c[T];      //scaling alpha[i][t]
        for(int i=0;i<N;i++)
            alpha[i][t] = c[t]*alpha[i][t];
    }
}

void betaPass(float *A[],float *B[],float *Pi[],int T,int ObsSeq[T],int N,int M,float *c,float *alpha[],float *beta[])
{
    for(int i=0;i<N;i++)
        beta[i][T-1] = c[T-1];
    
    for(int t=T-2;t>=0;t--)
    {
        for(int i=0;i<N;i++)
        {
            beta[i][t] = 0;
            for(int j=0;j<N;j++)
                beta[i][t] = beta[i][t]+((A[i][j]*B[j][ObsSeq[t+1]])*beta[j][t+1]);
            beta[i][t] = c[t]*beta[i][t];
        }
    }
}

void gamma_digammaPass(float *A[],float *B[],float *Pi[],int T,int ObsSeq[T],int N,int M,float *c,float *alpha[],float *beta[],float *gamma[], float ** digamma[])
{
    float denom;
    for(int t=0;t<T-1;t++)
    {
        denom = 0.0;
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
                denom += ((alpha[i][t])*(A[i][j]*B[j][ObsSeq[t+1]])*(beta[j][t+1]));
        }
        for(int i=0;i<N;i++)
        {
            gamma[i][t] = 0;
            for(int j=0;j<N;j++)
            {
                digamma[i][j][t] = ((alpha[i][t])*(A[i][j]*B[j][ObsSeq[t+1]])*(beta[j][t+1]))/denom;
                gamma[i][t] += digamma[i][j][t];
            }
        }
    }
    
    denom = 0.0;
    for(int i=0;i<N;i++)
        denom += alpha[i][T-1];
    for(int i=0;i<N;i++)
        gamma[i][T-1] = alpha[i][T-1]/denom;
}

void estimate(float *A[],float *B[],float *Pi[],int T,int ObsSeq[T],int N,int M,float *c,float *alpha[],float *beta[],float *gamma[], float ** digamma[])
{
    float numer = 0.0,
    denom = 0.0;
    for(int i=0;i<N;i++)
        Pi[0][i] = gamma[i][0];
    
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            numer = 0.0;
            denom = 0.0;
            for(int t=0;t<T-1;t++)
            {
                numer += digamma[i][j][t];
                denom += gamma[i][t];
            }
            A[i][j] = numer/denom;
        }
    }
    
    
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            numer = 0.0;
            denom = 0.0;
            for(int t=0;t<T;t++)
            {
                if(ObsSeq[t]==j)
                    numer += gamma[i][t];
                denom += gamma[i][t];
            }
            
            B[i][j] = numer/denom;
        }
    }
}

void computeProb(int T, float *c, float logProb)
{
    logProb = 0;
    for(int i=0;i<T;i++)
        logProb += log(c[i]);
    logProb = -logProb;
}

bool computeCheck(float logProb)
{
    bool result = false;
    iters += 1;
    float delta = fabs(logProb - old_log_prob);
    if(iters<minIters || delta>improve_threshold)
    {
        old_log_prob = logProb;
        result = true;
    }
    
    return result;
}


int main(int argc, const char * argv[]) {
    
    iters = 0;
    minIters = 20;
    improve_threshold = 0.01;
    old_log_prob = -INFINITY;
    int N,
    M=27,
    T=50000;
    
    int ObsSeq[T];
    
    if(argc>1)
        N = atoi(argv[1]);
    else
        N=2; //default case
    
    float *A[N], *B[N], *Pi[1], *alpha[N], *beta[N], *gamma[N], **digamma[N], *c, logProb = 0;
    
    Pi[0] = (float *)malloc(N * sizeof(float));
    c = (float *)malloc(T * sizeof(float));
    
    for(int i=0;i<N;i++)
    {
        
        A[i] = (float *)malloc(N * sizeof(float));
        B[i] = (float *)malloc(M * sizeof(float));
        alpha[i] = (float *)malloc(T * sizeof(float));
        beta[i] = (float *)malloc(T * sizeof(float));
        gamma[i] = (float *)malloc(T * sizeof(float));
        digamma[i] = (float **)malloc(N * sizeof(float));
        for (int j=0;j<N;j++)
            digamma[i][j] = (float *)malloc(T * sizeof(float));
        
    }
    
    Pi[0] = (float *)malloc(N * sizeof(float));
    c = (float *)malloc(T * sizeof(float));
    
    generateMatrix(N,N,A);      //generate stochastic A matrix
    generateMatrix(N,M,B);      //generate stochastic B matrix
    generateMatrix(1,N,Pi);     //generate stochastic Pi matrix
    generateObsSeq(T,ObsSeq);  //generate Observation Sequence
    
    
    printf("\nInitial A Matrix");
    for(int i=0;i<N;i++)
    {
        printf("\n");
        for(int j=0;j<N;j++)
        {
            printf("%f\t",A[i][j]);
        }
    }
    
    printf("\nInitial B Matrix\n");
    for(int i=0;i<M;i++)
    {
        printf("\n");
        for(int j=0;j<N;j++)
        {
            printf("%f\t",B[j][i]);
        }
    }
    
    printf("\nInitial Pi Matrix\n");
    for(int i=0;i<1;i++)
    {
        printf("\n");
        for(int j=0;j<N;j++)
        {
            printf("%f\t",Pi[i][j]);
        }
    }
    
    bool flag = true;
    
    while(flag)
    {
        alphaPass(A,B,Pi,T,ObsSeq,N,M,c,alpha);
        betaPass(A,B,Pi,T,ObsSeq,N,M,c,alpha,beta);
        gamma_digammaPass(A,B,Pi,T,ObsSeq,N,M,c,alpha,beta,gamma,digamma);
        estimate(A,B,Pi,T,ObsSeq,N,M,c,alpha,beta,gamma,digamma);
        computeProb(T,c,logProb);
        flag = computeCheck(logProb);
    }
    
    printf("\nFinal A Matrix\n");
    for(int i=0;i<N;i++)
    {
        printf("\n");
        for(int j=0;j<N;j++)
        {
            printf("%f\t",A[i][j]);
        }
    }
    
    printf("\nFinal B Matrix\n");
    for(int i=0;i<M;i++)
    {
        printf("\n");
        for(int j=0;j<N;j++)
        {
            printf("%f\t",B[j][i]);
        }
        printf("%c",97+i);
    }
    
    printf("\nFinal Pi Matrix\n");
    for(int i=0;i<1;i++)
    {
        printf("\n");
        for(int j=0;j<N;j++)
        {
            printf("%f\t",Pi[i][j]);
        }
    }
    
    free(A);free(B);free(alpha);free(beta);free(gamma);free(digamma);free(Pi[0]);free(c);
    
    return 0;
}





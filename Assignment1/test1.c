#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <string.h>

int iters,min_iters;
float old_log_prob, threshold;

int random_stochastic(int n, int m, float *mat[]) //n and m are the dimensions of the matrix
{
    int sum = 0, i, j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<m;j++)
            sum += mat[i][j] = rand();
        for(j=0;j<m;j++)
            mat[i][j] /= sum;
        sum = 0;
    }
    return 0;
}

int get_obs_seq(int T,int o[T])
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
                    o[t++] = (int)ch%97;
                else if(ch>=65 && ch<=90)
                    o[t++] = ((int)ch+32)%97;
                else if(ch==32)
                    o[t++] = 26;
                
            }
            memset(file_path,0,strlen(file_path));
            fclose(fp);
        }
    }
    return 0;
}

int alpha_pass(int N,int M,int T,float *A[],float *B[],float *pi[],int O[T],float *c,float *alpha[])
{
    int i,j,t;
    c[0] = 0;
    for(i=0;i<N;i++)
    {
        alpha[i][0] = pi[0][i] * B[i][O[0]];
        c[0] = c[0]+alpha[i][0];
    }
    c[0] = 1/c[0];
    for(i=0;i<N;i++)
        alpha[i][0] = c[0] * alpha[i][0];
    for(t=1;t<T;t++)
    {
        c[t] = 0;
        for(i=0;i<N;i++)
        {
            alpha[i][t] = 0;
            for(j=0;j<N;j++)
                alpha[i][t] = alpha[i][t] + (alpha[j][t-1] * A[j][i]);
            alpha[i][t] = alpha[i][t] * B[i][O[t]];
            c[t] = c[t]+alpha[i][t];
        }
        c[t] = 1/c[t];
        for(i=0;i<N;i++)
            alpha[i][t] = c[t]*alpha[i][t];
    }
    //for(i=0;i<T;i++)
    //    printf("%f\t",c[i]);
    return 0;
}

int beta_pass(int N,int M,int T,float *A[],float *B[],float *pi[],int O[T],float *c,float *alpha[],float *beta[])
{
    int i,j,t;
    for(i=0;i<N;i++)
    {
        beta[i][T-1] = c[T-1];
        //printf("beta%f\t",beta[i][T-1]);
    }
    
    for(t=T-2;t>-1;t--)
    {
        for(i=0;i<N;i++)
        {
            beta[i][t] = 0;
            for(j=0;j<N;j++)
            {
                beta[i][t] = beta[i][t] + ((A[i][j] * B[j][O[t+1]]) * beta[j][t+1]);
                //printf("beta-%f\n",beta[i][t]);
            }
            beta[i][t] = beta[i][t] * c[t];
        }
    }
    return 0;
    
}

int gamma_digamma(int N,int M,int T,float *A[],float *B[],float *pi[],int O[T],float *c,float *alpha[],float *beta[],float *gamma[], float **digamma[])
{
    int i,j,t;
    float denom;
    for(t=0;t<T-1;t++)
    {
        denom = 0;
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
                denom = denom+((alpha[i][t]) * (A[i][j]*B[j][O[t+1]]) * (beta[j][t+1]));
        }
        
        
        for(i=0;i<N;i++)
        {
            gamma[i][t] = 0;
            for(j=0;j<N;j++)
            {
                digamma[i][j][t] = (alpha[i][t] * A[i][j] * B[j][O[t+1]] * beta[j][t+1])/denom;
                gamma[i][t] += digamma[i][j][t];
            }
        }
    }
    denom = 0;
    for(i=0;i<N;i++)
        denom += alpha[i][T-1];
    return 0;
    
}

int estimate(int N,int M,int T,float *A[],float *B[],float *pi[],int O[T],float *c,float *alpha[],float *beta[],float *gamma[], float **digamma[])
{
    int i,j,t;
    float numer,denom;
    
    //pi
    for(i=0;i<N;i++)
        pi[0][i] = gamma[i][0];
    
    //A
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            numer = denom = 0;
            for(t=0;t<T-1;t++)
            {
                numer = numer+digamma[i][j][t];
                denom = denom+gamma[i][t];
            }
            A[i][j] = numer/denom;
        }
    }
    
    //B
    for(i=0;i<N;i++)
    {
        for(j=0;j<M;j++)
        {
            numer = denom = 0;
            for(t=0;t<T;t++)
            {
                if(O[t] == j)
                    numer = numer+gamma[i][t];
                denom = denom+gamma[i][t];
            }
            B[i][j] = numer/denom;
        }
    }
    return 0;
}

int probability(int T,float *c,float log_prob)
{
    int i;
    for(i=0;i<T;i++)
        log_prob = log_prob+log(c[i]);
    log_prob = -log_prob;
    return 0;
}

int decision(float log_prob)
{
    iters=iters+1;
    int flag=0;
    float delta = fabs(log_prob-old_log_prob);
    if(iters<min_iters || delta>=threshold)
    {
        old_log_prob=log_prob;
        flag = 1;
    }
    return flag;
}

int main(int argc, char *argv[])
{
    iters = 0;
    min_iters = 200;
    threshold = 0.001;
    old_log_prob=-INFINITY;
    
    
    int N,M=27,T=50000,i,j;
    int Obs_seq[T];
    //get N as command line argument
    if(argc>1)
        N = atoi(argv[1]);
    else
        N = 2;
    //initialize A,B and pi matrices
    float *A[N], *B[N], *pi[1], *alpha[N], *beta[N], *c, *gamma[N], **digamma[N],log_prob=0;
    for (i=0; i<N; i++)
    {
        A[i] = (float *)malloc(N * sizeof(float));
        B[i] = (float *)malloc(M * sizeof(float));
        alpha[i] = (float *)malloc(T * sizeof(float));
        beta[i] = (float *)malloc(T * sizeof(float));
        gamma[i] = (float *)malloc(T * sizeof(float));
        digamma[i] = (float **)malloc(N * sizeof(float));
        for(j=0;j<N;j++)
            digamma[i][j] = (float *)malloc(T * sizeof(float));
    }
    pi[0] = (float *)malloc(N * sizeof(float));
    c = (float *)malloc(T * sizeof(float));
    
    random_stochastic(N,N,A);
    random_stochastic(N,M,B);
    random_stochastic(1,N,pi);
    get_obs_seq(T,Obs_seq);
    
    printf("\nINITIAL RANDOMIZED ROW-STOCHASTIC MATRICES\n");
    printf("\nState transition matrix(A):\n");
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("%f\t",A[i][j]);
        }
        printf("\n");
    }
    printf("\nTranspose of Observation matrix(B):\n");
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("%f\t",B[j][i]);
        }
        printf("\n");
    }
    printf("\nInitial State Distribution matrix(pi)\n");
    for(i=0;i<1;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("%f\t",pi[i][j]);
        }
        printf("\n");
    }
    int flag=1;
    printf("\nPatience is Key!!!\n");
    while(flag)
    {
        alpha_pass(N,M,T,A,B,pi,Obs_seq,c,alpha);
        beta_pass(N,M,T,A,B,pi,Obs_seq,c,alpha,beta);
        gamma_digamma(N,M,T,A,B,pi,Obs_seq,c,alpha,beta,gamma,digamma);
        estimate(N,M,T,A,B,pi,Obs_seq,c,alpha,beta,gamma,digamma);
        probability(T,c,log_prob);
        flag = decision(log_prob);
    }
    printf("\nFinally after %d iterations....\n",iters);
    printf("\nState transition matrix(A):\n");
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("%f\t",A[i][j]);
        }
        printf("\n");
    }
    printf("\nTranspose of Observation matrix(B):\n");
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("%f\t",B[j][i]);
        }
        printf("%c",97+i);
        printf("\n");
    }
    printf("\nInitial State Distribution matrix(pi)\n");
    for(i=0;i<1;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("%f\t",pi[i][j]);
        }
        printf("\n");
    }
    free(A);free(B);free(alpha);free(beta);free(gamma);free(digamma);free(pi[0]);free(c);
    return 0;
}






#include "libraries_definitions.h"

using namespace std;

int gaussPivot (int, double **);                         
void substitucio (double **, int );
void canviarFiles (double **, int, int);

int gaussPivot (int n, double **mat) {                                  /* Retorna el nombre de canvis de files */
    int j, i, k, fimax, comptador=0;
    double m=0, max;                                                    /* m és el multiplicador */

    for(j=0; j<n-1; j++) {                                              /* Ho fem per cada columna */
        max=fabs(mat[j][j]);                                            /* Inicialitzem maxim a l'element de la diagonal */
        fimax=j;                                            
        for(k=j; k<n; k++)                                              /* Busquem el nou màxim i la fila */
            if(fabs(mat[k][j])>max) {
                max=fabs(mat[k][j]);
                fimax=k;
            }

        if(fimax!=j)                                                    /* Mirem si s'ha canviat alguna fila */
            comptador++;

        canviarFiles(mat, j, fimax);

        for(i=j+1; i<n; i++) {
            m=mat[i][j]/mat[j][j];                                      /* Trobem el multiplicador de cada coeficient */
            for(k=j; k<=n; k++) {                                       /* Modifiquem la matriu */
                mat[i][k]-= (m * mat[j][k]);    
                if (fabs(mat[i][k])<1.e-8)                              /* Si el nombre és molt petit posem un zero per estètica */
                    mat[i][k]=0;    
            }
        } 
    }
    return comptador;
}

//Function that recives a triangulated matrix and returns the solution of the system
void substitucio (double **matU, int n) {                   
    int i, j;
    double sum;
    for (i=n-1; i>=0; i--) {
        sum=0;
        for(j=n-1; j>i; j--) {
            sum+=(matU[j][i]*matU[n][j]);                   
        }
        matU[n][i]=(matU[n][i]-sum)/matU[i][i];
    }
    return;
}

//Function that changes two given rows of a matrix
void canviarFiles (double **mat, int fila1, int fila2) {
    double *aux;
    aux = mat[fila1];
    mat[fila1] = mat[fila2];
    mat[fila2] = aux;
    return;
}
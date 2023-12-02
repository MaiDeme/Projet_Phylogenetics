/*--------------------------------
           Master BIBS
     Universite Paris-Saclay
   Projet MiniInfo 1 2023-2024

Sujet propose par George Marchment
      george.marchment@lisn.fr
----------------------------------*/

#include "utils.h"

const int GAP = -5;
const int similarity_matrix[NB_BASE][NB_BASE] = {{10, -1, -3, -4, GAP},
                                                 {-1, 7, -5, -3, GAP},
                                                 {-3, -5, 9, 0, GAP},
                                                 {-4, -3, 0, 8, GAP},
                                                 {GAP, GAP, GAP, GAP, 0}};

/*----------------------------------------------------------------------------------
Fonctions permettant de calculer le score d'alignement entre 2 chaînes de caracteres
------------------------------------------------------------------------------------*/

/*
Input : Un caractere
Output : Un entier
Main : Fonction qui prend un caractere et retourne l'index correspondant (voir l'enum Base_Azotee)
*/
int get_val_base(char a)
{
    int index;
    switch (a)
    {
    case 'A':
        index = 0;
        break;
    case 'C':
        index = 1;
        break;
    case 'G':
        index = 2;
        break;
    case 'T':
        index = 3;
        break;
    case '-':
        index = 4;
        break;
    }
    return index;
    // TODO
}

/*
Input : Deux caracteres
Output : Un entier
Main : Fonction qui prend en entier 2 caracteres et qui retourne
       le score entre les 2 caracteres suivant la matrice de similarite
*/
int similarity_score(char ch1, char ch2)
{
    int i1 = get_val_base(ch1);
    int i2 = get_val_base(ch2);
    return similarity_matrix[i1][i2];
    // TODO
}

/*
Input : Deux chaînes de caracteres s
Output : Un entier
Main : Fonction qui prend en entier 2 chaînes de caracteres et qui retourne
       le score d'alignement entre les deux chaînes
*/
int score_alignement(char *alignement1, char *alignement2)
{
    int score = 0;
    for (int i = 0; alignement1[i] != '\0'; i++)
    {
        score += similarity_score(alignement1[i], alignement2[i]);
    }
    return score;
    // TODO
}

/*
Input : Deux chaînes de caracteres et un score
Output : None
Main : Procedure qui prend en entier 2 chaînes de caracteres et un score
       Et qui fait un bel affichage montrant l'alignement et le score
*/
void print_quality_alignement(char *ali1, char *ali2, int score)
{
    printf("Le score d'alignement : %d\n", score);
    printf("\t%s\n", ali1);
    printf("\t%s\n", ali2);
    // TODO
}

/*----------------------------
Algorithme de Needleman-Wunsch
------------------------------*/

/*
Input : 2 entiers et une matrice
Output : None
Main : Procedure qui Initialise la matrice M
*/
void initialise_M(int n, int m, int M[][m])
{
    for (int i = 0; i <n; i++)
    {
        M[i][0] =GAP*i;
    }

    for (int j = 0; j <m; j++)
    {
        M[0][j] =GAP*j ;
    }

    // TODO
}

/*
Input : 2 entiers et une matrice
Output : None
Main : Procedure qui Initialise la matrice T
*/
void initialise_T(int n, int m, char T[][m])
{
    for (int i = 0; i <n; i++)
    {
        T[i][0] = 'u';
    }

    for (int j = 0; j <m; j++)
    {
        T[0][j] = 'l';
    }

    T[0][0] = 'o';
}

/*
Input : un entier
Output : un caractere
Main : Fonction qui retourne le caractere correspondant au maximum index de la formule Mij
*/
char symbole(int entier)
{
    if (entier == 0)
    {
        return 'd';
    }
    else if (entier == 1)
    {
        return 'l';
    }
    return 'u';
}

/*
Input : une chaîne de caracteres
Output : None
Main : Procedure qui inverse une chaîne de caracteres
*/
void reverse_string(char *str)
{
    for (int i = 0; i < (strlen(str)-1)/2; i++)
    {
        int j = strlen(str) - 1 - i;
        char _ = str[j];
        str[j] = str[i];
        str[i] = _;
    }
}

/*
Input : - 2 Sequences
        - 2 index
        - 2 tailles
        - une matrice de taille n*m
        - Deux pointeurs sur des entiers
Output : None
Main : Procedure qui applique la formule Mij et qui sauvegarde
       la plus grande valeur dans le premier pointeur et l'index
       correspondant (0, 1 ou 2) dans le deuxieme.
*/
void fonction_Mij(Sequence *s1, Sequence *s2, int i, int j, int n, int m, int M[][m], int *max, int *index)
{
    int diag = M[i - 1][j - 1] + similarity_score(s1->seq[i-1], s2->seq[j-1]);
    int left = M[i][j - 1] + similarity_score(s1->seq[i-1], '-');
    int up = M[i - 1][j] + similarity_score(s2->seq[j-1], '-');
    *max = diag;
    *index = 0;

    if (*max <= left)
    {
        *max = left;
        *index = 1;
    }
    if (*max <= up)
    {
        *max = up;
        *index = 2;
    }
}

/*
Input : - 2 Sequences
        - Deux chaînes de caracteres
Output : None
Main : Procedure qui applique l'algorithme Needleman-Wunsch
       sur les 2 sequences et qui sauvegarde les alignements
       dans les 2 pointeurs
*/
void needleman_wunsch(Sequence seq1, Sequence seq2, char *alignement1, char *alignement2)
{
    // on declare et initialise les matrices M et T
    int n = strlen(seq1.seq);
    int m = strlen(seq2.seq);
    int M[n+1][m+1];
    initialise_M(n+1, m+1, M);
    char T[n+1][m+1];
    initialise_T(n+1, m+1, T);
    int max=0;
    int index=0;
    for (int i = 1; i < n+1; i++)
    {
        for (int j = 1; j < m+1; j++)
        {

            fonction_Mij(&seq1, &seq2, i, j, n+1, m+1, M, &max, &index);
            M[i][j] = max;
            T[i][j] = symbole(index);
        }
    }
    int score = M[n][m];
    int i = n;
    int j = m;
    while (T[i][j] != 'o')
    {
        if (T[i][j] == 'd')
        {
            appendString(alignement1,seq1.seq[i-1]);
            appendString(alignement2,seq2.seq[j-1]);

            i--;
            j--;
        }
        else
        {
            if (T[i][j] == 'l')
            {
                appendString(alignement1, '-');
                appendString(alignement2,seq2.seq[j-1]);
                j--;
            }
            else
            {
                appendString(alignement1, seq1.seq[i-1]);
                appendString(alignement2,'-');
                i--;
            }
        }

    }

    reverse_string(alignement1);
    reverse_string(alignement2);
      // TODO
}
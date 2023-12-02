/*--------------------------------
           Master BIBS
     Universite Paris-Saclay
   Projet MiniInfo 1 2023-2024

Sujet propose par George Marchment
      george.marchment@lisn.fr
----------------------------------*/

#include "utils.h"
#include <math.h>

/*
Input : Deux Sequences
Output : Float
Main : Fonction qui retourne la distance entre deux sequences 
*/
float distance(Sequence seq1, Sequence seq2) {
    float nb_substitution=0.;
    float nb_nucleotide=0.;
    for (int i=0; seq1.seq[i]!='\0';i++){
        if (seq1.seq[i]!='-' && seq2.seq[i]!='-'){
            nb_nucleotide++;
            if(seq1.seq[i]!=seq2.seq[i]){
                nb_substitution++;
            }
        }
    }

    return nb_substitution/nb_nucleotide;
}

/*
Input : Float
Output : Float
Main : Fonction qui applique la correction de Jukes-Cantor
*/
float jukes_cantor(float x) {
    return -(3./4)*logf(1-(4./3)*x);
}

/*-------------------------------------------------
Fonctions de manipulation de la matrice de distance
---------------------------------------------------*/

/*
Input : un entier et Une matrice de float
Output : None
Main : Procedure qui initialise une matrice à une matrice nulle
*/
void initialise_matrice(int entries, float matrice_distance[][entries]) {
    for (int i=0;i<entries;i++){
        for (int j=0;j<entries;j++){
            matrice_distance[i][j]=0.;
        }

    }
}

/*
Input : Deux entiers et Une matrice de float
Output : None
Main : Procedure qui print une matrice
*/
void print_matrix_float(int n, int m, float matrix[][m]) {
    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            if (j>=i){
                printf(".        ");
            }else{
                printf("%f ",matrix[i][j]);
            }
        }
        printf("\n");

    }
}

/*
Input : entier, matrice de float et une liste de sequence
Output : None
Main : Procedure qui remplit la matrice avec la distance entre les sequences
*/
void fill_distance_matrix(int entries, float matrice_distance[][entries], Sequence sequences[]) {
    for (int i=0;i<entries;i++){
        for (int j=0;j<entries;j++){
            matrice_distance[i][j]=jukes_cantor(distance(sequences[i],sequences[j]));;
        }

    }
}

/*
Input : Un fichier
Output : None
Main : Fonction qui prend une adresse en entree et qui calcule et affiche la matrice de distance correspondant aux sequences
*/
void show_distance_matrix(char* file_aligne){
    int nb_entries = get_number_entries(file_aligne);
    Sequence tab_sequences_aligne[nb_entries];
    parse_file(file_aligne, tab_sequences_aligne);

    float matrice_distance[nb_entries][nb_entries];
    initialise_matrice(nb_entries, matrice_distance);
    fill_distance_matrix(nb_entries, matrice_distance, tab_sequences_aligne);
    printf("La matrice de distance calculee pour le fichier '%s' :\n", file_aligne);
    print_matrix_float(nb_entries, nb_entries, matrice_distance);
    printf("\n");
}

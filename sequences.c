/*--------------------------------
           Master BIBS
     Universite Paris-Saclay
   Projet MiniInfo 1 2023-2024
 
Sujet propose par George Marchment
      george.marchment@lisn.fr
----------------------------------*/

#include "utils.h"

/*------------------------------------------------
Fonctions de manipulation de chaînes de caracteres
--------------------------------------------------*/

/*
Input : une chaîne de caracteres et un caractere
Output : None
Main : Ajoute le caractere à la fin de la chaîne de caracteres
*/
void appendString(char* str, char c) {
    int length = strlen(str);
    str[length] = c;
    str[length + 1] = '\0';
}

/*
Input : une chaîne de caracteres
Output : None
Main : Initialise la chaîne de caracteres comme une chaîne vide
*/
void set_empty_string(char* str) {
    str[0] = '\0';
}



/*------------------------------------
Fonctions de manipulation de sequences
--------------------------------------*/

/*
Input : une sequence et une chaîne de caracteres
Output : None
Main : Initialise la chaîne de caracteres en tant qu'ID pour la sequence
*/
void set_id(Sequence* sequence, char* ID) {
    int length = strlen(ID);
    for (int i = 0; i < length; i++) {
        sequence->ID[i] = ID[i];
    }
    sequence->ID[length] = '\0';
}

/*
Input : une sequence et une chaîne de caracteres
Output : None
Main : Initialise la chaîne de caracteres en tant que sequence pour l'objet sequence
*/
void set_seq(Sequence* sequence, char* seq) {
    int length = strlen(seq);
    for (int i = 0; i < length; i++) {
        sequence->seq[i] = seq[i];
    }
    sequence->seq[length] = '\0';
}

/*
Input : une sequence et un entier
Output : Un caractere
Main : Retourne le nucleotide correspondant à l'index
*/
char get_nucleotide(Sequence* sequence, int i) {
    return sequence->seq[i];
}

/*
Input : une sequence et deux chaînes de caracteres
Output : None
Main : Initialise l'ID et la sequence pour l'objet sequence
*/
void set_sequence(Sequence* sequence, char* ID, char* seq) {
    set_id(sequence, ID);
    set_seq(sequence, seq);
}

/*
Input : une sequence
Output : None
Main : Effectue l'affichage de la sequence
*/
void affiche_sequence(Sequence* sequence) {
    printf("\t* ID : %s\n",sequence->ID);
    printf("\t* Sequence : %s\n",sequence->seq);
    //DONE
}


/*-------------------------------
Fonctions de manipulation de fichiers
---------------------------------*/

/*
Input : Adresse d'un fichier
Output : Nombre de sequences dans le fichier
Main : Compte le nombre de sequences presentes dans le fichier
*/
//https://www.geeksforgeeks.org/c-program-to-read-contents-of-whole-file/
int get_number_entries(char* address) {
    int num = 0;
    FILE* ptr;
    char ch;
    ptr = fopen(address, "r");
    if (NULL == ptr) {
        printf("Le fichier ne peut pas être ouvert\n");
    }
    do {
        ch = fgetc(ptr);
        if (ch == '>') {
            num += 1;
        }
    } while (ch != EOF);
    fclose(ptr);
    return num;
}

/*
Input : Adresse d'un fichier
Output : Une chaîne de caracteres
Main : Fonction qui lit un fichier et retourne la chaîne de caracteres
*/
//https://stackoverflow.com/a/15713419
char* readFile(char* fileName) {
    FILE* file = fopen(fileName, "r");
    char* code;
    size_t n = 0;
    int c;

    if (file == NULL) return NULL; // Impossible d'ouvrir le fichier
    fseek(file, 0, SEEK_END);
    long f_size = ftell(file);
    fseek(file, 0, SEEK_SET);
    code = malloc(f_size);

    while ((c = fgetc(file)) != EOF) {
        code[n++] = (char)c;
    }
    code[n] = '\0';
    return code;
}

/*
Input : - une chaîne de caracteres correspondant au code
        - un entier correspondant à la position start, à partir de laquelle on va commencer à chercher
        - Une sequence qu'on va set à partir de l'information trouvee dans le code
Output : - un entier correspondant à la position où on a fini de "lire"
Main : Fonction qui prend le code ainsi qu'une position start, elle va chercher la prochaine sequence 
       definie apres le start et set l'objet sequence avec cette information.
       Elle retourne ensuite la derniere position qui est lu. Pour le cas de la derniere
       sequence elle retourne -1 
*/
int extract_next_sequence(char* code, int start, Sequence* sequence) {
    char seq[seq_MAX_LENGTH];
    set_empty_string(seq);
    char ID[ID_MAX_LENGTH];
    set_empty_string(ID);
    int i=start;

    //on trouve le > du fasta
    while (code[i]!='>'){
        i++;
    }
    i++;
    //tant qu'on est sur la ligne du > on avance dans les i
    while(code[i]!='\n'){
        appendString(ID,code[i]);
        i++;
    }
    i++;

    while (code[i]!='>' && code[i]!='\0'){
        if (code[i]!='\n'){
            appendString(seq,code[i]);

        }
        i++;
    }
    set_sequence(sequence,ID,seq);

    if (code[i]=='\0'){
        return -1;
    }else{
        return i-1;
    }
    //DONE
}

/*
Input : Adresse d'un fichier et une liste de sequences
Output : None
Main : Fonction qui lit un fichier, remplit la liste avec les sequences trouvees dans le fichier
*/
void parse_file(char* address, Sequence tab_sequences[]) {
    char *data=readFile(address);
    int nseq_=get_number_entries(address);
    int start=0;

    for (int i=0;i<nseq_;i++){
        start=extract_next_sequence(data,start,&tab_sequences[i]);
    }
    //TODO

}

/*
Input : Adresse d'un fichier
Output : None
Main : Parse le fichier, sauvegarde les sequences dans une liste. 
       Puis affiche une sequence une par une.
*/
void show_sequences_file(char* file){
    int nb_entries = get_number_entries(file);
    printf("Nombre d'entrees : %d\n", nb_entries);
    printf("\n");
    
    Sequence tab_sequences[nb_entries];
    parse_file(file, tab_sequences);

    for (int i = 0; i < nb_entries; i++) {
        affiche_sequence(&tab_sequences[i]);
        printf("\n");
    }
}

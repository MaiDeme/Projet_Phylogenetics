#include "utils.h"

/*
--------------------------------
           Master BIBS
     Universite Paris-Saclay
   Projet MiniInfo 1 2023-2024

Sujet propose par George Marchment
      george.marchment@lisn.fr
----------------------------------
*/

int max(int a, int b)
{
    if (a > b)
        return a;
    return b;
}

int min(int a, int b)
{
    if (a < b)
        return a;
    return b;
}

/*
---------------------------
Fonctions de matrice de float
-----------------------------*/

/*
Input : un entier, 2 matrices de float
Output : None
Main : Procedure qui copie l'ancienne matrice dans le nouveau
*/
void set_copy(int entries, float new_tab[][entries], float old_tab[][entries])
{
    for (int i = 0; i < entries; i++)
    {
        for (int j = 0; j < entries; j++)
        {
            new_tab[i][j] = old_tab[i][j];
        }
    }
}

/*
Input : - 2 entiers
        - une matrice de float
        - un pointeur sur un float
        - 2 pointeurs sur des entiers
Output : None
Main : Procedure qui trouve la valeur min dans la matrice Inferieur et stocke cette valeur
       dans le pointeur et ces index dans les 2 autres pointeurs
*/
void find_min_index_distance_matrix(int entries, int nb_noeud, float matrice_distance[][entries], float *min, int *i_min, int *j_min)
{
    *min = matrice_distance[1][0];
    *i_min = 1;
    *j_min = 0;
    for (int i = 1; i < nb_noeud; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (*min > matrice_distance[i][j])
            {
                *min = matrice_distance[i][j];
                *i_min = i;
                *j_min = j;
            }
        }
    }
    // TODO
}

/*--------------------------------
Fonctions de manipulation de noeud
----------------------------------*/
/*
Input : Un Noeud
Output : Entier
Main : Fonction qui retourne 1 si un noeud est une feuille, 0 sinon
*/
int est_feuille(Noeud *e)
{
    if (e->suivant_left == NULL && e->suivant_right == NULL)
    {
        return 1;
    }
    return 0;
}

/*
Input : None
Output : pointeur sur un noeud
Main : qui cree un noeud vide
*/
Noeud *new_noeud()
{
    Noeud *n = (Noeud *)malloc(sizeof(Noeud));
    n->suivant_left = NULL;
    n->suivant_right = NULL;
    n->valeur = NULL;
    n->nb_noeud = 1;
    return n;
}

/*
Input : pointeur sur un noeud
Output : pointeur sur un noeud
Main : qui cree une copie du noeud, en copiant chaque element un par un
*/
Noeud *create_copy(Noeud *e)
{
    Noeud *n = new_noeud();
    n->valeur = e->valeur;
    n->nb_noeud = e->nb_noeud;

    if (e->suivant_left != NULL)
    {
        Noeud *left = create_copy(e->suivant_left);
        n->suivant_left = left;
    }

    if (e->suivant_right != NULL)
    {
        Noeud *right = create_copy(e->suivant_right);
        n->suivant_right = right;
    }
    return n;
}

/*-----------------------------------
Fonctions de manipulation d'affichage
-------------------------------------*/

/*
Input : pointeur sur un arbre
Output : None
Main : procedure qui affiche un arbre
*/
void new_affichage(Arbre *a)
{
    // fonction recursive TODO
    Noeud *e = a->tete;

    printf("\n");
}

/*
Input : pointeur sur un noeud
Output : None
Main : procedure qui affiche le noeud (d'une maniere plate)
*/
void afficher_elem_plat(Noeud *e)
{
    if (e == NULL)
        return;

    if (est_feuille(e))
        printf("%s", e->valeur);
    else
    {
        printf("(");
        afficher_elem_plat(e->suivant_left);
        printf(", ");
        afficher_elem_plat(e->suivant_right);
        printf(")");
    }
}

/*
Input : pointeur sur un arbre
Output : None
Main : procedure qui affiche un arbre (d'une maniere plate)
*/
void afficher_arbre_plat(Arbre *a)
{
    Noeud *e = a->tete;
    printf("\n");
    afficher_elem_plat(e);
    printf("\n");
}

/*-----------------------------------------
Fonctions de manipulation de liste de noeud
-------------------------------------------*/

/*
Input : pointeur d'une Liste de Noeud
Output : entier
Main : Fonction qui retourne le nombre d'elements dans la liste
*/
int get_nb_noeuds(List_Noeuds *list)
{
    return list->nb_elements;
}

/*
Input : pointeur d'une Liste de Noeud et un entier
Output : Pointeur sur un element
Main : Fonction qui retourne l'element qui se trouve à l'index i
*/
Element *get_element(List_Noeuds *list, int i)
{
    if (i < get_nb_noeuds(list))
    {
        Element *temp_ele = list->head;
        for (int k = 0; k < i; k++)
        {
            temp_ele = temp_ele->next;
        }
        return temp_ele;
    }
}

/*
Input : pointeur d'une Liste de Noeud et un entier
Output : Pointeur sur un noeud
Main : Fonction qui retourne le noeud stocke dans l'element qui se trouve à l'index i
*/
Noeud *get_noeud_from_list(List_Noeuds *list, int i)
{
    return get_element(list, i)->data;
}

/*
Input : pointeur d'une Liste de Noeud et un pointeur sur un noeud
Output : None
Main : Procedure qui cree un nouvel element, stocke le noeud dans l'element et qui
       ajoute l'element à la fin de la liste, en faisant attention de bien modifier la
       taille de la liste.
*/
void add_Noeud(List_Noeuds *list, Noeud *n)
{
    Element *temp = get_element(list, get_nb_noeuds(list) - 1);
    Element *to_add = (Element *)malloc(sizeof(Element));
    to_add->next = NULL;
    to_add->data = n;
    temp->next = to_add;

    list->nb_elements++;
}

/*
Input : pointeur d'une Liste de Noeud et 2 entiers
Output : pointeur sur une liste
Main : Fonction qui regroupe ensemble les deux noeuds qui se trouvent à l'index i et j
       et place ce nouveau noeud dans un element qui se place au debut d'une
       nouvelle liste. Ainsi les autres elements de la liste originale sont concatenes
       à la nouvelle dans le même ordre. Ainsi taille(nouvelle) = taille(ancienne)-1.
       Elle retourne le pointeur de cette nouvelle liste.
*/
List_Noeuds *group_together(List_Noeuds *list_param, int i, int j)
{
    // le nouveau noeud qui contient les noeuds i et j
    Noeud *n = new_noeud();
    n->suivant_left = create_copy(get_noeud_from_list(list_param, i));
    n->suivant_right = create_copy(get_noeud_from_list(list_param, j));
    n->nb_noeud = 2;


    //copie de la liste d'avant
    List_Noeuds list_copy=*list_param;

    // Nouvelle liste
    
    List_Noeuds new_list={.head=&(Element){.next=NULL,.data=n},.nb_elements=1};

   
    // on concatene avec l'ancienne liste
    int r = 0;
    while (r < list_param->nb_elements)
    {
        if (r != i && r != j)
        {
            add_Noeud(&new_list, create_copy(get_noeud_from_list(list_param, r)));
        }
        r++;
    }
    return &new_list;
}

/*------------------
Fonctions pour UPGMA
--------------------*/

/*
Input : - taille de la matrice (entier)
        - pointeur d'une Liste de Noeud
        - matrice de float
        - 3 entiers
Output : float
Main : Fonction qui calcule la valeur d'une nouvelle cellule dans la matrice à partir de la formule d_ij,k
*/
float calcule_new_cell(int entries, List_Noeuds *list, float matrice_distance[][entries], int i, int j, int k)
{
    // TODO
    int ni = get_noeud_from_list(list, i)->nb_noeud;
    int nj = get_noeud_from_list(list, j)->nb_noeud;

    float s1 = matrice_distance[i][k] * ni / (ni + nj);
    float s2 = matrice_distance[j][k] * nj / (ni + nj);
    return s1 + s2;
}

/*
Input : - taille de la matrice (entier)
        - pointeur d'une Liste de Noeud
        - matrice de float
Output : pointeur d'une Liste de Noeud
Main : Fonction qui effectue une etape de l'algorithme de UPGMA et qui retourne un pointeur
      sur une nouvelle liste de noeuds
*/
List_Noeuds *fuse_matrice_upgma(int entries, List_Noeuds *list, float matrice_distance[][entries])
{
    // TODO
    float min;
    int i_min, j_min;
    int nb_noeuds = get_nb_noeuds(list);
    find_min_index_distance_matrix(entries, nb_noeuds, matrice_distance, &min, &i_min, &j_min);
    nb_noeuds--;

    float new_matrix[nb_noeuds][nb_noeuds];
    set_copy(nb_noeuds, new_matrix, matrice_distance);

    for (int i = 1, s = 1; i < nb_noeuds; i++)
    {
        for (int j = 0, t = 0; j < i; j++)
        {
            if (j == 0)
            {
                //on remplie la première
                new_matrix[i][0] = calcule_new_cell(entries, list, matrice_distance, i_min, j_min, j);
            }
            else
            {
                if (i < i_min && j < j_min)
                {
                    // tant qu'on a pas attend les lignes des se qui sont mainetnna groupées
                    //  on peut recopier sans s'occuper des indices
                    new_matrix[i][j] = matrice_distance[i - 1][j - 1];
                }
                else
                {
                    if (i == i_min)
                    {
                        s++;
                    }

                    if (j == j_min)
                    {
                        t++;
                    }
                }
                new_matrix[i][j] = matrice_distance[s][t];
            }
            t++;
        }
        s++;
    }

    // on copie dans l'ancienne matrice pour gerer les pointeurs
    set_copy(nb_noeuds, matrice_distance, new_matrix);

    return group_together(list, i_min, j_min);
}

/*
Input : - taille de la matrice (entier)
        - pointeur d'une Liste de Noeuds
        - matrice de float
Output : Un arbre
Main : Fonction qui effectue l'algorithme de UPGMA et qui retourne un arbre
*/
Arbre UPGMA(int entries, List_Noeuds *list, float matrice_distance[][entries])
{

    int nb_noeuds = get_nb_noeuds(list);
    while (nb_noeuds > 1)
    {
        list = fuse_matrice_upgma(entries, list, matrice_distance);
        nb_noeuds = get_nb_noeuds(list);
    }

    Arbre a;
    a.tete = list->head->data;

    return a;
}

/*-----------------------------
Fonctions pour Neighbor Joining
-------------------------------*/

/*
Input : - taille de la matrice (entier)
        - Nombre de noeuds
        - Liste de float
        - matrice de float
Output : None
Main : Procedure qui calcule les S de la matrice de distance donnee en parametre en suivant les regles definies dans le lien
*/
void calcul_S(int entries, int nb_noeuds, float S[nb_noeuds], float matrice[][entries])
{
    // TODO
    for (int i = 0; i < entries; i++)
    {
        for (int j = 0; j < entries; j++)
        {
            S[i] += matrice[j][i];
        }
        S[i] /= (entries - 2);
    }
}

/*
Input : - taille de la matrice (entier)
        - Nombre de noeuds
        - Liste de float
        - matrice de float
        - pointeur sur un float
        - 2 pointeurs sur des entiers
Output : None
Main : Procedure qui calcule le plus petit Mij et qui stocke les i et j dans les pointeurs
*/
void calcule_pair_Mij(int entries, int nb_noeuds, float S[nb_noeuds], float matrice[][entries], float *min_val, int *i_min, int *j_min)
{
    // TODO
    // on va stocker les Mij dans une nouvelle matrice
    float new_matrix_pair[entries][entries];
    for (int i = 0; i < entries; i++)
    {
        for (int j = 0; j < entries; j++)
        {
            new_matrix_pair[i][j] = matrice[i][j] - S[i] - S[j];
        }
    }

    find_min_index_distance_matrix(entries, nb_noeuds, new_matrix_pair, min_val, i_min, j_min);
}

/*
Input : - taille de la matrice (entier)
        - pointeur d'une Liste de Noeud
        - matrice de float
Output : pointeur d'une Liste de Noeud
Main : Fonction qui effectue une etape de l'algorithme de Neighbor Joining et qui retourne un pointeur
      sur une nouvelle liste de noeuds
*/
List_Noeuds *fuse_matrice_NJ(int entries, List_Noeuds *list, float matrice_distance[][entries])
{
    // TODO
    float min;
    int i_min, j_min;

    find_min_index_distance_matrix(entries, get_nb_noeuds(list), matrice_distance, &min, &i_min, &j_min);
    // TODO trouver un moyen de gérer la matrice de distance

    return group_together(list, i_min, j_min);
}

/*
Input : - taille de la matrice (entier)
        - pointeur d'une Liste de Noeud
        - matrice de float
Output : Un arbre
Main : Fonction qui effectue l'algorithme de Neighbor Joining et qui retourne un arbre
*/
Arbre Neighbor_Joining(int entries, List_Noeuds *list, float matrice_distance[][entries])
{

    int nb_noeuds = get_nb_noeuds(list);
    while (nb_noeuds > 1)
    {
        list = fuse_matrice_NJ(entries, list, matrice_distance);
        nb_noeuds = get_nb_noeuds(list);
    }

    Arbre a;
    a.tete = list->head->data;
    return a;
}

//---------------------------------------------------------------------------------------

/*
Input : - une chaîne de caracteres correspondant à une adresse
        - l'algorithme de construction voulu, un caractere (U ou N)
Output : None
Main : Une procedure qui prend une adresse d'un fichier, le parse, extrait les sequences
       alignees et ensuite construit et affiche l'arbre phylogenetique correspondant
       Soit avec l'algo UPGMA ou Neighbor joining.
*/
void show_tree(char *file_aligne, char Algorithme)
{
    int nb_entries = get_number_entries(file_aligne);
    Sequence tab_sequences_aligne[nb_entries];
    parse_file(file_aligne, tab_sequences_aligne);

    float matrice_distance[nb_entries][nb_entries];
    initialise_matrice(nb_entries, matrice_distance);
    fill_distance_matrix(nb_entries, matrice_distance, tab_sequences_aligne);

    List_Noeuds list;
    list.head = NULL;
    Noeud *n = new_noeud();
    n->valeur = tab_sequences_aligne[0].ID;

    Element *temp = (Element *)malloc(sizeof(Element));
    temp->next = NULL;
    temp->data = n;
    list.head = temp;
    list.nb_elements = 1;

    for (int i = 1; i < nb_entries; i++)
    {
        Noeud *n = new_noeud();
        n->valeur = tab_sequences_aligne[i].ID;
        add_Noeud(&list, n);
    }

    Arbre a;
    if (Algorithme == 'U')
    {
        printf("Arbre UPGMA construit pour le fichier '%s' :\n", file_aligne);
        a = UPGMA(nb_entries, &list, matrice_distance);
    }
    else
    {
        printf("Arbre NEIGHBOR JOINING construit pour le fichier '%s' :\n", file_aligne);
        a = Neighbor_Joining(nb_entries, &list, matrice_distance);
    }
    afficher_arbre_plat(&a);
}

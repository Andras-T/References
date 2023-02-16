#ifndef SOURCE_H
#define SOURCE_H

typedef struct ListaElem
{
    char nehezseg[10]; //nehézségi szint
    char kerdes[120];
    char A[100];
    char B[100];
    char C[100];
    char D[100];
    char M[10]; //megoldás
    struct ListaElem *kov;
} ListaElem;

typedef struct Seged
{
    int szint;    //kiválasztott szint
    int gep;      //gép segítsége
    int skip;     //skip segítség
    int kilep;    //kilépés jelzése
    int pontszam; //elért pontszám
} Seged;

#endif
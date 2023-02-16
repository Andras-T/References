#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "source.h"
#include "listakezelo.h"

bool megfelel(char *buf, int szint)
{
    long szam = strtol(buf, 0, 10); //sor elején lévő szám
    if (szam == szint)
    {
        return true;
    }

    return false;
}

ListaElem *hozzafuz(ListaElem *eleje)
{
    ListaElem *uj;
    uj = (ListaElem *)malloc(sizeof(ListaElem));
    if (uj == NULL)
    {
        printf("Nem sikerült memóriát foglalni");
        return 0;
    }

    uj->kov = NULL;

    if (eleje == NULL)
    {
        // üres listánál ez lesz az első elem
        eleje = uj;
    }
    else
    {
        // ha nem üres a lista, az utolsó után fűzzük
        ListaElem *mozgo = eleje;
        while (mozgo->kov != NULL)
        {
            mozgo = mozgo->kov;
        }
        mozgo->kov = uj;
    }
    return eleje;
}

ListaElem *utolso(ListaElem *eleje)
{
    while (eleje->kov != NULL)
    {
        eleje = eleje->kov;
    }

    return eleje; //utolsó elem
}

void szet_valogat(char *buf, ListaElem *aktualis)
{
    int i = 0;
    char eddig[] = ";";
    char *ptr = strtok(buf, eddig);
    strcpy(aktualis->nehezseg, ptr);

    while (ptr)
    {
        switch (i)
        {
        case 1:
            strcpy(aktualis->kerdes, ptr);
            break;
        case 2:
            strcpy(aktualis->A, ptr);
            break;
        case 3:
            strcpy(aktualis->B, ptr);
            break;
        case 4:
            strcpy(aktualis->C, ptr);
            break;
        case 5:
            strcpy(aktualis->D, ptr);
            break;
        case 6:
            strcpy(aktualis->M, ptr);
            break;
        default:
            break;
        }
        ptr = strtok(NULL, eddig);
        i++;
    }
}

void felszabadit(ListaElem *eleje)
{
    ListaElem *iter;
    iter = eleje;
    while (iter != NULL)
    {
        ListaElem *kov = iter->kov; /* következő elem */
        free(iter);
        iter = kov;
    }
}
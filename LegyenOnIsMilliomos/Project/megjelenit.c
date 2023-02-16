#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "source.h"
#include "megjelenit.h"

void megjelenit(ListaElem *aktualis, Seged *egyeb)
{
    printf("\n%s\n\n   ■ A: %s\n   ■ B: %s\n   ■ C: %s\n   ■ D: %s\n\nKérem válasszon (A,B,C,D)\n\n",
           aktualis->kerdes, aktualis->A, aktualis->B, aktualis->C, aktualis->D);
    printf("Skip lehetőségek száma: %d    [1-es karakter]\nGép segítségének száma: %d    [2-es karakter]\n",
           egyeb->skip, egyeb->gep);
}

void ketto_megjelenit(ListaElem *aktualis, Seged *egyeb)
{
    char M = aktualis->M[0];
    int random = 0;
    printf("\n%s\n\n", aktualis->kerdes);
    random = rand() % 3 + 1;

    switch (M)
    {
    case 'A':
        if (random == 1)
        {
            printf("   ■ A: %s\n   ■ B:%s ", aktualis->A, aktualis->B);
        }
        else if (random == 2)
        {
            printf("   ■ A: %s\n   ■ C:%s ", aktualis->A, aktualis->C);
        }
        else
        {
            printf("   ■ A: %s\n   ■ D:%s ", aktualis->A, aktualis->D);
        }
        break;
    case 'B':
        if (random == 1)
        {
            printf("   ■ A: %s\n   ■ B:%s ", aktualis->A, aktualis->B);
        }
        else if (random == 2)
        {
            printf("   ■ B: %s\n   ■ C:%s ", aktualis->B, aktualis->C);
        }
        else
        {
            printf("   ■ B: %s\n   ■ D:%s ", aktualis->B, aktualis->D);
        }
        break;
    case 'C':
        if (random == 1)
        {
            printf("   ■ A: %s\n   ■ C:%s ", aktualis->A, aktualis->C);
        }
        else if (random == 2)
        {
            printf("   ■ B: %s\n   ■ C:%s ", aktualis->B, aktualis->C);
        }
        else
        {
            printf("   ■ C: %s\n   ■ D:%s ", aktualis->C, aktualis->D);
        }
        break;
    case 'D':
        if (random == 1)
        {
            printf("   ■ A: %s\n   ■ D:%s ", aktualis->A, aktualis->D);
        }
        else if (random == 2)
        {
            printf("   ■ B: %s\n   ■ D:%s ", aktualis->B, aktualis->D);
        }
        else
        {
            printf("   ■ C: %s\n   ■ D:%s ", aktualis->C, aktualis->D);
        }
        break;
    default:
        break;
    }
    printf("\n\nSkip lehetőségek száma: %d    [1-es karakter]\nGép segítségének száma: %d    [2-es karakter]\n\n",
           egyeb->skip, egyeb->gep);
}

void mentes(int pontszam)
{
    FILE *Pont = fopen("pontszam.txt", "r");
    int max = 0;
    fscanf(Pont, "%d", &max);

    if (max < pontszam)
    {
        FILE *Pont2 = fopen("pontszam.txt", "w");
        fprintf(Pont2, "%d", pontszam); //beírja a legmagasabb pontszámot
        fclose(Pont2);
    }
    printf("---------------------------\nEddigi legmagasabb pontszám:    %d\nAktuális pontszám:  %d\n", max, pontszam);
    fclose(Pont);
    return;
}
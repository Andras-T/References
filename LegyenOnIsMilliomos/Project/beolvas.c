#include <stdio.h>
#include <stdbool.h>
#include "source.h"
#include "beolvas.h"
#include "megjelenit.h"

int szint_beolvas(void)
{
    int max = 15, min = 1;
    int beolvasott;
    printf("Kérem adjon meg egy nehézségi szintet. [0-15]\n");

    while (true)
    {
        int sikerult = scanf("%d", &beolvasott);
        if (sikerult != 1) //sikertelen beolvasás
        {
            printf("Amit megadott, nem szám volt!\n");
            scanf("%*[^\n]"); //bemenet kiürítése
        }
        else if (beolvasott > max || beolvasott < min)
        {
            printf("A megadott szám nem volt megfelelő! (min: %d, max: %d)\n", min, max);
        }
        else //sikeres és megfelelő beolvasás
        {
            return beolvasott;
        }
    }
}

bool valasz(ListaElem *aktualis, Seged *egyeb)
{
    bool eredmeny;
    char V;
    printf("\nVálasz : ");
    scanf("%s", &V);

    if (V == '2' && egyeb->gep > 0)
    {
        egyeb->gep--;
        ketto_megjelenit(aktualis, egyeb);
        eredmeny = valasz(aktualis, egyeb);
        if (eredmeny)
            egyeb->pontszam--; //ketszeres pont hozzaadas ellen
    }
    else if (V == '1' && egyeb->skip > 0)
    {
        egyeb->skip--;
        egyeb->pontszam--; // ne kapjon pontot
        eredmeny = true;
    }
    else if (V == 'q')
    {
        egyeb->kilep = 1;
        eredmeny = false;
    }
    else
    {
        eredmeny = ellenoriz(aktualis, V);
    }
    if (eredmeny)
    {
        egyeb->pontszam++;
    }

    return eredmeny;
}

bool ellenoriz(ListaElem *aktualis, char valasz)
{
    char megoldas;
    megoldas = aktualis->M[0];
    if (valasz >= 'a')
    {
        valasz -= 'a' - 'A';
    }
    if (valasz == megoldas)
    {
        return true;
    }
    return false;
}

bool vege(void)
{
    char beolvasott;
    printf("Kilépéshez [q], folytatáshoz [c]\n");

    while (true)
    {
        scanf("%s", &beolvasott);
        if (beolvasott == 'q') //kilépés
        {
            return false;
        }
        else if (beolvasott == 'c') //folytatás
        {
            return true;
        }
        else
        {
            printf("Rossz bemenet!\nKilépéshez [q], folytatáshoz [c]\n");
        }
    }
}
#ifdef _WIN32
#include <windows.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "debugmalloc.h"

#include "source.h"
#include "listakezelo.h"
#include "megjelenit.h"
#include "beolvas.h"

void jatek_betolt(Seged *egyeb);                 // file beolvasása, játék indítása

//////////////////////////////////////////////////////////////////

int main(void)
{

#ifdef _WIN32 //ékezetes karakterek
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
#endif

    srand(time(0));
    bool ujra; //játék többszöri elindításához segéd ciklus feltétel
    printf("  ██╗     ███████╗ ██████╗██╗   ██╗███████╗███╗   ██╗     ██████╗ ███╗   ██╗    ██╗███████╗    ███╗   ███╗██╗██╗     ██╗     ██╗ ██████╗ ███╗   ███╗ ██████╗ ███████╗██╗\n  ██║     ██╔════╝██╔════╝╚██╗ ██╔╝██╔════╝████╗  ██║    ██╔═══██╗████╗  ██║    ██║██╔════╝    ████╗ ████║██║██║     ██║     ██║██╔═══██╗████╗ ████║██╔═══██╗██╔════╝██║\n  ██║     █████╗  ██║  ███╗╚████╔╝ █████╗  ██╔██╗ ██║    ██║   ██║██╔██╗ ██║    ██║███████╗    ██╔████╔██║██║██║     ██║     ██║██║   ██║██╔████╔██║██║   ██║███████╗██║\n  ██║     ██╔══╝  ██║   ██║ ╚██╔╝  ██╔══╝  ██║╚██╗██║    ██║   ██║██║╚██╗██║    ██║╚════██║    ██║╚██╔╝██║██║██║     ██║     ██║██║   ██║██║╚██╔╝██║██║   ██║╚════██║╚═╝\n  ███████╗███████╗╚██████╔╝  ██║   ███████╗██║ ╚████║    ╚██████╔╝██║ ╚████║    ██║███████║    ██║ ╚═╝ ██║██║███████╗███████╗██║╚██████╔╝██║ ╚═╝ ██║╚██████╔╝███████║██╗\n  ╚══════╝╚══════╝ ╚═════╝   ╚═╝   ╚══════╝╚═╝  ╚═══╝     ╚═════╝ ╚═╝  ╚═══╝    ╚═╝╚══════╝    ╚═╝     ╚═╝╚═╝╚══════╝╚══════╝╚═╝ ╚═════╝ ╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝\n\n");
    do
    {
        Seged egyeb = {0, 1, 1, 0, 0}; // Seged eszközök resetelése

        egyeb.szint = szint_beolvas();
        jatek_betolt(&egyeb);
        mentes(egyeb.pontszam);

        if (egyeb.kilep == 1) //játék közbeni kilépés
            ujra = false;
        else
            ujra = vege(); //játék végén kilépés vagy újraindítás

    } while (ujra); //újraindít vagy kilép
    return 0;
}

void jatek_betolt(Seged *egyeb)
{
    ListaElem *eleje = NULL;
    FILE *F = fopen("tablazat.txt", "r"); //ez tartalmazza a kérdéseket és válaszokat
    if (F == 0)
    {
        printf("Nem sikerült megnyitni a fájlt\n");
        return;
    }
    char buf[256];
    bool helyes = true; //amíg nincs rossz válasz
    do
    {
        fgets(buf, 256, F);

        if (megfelel(buf, egyeb->szint)) //szint egyezés esetén betölti a feladatot
        {
            buf[strlen(buf)] = '\0'; //levagja a \n-t
            eleje = hozzafuz(eleje);
            ListaElem *aktualis = utolso(eleje);

            szet_valogat(buf, aktualis);

            megjelenit(aktualis, egyeb);
            helyes = valasz(aktualis, egyeb); //válaszadás -> ellenőrzés
        }
    } while (helyes); //újabb kérdés vagy kilépés/új játék
    printf("\n  ╦  ╦╔═╗╔═╗╔═╗┬  \n  ╚╗╔╝║╣ ║ ╦║╣ │  \n   ╚╝ ╚═╝╚═╝╚═╝o \n");
    fclose(F);
    felszabadit(eleje);
}

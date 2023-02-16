#ifndef LISTAKEZELO_H
#define LISTAKEZELO_H

void szet_valogat(char *buf, ListaElem *adatok); //szetválogat structba
bool megfelel(char *buf, int szint);             //megfelelő nehézségü feldat ellenőrzése
ListaElem *hozzafuz(ListaElem *eleje);             //láncolt listát bővíti
void felszabadit(ListaElem *eleje);              //felszabadítja a dinamikusan foglalt adatokat
ListaElem *utolso(ListaElem *eleje);             //megkeresi a láncolt lista aktuális elemét

#endif
#ifndef BEOLVAS_H
#define BEOLVAS_H

int szint_beolvas(void);                        //szint választása
bool valasz(ListaElem *adatok, Seged *egyeb);   //beolvassa a választ
bool ellenoriz(ListaElem *adatok, char valasz); //helyes-e a válasz?
bool vege(void);                                //új játék vagy kilépés
#endif

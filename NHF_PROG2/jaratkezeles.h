#ifndef JARATKEZELES_H
#define JARATKEZELES_H

#include "jarat.h"
#include "vonatkezeles.h"

namespace Jaratkezeles{

/*!
 *  Kivalasztott vonat visszaadasa
 */
std::string getvonat(vectorClass <Vonat>& vonatok);

/*!
 *  Szetvalasztja a sztringet és visszaad egy "Datum"-ot
 */
Datum date_szetvalaszt(std::string ido);

/*!
 *  ido bevitele
 */
Datum gettime();

/*!
 *  uj jarat hozzaadasa
 */
void addjarat(vectorClass <Vonat>& vonatok, vectorClass <Jarat>& jaratok, Allomas& allomasok);

/*!
 *  megjelenites
 */
void megjelenit(vectorClass <Jarat>& jaratok, vectorClass <Vonat>& vonatok);

/*!
 *  utolso jarat torlese
 */
void deljarat(vectorClass <Jarat>& jaratok);

/*!
 *  a kivalasztott jarathoz foglal
 */
void foglal(Jarat& aktualis, vectorClass <Vonat>& vonatok);

}

#endif // !JARATKEZELES_H

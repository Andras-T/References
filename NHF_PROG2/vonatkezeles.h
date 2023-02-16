#ifndef VONATKEZELES_H
#define VONATKEZELES_H

#include "vector.h"
#include "vonat.h"


/*!
 *  elvalasztosor megjelenitese
 */
void elvalaszto();

namespace Vonatkezeles {

	/*!
	 *  uj vonat hozzaadasa
	 */
	void addvonat(vectorClass <Vonat>& vonatok, const char* filename = "vonatok.txt");

	/*!
	 *  torli az utolso vonatot
	 */
	void delvonat(vectorClass <Vonat>& vonatok, const char* filename = "vonatok.txt");

	/*!
	 *  uj kocsi hozzaadasa vonathoz
	 */
	void kocsi_add(vectorClass <Vonat>& vonatok, const char* filename = "vonatok.txt");

	/*!
	 *  torli a kivalasztott vonat utolso kocsijat
	 */
	void kocsi_del(vectorClass <Vonat>& vonatok, const char* filename = "vonatok.txt");

	/*!
	 *  Szetvalasztja a sztringet és visszaad egy "Vonat"-ot
	 */
	void szet_valogat(std::string buf, vectorClass <Vonat>& vonatok);

	/*!
	 *  vonatok beolvasasa fajlbol
	 */
	void vonat_betolt(vectorClass <Vonat>& vonatok, const char* filename = "vonatok.txt");

	/*!
	 *  vonatok es a hozza tartozo kocsik megjelenitese
	 */
	void megjelenit(vectorClass <Vonat>& vonatok);

}
#endif // !VONATKEZELES_H


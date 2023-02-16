#ifndef JARAT_H
#define JARAT_H

#include "vonat.h"
#include "allomas.h"
#include "datum.h"
#include "jegyek.h"
#include "diak.h"
#include "felnott.h"
#include "bicikli.h"

class Jarat
{
public:

	/*!
	 *  Konstruktor.
	 *
	 *      @param vonat 
	 *      @param indulas idopontja
	 *      @param erkezes idopontja
	 *      @param indulo allomas
	 *      @param erkezo allomas
	 */
	Jarat(const std::string vonat = " ", const Datum ind = 0, const Datum erk = 0, const std::string i_all = "", const std::string e_all = "");
	
	/*!
	 *  Destruktor.
	 */
	~Jarat() {}

	/*!
	 *  Operator= jarat.
	 *
	 *      @param masik jarat
	 *
	 *      @return jarat
	 */
	Jarat& operator=(const Jarat& other);

	/*!
	 *  lefoglal egy helyet a kivalasztott osztalyra megfelelo jegytipusra
	 *
	 *      @param kocsiosztaly
	 *      @param jegytipus 
	 */
	void foglal(const int k_osztaly, const int jegytipus, vectorClass <Vonat>& vonatok);

	/*!
	 *  jarat kiiro fuggveny
	 */
	void kiir(vectorClass <Vonat>* vonatok) const;

	/*!
	 *  jegyek kiirasa
	 */
	void jegykiir() const;

private:
	
	//Jarathoz tartozo vonat
	std::string vonat;

	//indulo/erkezeo allomas
	std::string i_allomas;
	std::string e_allomas;

	//indulasi/erkezesi ido
	Datum indulas;
	Datum erkezes;

	//jegyek: normal
	vectorClass <Jegyek> jegyek;
	// diakjegy
	vectorClass <DiakJegy> diakjegyek;
	//felnottjegy
	vectorClass <FelnottJegy> felnottjegyek;
	//biciklijegy
	vectorClass <Bicikli> biciklijegyek;
};

#endif // !JARAT_H

#ifndef VONAT_H
#define VONAT_H

#include"vector.h"
#include "kocsi.h"

class Vonat
{
public:

	/*!
	 *  Konstruktor.
	 *
	 *      @param nev 
	 */
	Vonat(const std::string nev = " ") : nev(nev), kocsik() {}

	/*!
	 *  Masolo konstruktor.
	 *
	 *      @param nev 
	 */
	Vonat(const char* nev) : nev(nev), kocsik() {}

	/*!
	 *  Destruktor.
	 */
	~Vonat();
	
	/*!
	 *  Operator= vonat.
	 *
	 *      @param masik vonat
	 *
	 *      @return vonat
	 */
	Vonat& operator=(const Vonat& other);
	
	//compare operator
	bool operator==(const Vonat other) const;

	/*!
	 *  indexelo operator
	 *
	 *      @param index
	 *
	 *      @return kocsi
	 */
	Kocsi& operator[](const size_t i) const;

	/*!
	 *  Atnevezes
	 *
	 *      @param nev
	 */
	void atnevez(const std::string name) { this->nev = name; }

	/*!
	 *  uj kocsi felvetele
	 *
	 *      @param max hely   
	 *      @param osztaly   
	 *      @param kocsi nev 
	 */
	void addkocsi(const int maxhely, const int osztaly, const int kocsi_nev);
	
	/*!
	 *  utolso kocsi eltávolitasa
	 */
	void delkocsi() { kocsik.pop(); }

	/*!
	 *  kiiro standard kimenetre (kocsikkal egyutt)
	 */ 
	void kiir() const;

	/*!
	 *  vonatnevet adja meg
	 *
	 *      @return vonatnev
	 */
	std::string getnev() const {return nev;}

	/*!
	 *  kiuriti a vonatot
	 */
	void kiurit();

	/*!
	 *  megadja a kocsik szamat
	 *
	 *      @return kocsik szama
	 */
	size_t getsize() const;

	/*!
	 *  visszaad egy kocsit a megadott osztallyal
	 *
	 *      @param osztaly 
	 *
	 *      @return kocsi
	 */
	Kocsi* getkocsi(const int osztaly) const;

private:	/*-------------------PRIVAT---------------------*/
	
	//vonat neve
	std::string nev;
	
	//Vonathoz tartozo kocsik
	vectorClass <Kocsi> kocsik;
};

/*!
 *  inserter
 */
std::ostream& operator<<(std::ostream& os, const Vonat& train);

#endif // !VONAT_H

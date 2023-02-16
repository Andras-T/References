#ifndef ALLOMAS_H
#define ALLOMAS_H

#include <fstream>
#include <string>
#include "vector.h"

class Allomas
{
private:

	//vasutallomasok
	vectorClass<std::string> allomasok;

public:

	/*!
	 *  Konstruktor.
	 *
	 *      @param melyik fajlbol toltse be
	 */
	Allomas(const char* filename = "allomasok.txt");

	/*!
	 *  Destruktor.
	 */
	~Allomas() {}

	
	/*!
	 *  visszaadja, hogy szerepel-e az allomasok kozott
	 *
	 *      @param all allomasnev
	 *
	 *      @return szerepel-e
	 */
	bool init(const std::string all) const;

	
	/*!
	 *  kibovíti egy ujabb allomassal;
	 *
	 *      @param uj allomasnev       
	 *      @param melyik fajlba bovitse
	 */
	void addAllomas(const std::string uj, const char* filename = "allomasok.txt");

	
	/*!
	 *  utolso allomas eltavolitasa
	 *
	 *      @param filename melyik fajlbol torolje
	 */
	void delAllomas(const char* filename = "allomasok.txt");

	//
	/*!
	 *  @return allomasok szama
	 */
	int givedb() const;	

	/*!
	 *  Operator[]s the allomas.
	 *
	 *      @param index
	 *
	 *      @return  i. allomasnev
	 */
	std::string operator[](const size_t i)
	{
		return allomasok[i];
	}
};

/*!
 *  inserter
 */
std::ostream& operator<<(std::ostream& os, Allomas& a);	//inserter

#endif // !ALLOMAS_H
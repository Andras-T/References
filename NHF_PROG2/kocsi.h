#ifndef KOCSI_H
#define KOCSI_H

#include <iostream>

#include <string> 

class Kocsi
{
public:

	/*!
	 *  Konstruktor.
	 *
	 *      @param maxhely 
	 *      @param osztaly 
	 *      @param nev     
	 */
	Kocsi(int maxhely = 0, int osztaly = 0, int nev = 0);

	/*!
	 *  Destruktor.
	 */
	~Kocsi() {}

	/*!
	 *  Operator= kocsi.
	 *
	 *      @param masik kocsi
	 *
	 *      @return kocsi
	 */
	Kocsi& operator=(const Kocsi& other);

	//
	/*!
	 *  kocsiosztalyt adja vissza
	 *
	 *      @return kocsiosztaly
	 */
	int getosztaly() const { return osztaly; }

	//
	/*!
	 *  kocsiszamot ad vissza
	 *
	 *      @return 
	 */
	int getkocsiszam() const { return nev; }

	/*!
	 *  max ferohelyet adja vissza
	 *
	 *      @return max ferohely
	 */
	int getmax() const { return maxhely; }

	//
	/*!
	 *  szabad helyek szamat adja vissza
	 *
	 *      @return szabad helyek
	 */
	int getszabad() const { return szabad; }

	/*!
	 *  szabad helyek csokkentese
	 */
	void subszabad();

	/*!
	 *  reseteli a kocsit
	 */
	void reset();

	/*!
	 *  a kocsi adatait adja vissza
	 *
	 *      @return kocsi adatai
	 */
	std::string kiiras() const;

private:
	int maxhely;
	int osztaly;
	int nev;
	int szabad;
};

/*!
 *  inserter
 */
std::ostream& operator<<(std::ostream& os, const Kocsi& k);

#endif // !KOCSI_H

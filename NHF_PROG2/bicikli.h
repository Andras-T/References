#ifndef BICIKLI_H
#define BICIKLI_H

#include "jegyek.h"

class Bicikli : public Jegyek
{
public:

	/*!
	 *  Konstruktor.
	 *
	 *      @param hely        
	 *      @param osztaly     
	 *      @param kocsi       
	 *      @param biciklijegy 
	 *      @param alapar      
	 */
	Bicikli(const int hely = 0, const int osztaly = 0, const int kocsi = 0, const double biciklijegy = 420, const double alapar = 500)
		: Jegyek(hely, osztaly, kocsi, alapar), biciklijegy(biciklijegy) {}

	/*!
	 *  Destruktor.
	 */
	~Bicikli() {}


	/*!
	 *  Operator= bicikli.
	 *
	 *      @param masik biciklijegy
	 *
	 *      @return biciklijegy
	 */
	Bicikli& operator=(const Bicikli& other);

	/*!
	 *  Jegyarat adja vissza
	 *
	 *      @return jegyar
	 */
	double getar() const { return biciklijegy; }

private:

	double biciklijegy;
};

/*!
 *  inserter
 */
std::ostream& operator<<(std::ostream& os, const Bicikli& a);

#endif // !BICIKLI_H

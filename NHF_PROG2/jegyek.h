#ifndef JEGYEK_H
#define JEGYEK_H

#include <iostream>

//ososztaly a tobbi jegytipusnak
class Jegyek
{
public:

	/*!
	 *  Konstruktor.
	 *
	 *      @param hely      
	 *      @param osztaly   
	 *      @param kocsiszam 
	 *      @param alapar    
	 */
	Jegyek(const int hely = 0, const int osztaly = 0, const int kocsiszam = 0, const double alapar = 500)
		:hely(hely), osztaly(osztaly), kocsiszam(kocsiszam), alapar(alapar) {}

	/*!
	 *  Destruktor.
	 */
	~Jegyek() {}
	
	/*!
	 *  Operator= jegyek.
	 *
	 *      @param masik jegy
	 *
	 *      @return 
	 */
	Jegyek& operator=(const Jegyek& other);

	/*!
	 *  alaparat adja vissza
	 *
	 *      @return alapar
	 */ 
	double getalap() const { return alapar; }

	/*!
	 *  jegyarat adja vissza
	 *
	 *      @return jegyar
	 */
	virtual double getar() const { return alapar; }

	/*!
	 *  Kocsi osztályat adja meg
	 *
	 *      @return osztaly
	 */
	int getosztaly() const { return osztaly; }

	/*!
	 *  Kocsiszamat adja meg
	 *
	 *      @return kocsiszam
	 */
	int getkocsiszam() const { return osztaly; }

	/*!
	 *  Lefoglalt helyet adja meg
	 *
	 *      @return hely
	 */
	int gethely() const { return hely; }

private:

	int hely;
	int osztaly;
	int kocsiszam;
	double alapar;
};

/*!
 *  inserter
 */
std::ostream& operator<<(std::ostream& os, const Jegyek& a);

#endif // !JEGYEK_H

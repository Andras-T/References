#ifndef FELNOTT_H
#define FELNOTT_H

//#include <iostream>
#include "jegyek.h"

class FelnottJegy : public Jegyek
{
public:

	/*!
	 *  Konstruktor.
	 *
	 *      @param hely    
	 *      @param osztaly 
	 *      @param kocsi   
	 *      @param felar   
	 *      @param alapar  
	 */
	FelnottJegy(const int hely = 0, const int osztaly = 0, const int kocsi = 0, const double felar = 50, const double alapar = 500)
		:Jegyek(hely, osztaly, kocsi, alapar), felar(felar) {}

	/*!
	 *  Destruktor.
	 */
	~FelnottJegy() {}

	/*!
	 *  Operator= felnottjegy.
	 *
	 *      @param masik felnottjegy
	 *
	 *      @return felnottjegy
	 */
	FelnottJegy& operator=(const FelnottJegy& other);

	/*!
	 *  jegyarat adja vissza
	 *
	 *      @return jegyar
	 */ 
	double getar() const { return felar + getalap(); }

private:

	double felar;
};

/*!
 *  inserter
 */
std::ostream& operator<<(std::ostream& os, const FelnottJegy& a);

#endif // !FELNOTT_H

#ifndef DIAK_H
#define DIAK_H

#include "jegyek.h"

class DiakJegy: public Jegyek
{
public:

	/*!
	 *  KonstruKtor.
	 *
	 *      @param hely       
	 *      @param osztaly    
	 *      @param kocsi      
	 *      @param kedvezmeny 
	 *      @param alapar     
	 */
	DiakJegy(const int hely = 0, const int osztaly = 0, const int kocsi = 0, const double kedvezmeny = 50, const double alapar = 500)
		: Jegyek(hely, osztaly, kocsi, alapar), kedvezmeny(kedvezmeny) {}
	
	/*!
	 *  Destruktor
	 */
	~DiakJegy() {}


	/*!
	 *  Operator= diakjegy.
	 *
	 *      @param masik diakjegy
	 *
	 *      @return diakjegy
	 */
	DiakJegy& operator=(const DiakJegy& other);

	/*!
	 *  jegyarat adja vissza
	 *
	 *      @return jegyar
	 */
	double getar() const;

private:

	//kedvezmeny %-ban
	double kedvezmeny;	
};

/*!
 *  inserter
 */
std::ostream& operator<<(std::ostream& os, const DiakJegy& a);

#endif // !DIAK_H

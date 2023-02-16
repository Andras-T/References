#ifndef DATUM_H
#define DATUM_H

#include <iostream>

class Datum
{
private:
	int ora;
	int perc;

public:

	//default ctor
	/*!
	 *  Konstruktor.
	 *
	 *      @param ora
	 *      @param perc
	 */
	Datum(const int h = 0, const int min = 0) : ora(h), perc(min) {}

	/*!
	 *  Masolo konstruKtor.
	 *
	 *      @param masik datum 
	 */
	Datum(const Datum& other) : ora(other.ora), perc(other.perc) {}

	/*!
	 *  Destruktor.
	 */
	~Datum() {}

	/*!
	 *  Operator= datum.
	 *
	 *      @param masik datum
	 *
	 *      @return datum
	 */
	Datum& operator=(const Datum& other)
	{
		ora = other.ora;
		perc = other.perc;
		return *this;
	}

	/*!
	 *  Visszaadja az orat
	 *
	 *      @return ora
	 */
	int get_h() const { return ora; }

	/*!
	 *  visszaadja a percet
	 *
	 *      @return perc
	 */
	int get_min() const { return perc; }
};

/*!
 *  inserter 
 */
std::ostream& operator<<(std::ostream& os, const Datum& a);

#endif // !DATUM_H

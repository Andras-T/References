#include <iostream>
#include "jegyek.h"
#include "memtrace.h"


Jegyek& Jegyek::operator=(const Jegyek& other)
{
	hely = other.hely;
	osztaly = other.osztaly;
	kocsiszam = other.kocsiszam;
	alapar = other.alapar;
	return *this;
}

std::ostream& operator<<(std::ostream& os, const Jegyek& a)
{
	return os << " Kocsisztaly: " << a.getosztaly() << "\n Kocsiszam: " << a.getkocsiszam() << "\n Hely: " << a.gethely() << std::endl;
}
#include "diak.h"
#include "memtrace.h"


DiakJegy& DiakJegy::operator=(const DiakJegy& other)
{
	Jegyek::operator=(other);
	kedvezmeny = other.kedvezmeny;
	return *this;
}


double DiakJegy::getar() const
{
	double sum = getalap() * (1 - kedvezmeny / 100);
	return sum;
}

std::ostream& operator<<(std::ostream& os, const DiakJegy& a)
{
	return os << Jegyek(a) << " Jegyar: " << a.getar() << std::endl;
}

#include "bicikli.h"
#include "memtrace.h"


Bicikli& Bicikli::operator=(const Bicikli& other)
{
	Jegyek::operator=(other);
	biciklijegy = other.biciklijegy;
	return *this;
}

std::ostream& operator<<(std::ostream& os, const Bicikli& a)
{
	return os << Jegyek(a) << " Jegyar: " << a.getar() << std::endl;
}
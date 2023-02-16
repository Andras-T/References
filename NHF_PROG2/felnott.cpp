#include "felnott.h"
#include "memtrace.h"


FelnottJegy& FelnottJegy::operator=(const FelnottJegy& other)
{
	Jegyek::operator=(other);
	felar = other.felar;
	return *this;
}

std::ostream& operator<<(std::ostream& os, const FelnottJegy& a)
{
	return os << Jegyek(a) << " Jegyar: " << a.getar() << std::endl;
}
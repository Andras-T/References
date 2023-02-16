#include "kocsi.h"
#include "memtrace.h"


Kocsi::Kocsi(int maxhely, int osztaly, int nev) :maxhely(maxhely), osztaly(osztaly), nev(nev), szabad(maxhely) {}

Kocsi& Kocsi::operator=(const Kocsi& other)
{
	maxhely = other.maxhely;
	osztaly = other.osztaly;
	nev = other.nev;
	szabad = other.szabad;
	return *this;
}

void Kocsi::subszabad()
{
	if (szabad < 1)
	{
		throw std::invalid_argument("Nincs tobb szabad hely");
	}
	else
	{
		szabad--;
	}
}

void Kocsi::reset()
{
	szabad = maxhely;
}

std::string Kocsi::kiiras() const
{
	std::string szoveg;
	szoveg = "\tKocsiszam: ";
	szoveg += std::to_string(nev);
	szoveg += ", Osztaly: ";
	szoveg += std::to_string(osztaly);
	szoveg += ", Maxhelyek szama: ";
	szoveg += std::to_string(maxhely);
	return szoveg;
}

std::ostream& operator<<(std::ostream& os, const Kocsi& k)
{
	return os << k.kiiras() << std::endl;
}

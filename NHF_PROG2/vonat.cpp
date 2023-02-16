#include "vonat.h"
#include "memtrace.h"


Vonat::~Vonat() {}

Vonat& Vonat::operator=(const Vonat& other)
{
	nev = other.nev;
	/*kocsik.torol();
	if (other.kocsik.size() != 0)
	{
		for (size_t i = 0; i < other.kocsik.size(); i++)
		{
			kocsik.push(other.kocsik[i]);
		}
	}*/
	kocsik = other.kocsik;
	return *this;
}

bool Vonat::operator==(const Vonat other) const
{
	return (this->nev == other.nev);
}

Kocsi& Vonat::operator[](const size_t i) const
{
	return kocsik[i];
}

void Vonat::addkocsi(const int maxhely, const int osztaly, const int kocsi_nev)
{
	Kocsi uj(maxhely, osztaly, kocsi_nev);
	kocsik.push(uj);
}

void Vonat::kiir() const
{
	std::cout << "\n Vonatnev: " << nev << std::endl;
	if (kocsik.size() != 0)

		std::cout << " Kocsik:\n" << kocsik;
}

Kocsi* Vonat::getkocsi(const int osztaly) const
{
	for (size_t i = 0; i < kocsik.size(); i++)
	{
		if (osztaly == kocsik[i].getosztaly() && kocsik[i].getszabad() > 0)
		{
			return &(kocsik[i]);
		}
	}
	throw std::invalid_argument("Nincs ilyen osztalyu kocsi vagy elfogytak a helyek");
}

void Vonat::kiurit()
{
	for (size_t i = 0; i < kocsik.size(); i++)
	{
		kocsik[i].reset();
	}
}

size_t Vonat::getsize() const
{
	return kocsik.size();
}

std::ostream& operator<<(std::ostream& os, const Vonat& train)
{
	os << train.getnev() << ": ";
	if (train.getsize() > 0)
	{
		for (size_t i = 0; i < train.getsize(); i++)
		{
			os << train[i].getkocsiszam() << ", " << train[i].getosztaly() << ", " << train[i].getmax() << "; ";
		}
	}
	os << "\n";
	return os;
}

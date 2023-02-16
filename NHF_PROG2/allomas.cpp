#include "allomas.h"
#include "memtrace.h"


Allomas::Allomas(const char* filename)
{
	// FILE megnyitasa
	std::ifstream MyFile(filename);	
	//ideiglenes tarolo
	std::string buf;	
	while (std::getline(MyFile, buf))	//sor olvasasa
		allomasok.push(buf);	//allomas hozzaadasa
	//FILE bezarasa
	MyFile.close();	
}

bool Allomas::init(const std::string all) const	
{
	return allomasok.init(all);
}

void Allomas::addAllomas(const std::string uj, const char* filename)	
{
	if (allomasok.init(uj))
		throw std::invalid_argument("Mar szerepel");
	allomasok.push(uj);
	std::ofstream MyFile(filename);
	for (size_t i = 0; i < allomasok.size(); i++)
	{
		MyFile << allomasok[i] << "\n";
	}
	MyFile.close();
}

void Allomas::delAllomas(const char* filename)
{
	if (allomasok.size() > 0)
	{
		allomasok.pop();
		std::ofstream MyFile(filename);
		for (size_t i = 0; i < allomasok.size(); i++)
		{
			MyFile << allomasok[i] << "\n";
		}
		MyFile.close();
	}
	else
		throw std::invalid_argument("Nincs tobb elem");
}

int Allomas::givedb() const	
{
	return allomasok.size();
}

std::ostream& operator<<(std::ostream& os, Allomas& a)	//inserter
{
	size_t hossz = a.givedb();
	for (size_t i = 0; i < hossz; i++)
	{
		os << " [" << i + 1 << "] " << a[i] << "\n";
	}
	return os;
}
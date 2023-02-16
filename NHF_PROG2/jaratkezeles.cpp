#include "jaratkezeles.h"
#include "memtrace.h"
#include "vonatkezeles.h"

std::string Jaratkezeles::getvonat(vectorClass <Vonat>& vonatok)
{
	//vonatok megjelenitese
	Vonatkezeles::megjelenit(vonatok);
	std::string nev;
	std::cout << " Valasszon vonatot: ";
	std::cin >> nev;
	//Vonat megkeresese
	for (size_t i = 0; i < vonatok.size(); i++)
	{
		if (vonatok[i] == nev)
		{
			if (vonatok[i].getsize() == 0)
			{
				throw std::invalid_argument("Nincs kocsi ehez a vonathoz");
			}
			return vonatok[i].getnev();
		}
	}
	throw std::invalid_argument("Nincs ilyen vonat");
}

Datum Jaratkezeles::date_szetvalaszt(std::string ido)
{
	//ideiglenes tarolo, elvalaszto
	std::string buffer, delimiter = ":";
	int hour = 0, min = 0;
	//szo vegenek a pozicioja
	size_t pos = 0;
	size_t index = 1;
	while ((pos = ido.find(delimiter)) != std::string::npos || index < 3)
	{
		switch (index)
		{
		case 1:	//ora
			buffer = ido.substr(0, pos);
			hour = atoi(buffer.c_str());
			break;
		case 2:	//perc
			buffer = ido.substr(0, pos);
			min = atoi(buffer.c_str());
			break;
		default:
			break;
		}
		//sztring roviditese
		ido.erase(0, pos + delimiter.length());
		index++;
	}
	if (hour > 23 || hour < 0 || min > 59 || min < 0)
	{
		throw " !!Rossz formatum!!\n ";
	}
	Datum time(hour, min);
	return time;
}

Datum Jaratkezeles::gettime()
{
	std::string ido;
	//elfogadhato formatum pl: 13:34
	std::cout << " Formatum: Ora:perc\n ";
	std::cin >> ido;
	Datum time;
	try
	{
		time = date_szetvalaszt(ido);
	}
	catch (const char*)
	{
		throw std::invalid_argument(" !!Rossz formatum!!\n");
	}

	return time;
}

void Jaratkezeles::addjarat(vectorClass <Vonat>& vonatok, vectorClass <Jarat>& jaratok, Allomas& allomasok)
{
	//vonat hozzaadasa jarathoz
	try
	{
		std::string uj = getvonat(vonatok);

		//indulo es erkezo allomasok hozzaadasa
		std::string indulo, erkezo;
		std::cout << "\n Allomasok:\n" << allomasok << "\n Indulo allomas: ";
		std::cin >> indulo;
		std::cout << " Erkezeo allomas: ";
		std::cin >> erkezo;
		if (!allomasok.init(indulo) || !allomasok.init(erkezo))
		{
			std::cout << " Nincs ilyen allomas\n";
			return;
		}

		//indulasi/erkezesi ido
		Datum ind, erk;
		std::cout << " Indulas:\n";
		ind = gettime();
		std::cout << " Erkezes:\n";
		erk = gettime();

		Jarat jarat(uj, ind, erk, indulo, erkezo);
		std::cout << "\n";
		//uj jarat hozzadasa a tarolohoz
		jaratok.push(jarat);
		jaratok[jaratok.size() - 1].kiir(&vonatok);
	}
	catch (std::invalid_argument kivetel)
	{
		std::cout << " !" << kivetel.what() << "!\n";
		return;
	}
}

void Jaratkezeles::megjelenit(vectorClass <Jarat>& jaratok, vectorClass <Vonat>& vonatok)
{
	std::cout << "\n Jaratok:\n";
	if (jaratok.size() == 0)	//Ha ures
	{
		std::cout << " -Meg nincs jarat-\n";
	}
	for (size_t i = 0; i < jaratok.size(); i++)
	{
		std::cout << "\n [" << i + 1 << "] ";
		jaratok[jaratok.size()-1-i].kiir(&vonatok);
	}
}

void Jaratkezeles::deljarat(vectorClass <Jarat>& jaratok)
{
	try
	{
		jaratok.pop();
	}
	catch (std::invalid_argument kivetel)	//Ha nincs tobb torolheto
	{
		std::cout << " !!" << kivetel.what() << "!!\n";
	}
}

void Jaratkezeles::foglal(Jarat& aktualis, vectorClass <Vonat>& vonatok)
{
	std::string szoveg = " Jegytipus: [0 = alapjegy] [1 = felnottjegy] [2 = diakjegy] [3 = biciklijegy] [4 = Megse]\n ";
	std::cout << szoveg;

	bool sikertelen = true;
	std::string jegytipus, osztaly;
	while (sikertelen == true)
	{
		std::cin >> jegytipus;
		if (jegytipus == "0" || jegytipus == "1" || jegytipus == "2" || jegytipus == "3")
		{
			std::cout << " Osztaly [1/2]\n ";
			std::cin >> osztaly;
			if (atoi(osztaly.c_str()) == 1 || atoi(osztaly.c_str()) == 2)
			{
				sikertelen = false;
				try
				{
					aktualis.foglal(atoi(osztaly.c_str()), atoi(jegytipus.c_str()), vonatok);
				}
				catch (std::invalid_argument kivetel)
				{
					std::cout << kivetel.what();
				}
			}
			else std::cout << "\n Helytelen formatum\n" << szoveg;
		}
		else if (jegytipus == "4") return;
		else std::cout << "\n Helytelen formatum\n" << szoveg;
	}
}
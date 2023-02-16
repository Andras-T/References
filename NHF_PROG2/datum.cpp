#include "datum.h"
#include "memtrace.h"


std::ostream& operator<<(std::ostream& os, const Datum& a)
{

	if (a.get_h() < 10)
	{
		os << "0" << a.get_h() << ":";
	}
	else
	{
		os << a.get_h() << ":";
	}
	if (a.get_min() < 10)
	{
		os << "0" << a.get_min();
	}
	else
	{
		os << a.get_min();
	}
	return os;
}
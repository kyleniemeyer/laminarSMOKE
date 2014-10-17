#if OPENSMOKE_USE_ISAT == 1
if(isatCheck == true)
{
	#include "chemistry_ISAT.H"
}
else
{
	#include "chemistry_DI.H"
}
#else
	#include "chemistry_DI.H"
#endif

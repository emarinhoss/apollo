
#include <wxcreator.h>

template <>
std::string WxCreator_AppendSuffix<float>( const std::string & baseString )
{
	return baseString + "_f";
}

template <>
std::string WxCreator_AppendSuffix<double>( const std::string & baseString )
{
	return baseString + "_d";
}

template <>
std::string WxCreator_AppendSuffix<int>( const std::string & baseString )
{
	return baseString + "_i";
}

std::string WxCreator_AppendRealSuffix( const std::string & baseString, const std::string & realType )
{
	if (realType=="double")
		WxCreator_AppendSuffix<double>( baseString );
	else if (realType=="float")
		return WxCreator_AppendSuffix<float>(baseString);
	else if (realType=="int")
		return WxCreator_AppendSuffix<int>(baseString);
	// unmatched do nothing
	return baseString;
}

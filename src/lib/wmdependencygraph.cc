#include <wmdependencygraph.h>

const char * SubregionCompletionType::edgeTypeEnumText[] = { "all", "ghost", "boundary" };
unsigned SubregionCompletionType::numEdgeTypes = sizeof(edgeTypeEnumText) / sizeof(char *);
#pragma once

#ifdef LIBRARY_EXPORTS
#    define LIBRARY_API __declspec(dllexport)
#else
#    define LIBRARY_API __declspec(dllimport)
#endif

extern "C"
{
	LIBRARY_API void __cdecl Initialize();
	LIBRARY_API bool __cdecl GetNextFrame(
		int* textureFormat, /* rs2_format enum */
		int* textureWidth,
		int* textureHeight,
		const void** texturePixels,
		int* pointsSize,
		float** pointsX,
		float** pointsY,
		float** pointsZ,
		float** textureU,
		float** textureV);
}

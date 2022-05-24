#pragma once

#include <librealsense2/rs.hpp> // Include RealSense Cross Platform API

class RSPC
{
public:
	static RSPC& Instance();

private:
	RSPC(): _pointsX(NULL),
		_pointsY(NULL),
		_pointsZ(NULL),
		_textureU(NULL),
		_textureV(NULL)
	{
	
	}                    
	// Declare RealSense pipeline, encapsulating the actual device and sensors
	rs2::pipeline pipe;

	// Declare pointcloud object, for calculating pointclouds and texture mappings
	rs2::pointcloud pc;

	// We want the points object to be persistent so we can display the last cloud when a frame drops
	rs2::points points;

	void DeallocateArray(float* & points)
	{
		if (points != NULL)
		{
			delete[] points;
			points = NULL;
		}
	}

	void Deallocate()
	{
		DeallocateArray(_pointsX);
		DeallocateArray(_pointsY);
		DeallocateArray(_pointsZ);
		DeallocateArray(_textureU);
		DeallocateArray(_textureV);
	}

public:
	RSPC(RSPC const&) = delete;
	void operator=(RSPC const&) = delete;

	float* _pointsX;
	float* _pointsY;
	float* _pointsZ;
	float* _textureU;
	float* _textureV;

	void Initialize();
	bool GetNextFrame(
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
};

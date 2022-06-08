#include "RSPCSensor.h"
#include "rs_frame.hpp"
#include "RSPCFrame.h"

RSPCSensor& RSPCSensor::Instance()
{
	static RSPCSensor    instance; // Guaranteed to be destroyed.
						  // Instantiated on first use.
	return instance;
}

void RSPCSensor::Initialize()
{
	// Start streaming with default recommended configuration
	try {
		pipe.start();
	}
	catch (const std::exception ex)
	{
	}
}

//RSPCFrame
bool RSPCSensor::GetNextFrame(
	int* textureFormat, /* rs2_format enum */
	int* textureWidth,
	int* textureHeight,
	const void** texturePixels,
	int* pointsSize,
	float** pointsX,
	float** pointsY,
	float** pointsZ,
	float** textureU,
	float** textureV)
{
	Deallocate();
	*pointsSize = 0;

	// Wait for the next set of frames from the camera
	auto frames = pipe.wait_for_frames();

	auto color = frames.get_color_frame();

	// For cameras that don't have RGB sensor, we'll map the pointcloud to infrared instead of color
	if (!color)
		color = frames.get_infrared_frame();

	// Tell pointcloud object to map to this color frame
	pc.map_to(color);

	auto depth = frames.get_depth_frame();

	// Generate the pointcloud and texture mappings
	points = pc.calculate(depth);
	if (!color)
	{
		return false;
	}

	*textureFormat = static_cast<int>(color.get_profile().format());
	*textureWidth = color.get_width();
	*textureHeight = color.get_height();
	*texturePixels = color.get_data();
	if (!points)
	{
		return false;
	}

	*pointsSize = points.size();
	const rs2::vertex* v = points.get_vertices();
	const rs2::texture_coordinate* tc = points.get_texture_coordinates();
	_pointsX = *pointsX = new float[*pointsSize];
	_pointsY = *pointsY = new float[*pointsSize];
	_pointsZ = *pointsZ = new float[*pointsSize];
	_textureU = *textureU = new float[*pointsSize];
	_textureV = *textureV = new float[*pointsSize];

	for (int i = 0; i < *pointsSize; ++i)
	{
		(*pointsX)[i] = v[i].x;
		(*pointsY)[i] = v[i].y;
		(*pointsZ)[i] = v[i].z;
		(*textureU)[i] = tc[i].u;
		(*textureV)[i] = tc[i].v;
	}

	return true;
}

#pragma once

#include "RSPCSensor.h"

class FrameExtratManager 
{
	 

public:
	void ExtractNextFrame()
	{ 
		RSPCSensor::Instance().GetNextFrame();

	}



};
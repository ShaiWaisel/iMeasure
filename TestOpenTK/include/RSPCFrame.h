#pragma once

class Frame {
public:

	Frame(int textheight,
		int texwidth,
		const void* texPixels,
		int pointsSize,
		const float* pointsX,
		const float* pointsY,
		const float* pointsZ,
		const float textureU,
		const float textureV): 
		_textheight(textheight),
		_texwidth(texwidth),
		_texPixels(texPixels),
		_pointsSize(pointsSize),
		_pointsX(pointsX),
		_pointsY(pointsY),
		_pointsZ(pointsZ),
		_textureU(textureU),
		_textureV(textureV)
	{
	}

private:

	int _textheight;
	int _texwidth;
	const void* _texPixels;
	int _pointsSize;
	const float* _pointsX;
	const float* _pointsY; 
	const float* _pointsZ; 
	const float _textureU; 
	const float _textureV;
};
#ifndef EXTRA_FUNCTIONS
#define EXTRA_FUNCTIONS

struct ColorRGB
{
	float red;
	float green;
	float blue;
};

struct ColorHSV
{
	int hue; //degrees [0, 360]
	float saturation;
	float value;
};

float power(float base, int exponent);
ColorHSV ConvertFromRGBToHSV(ColorRGB);
ColorRGB ConvertFromHSVToRGB(ColorHSV);

#endif
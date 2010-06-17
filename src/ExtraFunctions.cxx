#include "vtkKWProstateErrorMapRenderingWidget.h"

//These are useful functions which aren't tied to any
//particular class

float power(float base, int exponent)
{
	float result = 1;

	for(int i=0; i<exponent; i++)
	{
		result *= base;
	}

	return result;
}

ColorHSV ConvertFromRGBToHSV(ColorRGB colorRGB)
{
	ColorHSV colorHSV;

	float max = colorRGB.red;
	float min = colorRGB.red;

	if(colorRGB.green > max)
	{
		max = colorRGB.green;
	}
	if(colorRGB.green < min)
	{
		min = colorRGB.green;
	}

	if(colorRGB.blue > max)
	{
		max = colorRGB.blue;
	}
	if(colorRGB.blue < min)
	{
		min = colorRGB.blue;
	}

	//find hue
	if(max == min)
	{
		colorHSV.hue = 0;
	}
	else if(max == colorRGB.red)
	{
		colorHSV.hue = ((int)(60 * ((float)(colorRGB.green - colorRGB.blue)/(float)(max - min)) + 360))%360;
	}
	else if(max == colorRGB.green)
	{
		colorHSV.hue = (int)(60 * ((float)(colorRGB.blue - colorRGB.green)/(float)(max - min)) + 120);
	}
	else //max == colorRGB.blue
	{
		colorHSV.hue = (int)(60 * ((float)(colorRGB.red - colorRGB.green)/(float)(max - min)) + 240);
	}

	//find saturation
	if(max == 0)
	{
		colorHSV.saturation = 0;
	}
	else
	{
		colorHSV.saturation = (float)(max - min)/(float)min;
	}

	//find value
	colorHSV.value = max;

	return colorHSV;
}

ColorRGB ConvertFromHSVToRGB(ColorHSV colorHSV)
{
	ColorRGB colorRGB;

	int h2 = ((int)floor((float)colorHSV.hue/(float)60))%6;
	float f = ((float)colorHSV.hue/(float)60) - (int)floor((float)colorHSV.hue/(float)60);
	float v = colorHSV.value;
	float p = v * (1 - colorHSV.saturation);
	float q = v * (1 - f * colorHSV.saturation);
	float t = v * (1 - (1 - f)*colorHSV.saturation);

	if(h2 == 0)
	{
		colorRGB.red = v;
		colorRGB.green = t;
		colorRGB.blue = p;
	}
	else if(h2 == 1)
	{
		colorRGB.red = q;
		colorRGB.green = v;
		colorRGB.blue = p;
	}
	else if(h2 == 2)
	{
		colorRGB.red = p;
		colorRGB.green = v;
		colorRGB.blue = t;
	}
	else if(h2 == 3)
	{
		colorRGB.red = p;
		colorRGB.green = q;
		colorRGB.blue = v;
	}
	else if(h2 == 4)
	{
		colorRGB.red = t;
		colorRGB.green = p;
		colorRGB.blue = v;
	}
	else //h2 == 5
	{
		colorRGB.red = v;
		colorRGB.green = p;
		colorRGB.blue = q;
	}

	return colorRGB;
}
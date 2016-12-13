// Image Class (Written by Benoit Huet 06/97)
// 
// 
//

#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include "image_util.h"
#include "Image.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef VERBOSE
#define VERBOSE 0
#endif

Image::Image(char * fname)
{
	unsigned char tmpChar;
	char line[256];
	int tmpInt;	
	int MaxValue;

	ifstream ImgFile(fname);
	if(!ImgFile) {
		cerr << "Problem reading file (not found)... Aborting" << endl;
		exit(-1);
	}

	ImgFile >> tmpChar;
	if(tmpChar!='P') {
		cerr << "This file is not a standard IMAGE file... Aborting" << endl;
		exit(-1);
	}

	ImgFile >> tmpInt;
	switch(tmpInt) {
	case mono_ascii: 
	  format = mono_ascii;
	  type = mono;
	  break;
	case gray_ascii: 
	  format = gray_ascii;
	  type = gray;
	  break;
	case color_ascii: 
	  format = color_ascii;
	  type = color;
	  break;
	case mono_raw: 
	  format = mono_raw;
	  type = mono;
	  break;
	case gray_raw: 
	  format = gray_raw;
	  type = gray;
	  break;
	case color_raw: 
	  format = color_raw;
	  type = color;
	  break;
	default: 
	  format = 0;
	  break;
	}

	while((ImgFile >> tmpChar) && (tmpChar=='#') && ImgFile.good()) {
	  ImgFile.getline(line,256);
	}
	ImgFile.putback(tmpChar);
	ImgFile >> width >> height;
	
	if(VERBOSE) cout << "Image info: width=" << width << " height=" << height << endl;
	
	if(format==mono_ascii || format==mono_raw)
		MaxValue = 1;
	else ImgFile >> MaxValue;

	if(VERBOSE) cout << "Image Info: MaxValue=" << MaxValue << endl;

	if(format==color_ascii || format==color_raw) {
	  Alloc2D(R,unsigned char,height,width);
	  Alloc2D(G,unsigned char,height,width);
	  Alloc2D(B,unsigned char,height,width);

	  if(format==color_ascii) {
	    double factor = 1;	   
	    if(MaxValue!=255)
	      factor = (double)((double)255/(double)MaxValue);

	    for(unsigned int y=0; y < height; y++)
	      for(unsigned int x=0; x < width; x++) {
		ImgFile >> tmpInt;
		if(tmpInt<256) 
		  R[y][x]=(char)((int)(factor*tmpInt));
		else {
		  cerr << "Warning: Value (" << tmpInt << ") should be either 0 or 255... Defaulting to 0" << endl; 
		  R[y][x]='0';
		} 
		ImgFile >> tmpInt;
		if(tmpInt<256) 
		  G[y][x]=(char)((int)factor*tmpInt);
		else {
		  cerr << "Warning: Value (" << tmpInt << ") should be either 0 or 255... Defaulting to 0" << endl; 
		  G[y][x]='0';
		} 
		ImgFile >> tmpInt;
		if(tmpInt<256) 
		  B[y][x]=(char)((int)factor*tmpInt);
		else {
		  cerr << "Warning: Value (" << tmpInt << ") should be either 0 or 255... Defaulting to 0" << endl; 
		  B[y][x]='0';
		} 
	      }
	  }
	  else {
	    if(VERBOSE) cout << "Reading Color Raw image" << endl;
	    ImgFile.get((char &) tmpChar);
	    for(unsigned int y=0; y < height; y++)
	      for(unsigned int x=0; x < width; x++) {
		ImgFile.get((char &) R[y][x]);
		ImgFile.get((char &) G[y][x]);
		ImgFile.get((char &) B[y][x]);
	      }
	  }
	}
	else {
	  Alloc2D(R,unsigned char,height,width);
	  G = R;
	  B = R;

	  if(format==mono_ascii) {
	    if(VERBOSE) cout << "Reading BW Ascii image" << endl;
	    // 0 -> white and 1 -> black (file format!)
	    for(unsigned int y=0; y < height; y++)
	      for(unsigned int x=0; x < width; x++) {
		ImgFile >> tmpInt;
		if(tmpInt==0)
		  R[y][x]=(unsigned char)255;
		else if(tmpInt!=1) {
		  cerr << "Warning: Value (" << tmpInt << ") should be either 0 or 1... Defaulting to 1" << endl; 
		  R[y][x]=(unsigned char)0;
		}
	      }
	  }
	  else if(format==mono_raw) {
	    if(VERBOSE) cout << "Reading BW Raw image" << endl;
	    for(unsigned int y=0; y < height; y++)
	      for(unsigned int x=0; x < width;x+=8) {
		ImgFile.get((char &) tmpChar);
		//ImgFile >> tmpChar;
		//char dummy;
		//cin >> dummy;

		//cout << endl << (int)tmpChar << "--";
		unsigned char mask=128;
		for(unsigned int bit=0; bit < 8 && (x+bit)<width; bit++, mask/=2)
		  {
		    R[y][x+bit] = (tmpChar & mask) ? (unsigned char)0 : (unsigned char)255;
		    //(tmpChar & mask) ? cout << " 0" : cout << " 1";
		  }
	      }
	  }
	  else if(format==gray_ascii) {
	    if(VERBOSE) cout << "Reading Gray Ascii image" << endl;
	    double factor = (double)255/(double)MaxValue;

	    for(unsigned int y=0; y < height; y++)
	      for(unsigned int x=0; x < width; x++) {
		ImgFile >> tmpInt;
		if(tmpInt<256) 
		  R[y][x]=(char)((int)(factor*tmpInt));
		else {
		  cerr << "Warning: Value (" << tmpInt << ") should be either 0 or 255... Defaulting to 0" << endl; 
		  R[y][x]='0';
		}
	      }
	  }
	  else if(format==gray_raw) {
	    if(VERBOSE) cout << "Reading Gray Raw image" << endl;
	    ImgFile.get((char &) tmpChar);
	    for(unsigned int y=0; y < height; y++)
	      for(unsigned int x=0; x < width; x++) {
		ImgFile.get((char &) R[y][x]);
	      }
	  }
	}
	ImgFile.close();
	filename = strdup(fname);
}

Image::Image(unsigned int Width,unsigned int Height,int Type)
{
  filename = NULL;
  width = Width;
  height = Height;
  type = Type;
  format = unknown;

  if(type==color) {
    Alloc2D(R,unsigned char,height,width);
    Alloc2D(G,unsigned char,height,width);
    Alloc2D(B,unsigned char,height,width);
  }
  else {
    Alloc2D(R,unsigned char,height,width);
    G = R;
    B = R;
  }
}

Image::~Image()
{
	if(type==color) {
		if(R) Free2D(R);
		if(G) Free2D(G);
		if(B) Free2D(B);
	}
	else {
		if(R) Free2D(R);
		G = NULL;
		B = NULL;
	}

	if(filename) delete filename;
  
}

int Image::read(char * fname)
{
  // DOES NOTHING!
	return 1;
}

int Image::write(int saveformat, char * fname)
{
  static char * ImageInfo = "# Created using EyeViews (Written by Benoit Huet 06/97)";
  int MaxLineLenght =70;
  int linepos=0;

  ofstream ImgFile(fname);
  if(!ImgFile) {
    cerr << "Error" << endl;
    exit(-1);
  }
  if(saveformat==mono_ascii)
    {
      ImgFile << "P1" << endl << ImageInfo << endl << width << " " << height << endl;
      for(unsigned int y=0; y < height; y++)
	for(unsigned int x=0; x < width; x++) {
	  if(R[y][x]<127) ImgFile << 1;
	  else ImgFile << 0;
	  if(++linepos % (MaxLineLenght/2))
	    ImgFile << " ";
	  else
	    ImgFile << endl;
	}
    }
  else if(saveformat==mono_raw)
    {
      ImgFile << "P4" << endl << ImageInfo << endl << width << " " << height << endl;
      for(unsigned int y=0; y < height; y++)
	for(unsigned int x=0; x < width;) {      
	  char tmpChar=0;
	  for(int mask=128, bit=0; x < width && bit < 8; mask=mask/2, bit++,  x++)
	    if(R[y][x]<127) tmpChar |= mask;
	  
	  ImgFile << tmpChar; 
	  //if(!(++linepos % MaxLineLenght))
	  //  ImgFile << endl;
	}
    }
  else if(saveformat==gray_ascii) {
    ImgFile << "P2" << endl << ImageInfo << endl << width << " " << height << endl << 255 << endl;
      for(unsigned int y=0; y < height; y++)
	for(unsigned int x=0; x < width; x++) {
	  ImgFile.width(3);
	  ImgFile << (int)R[y][x];
	  linepos+=3;
	  if((((double) linepos+4)/(double)MaxLineLenght)>1)
	    {ImgFile << endl;linepos=0;}
	  else {ImgFile << " "; linepos++;}
	}
  }
  else if(saveformat==gray_raw) {
    ImgFile << "P5" << endl << ImageInfo << endl << width << " " << height << endl << 255 << endl;
      for(unsigned int y=0; y < height; y++)
	  for(unsigned int x=0; x < width; x++)
	    ImgFile << R[y][x];
  }
  else if(saveformat==color_ascii) {
    ImgFile << "P3" << endl << ImageInfo << endl << width << " " << height << endl << 255 << endl;
      for(unsigned int y=0; y < height; y++)
	for(unsigned int x=0; x < width; x++) {
	  ImgFile.width(3);
	  ImgFile << (int)R[y][x];
	  ImgFile << " ";
	  ImgFile.width(3);
	  ImgFile << (int)G[y][x];
	  ImgFile << " ";
	  ImgFile.width(3);
	  ImgFile << (int)B[y][x];
	  linepos+=11;
	  if((((double) linepos+12)/(double)MaxLineLenght)>1)
	    {ImgFile << endl;linepos=0;}
	  else {ImgFile << " "; linepos++;}
	}
  }
  else { //saveformat==color_raw
    ImgFile << "P6" << endl << ImageInfo << endl << width << " " << height << endl << 255 << endl;
      for(unsigned int y=0; y < height; y++)
	for(unsigned int x=0; x < width; x++) {
	  ImgFile << R[y][x];
	  ImgFile << G[y][x];
	  ImgFile << B[y][x];
	}
  }
  ImgFile.close();
  return(1);
}

//
// Routines doing linear interpolation for pixel values 
// at floating point position
//

unsigned char Image::get(float x, float y)
{
unsigned int x1=(unsigned int)x,x2=x1+1; if (x2>=width) x2=width-1;
unsigned int y1=(unsigned int)y,y2=y1+1; if (y2>=height) y2=height-1;
//if(x1==3 && y1==0)
// {
//   cout << "x1=" << x1 << " x2=" << x2 << " -- " << " y1=" << y1 << " y2=" << y2 << endl;
// }
float tx=x-(float)x1,itx=1.0f-tx;
float tmp1=((float)get(x1,y1)*itx)+((float)get(x2,y1)*tx);
float tmp2=((float)get(x1,y2)*itx)+((float)get(x2,y2)*tx);
float ty=y-(float)y1;

// if(x1==3 && y1==0)
//   {
//   cout << x << " - " << y <<" - TROUBLE - tx=" << tx << " itx=" << itx << " ty=" << ty << " return="<< (float)(tmp1*(1.0f-ty)+tmp2*ty) << endl;

//   cout << "1=" << (int)get(x1,y1) << " 2=" << (int)get(x2,y1) << " 3=" << (int)get(x1,y2)<< " 4=" << (int)get(x2,y2) << endl;
//   }
return ((unsigned char)(tmp1*(1.0f-ty)+tmp2*ty));
}

unsigned char Image::getR(float x, float y)
{
unsigned int x1=(int)x,x2=x1+1; if (x2>=width) x2=width-1;
unsigned int y1=(int)y,y2=y1+1; if (y2>=height) y2=height-1;

float tx=x-(float)x1,itx=1.0f-tx;
float tmp1=((float)getR(x1,y1)*itx)+((float)getR(x2,y1)*tx);
float tmp2=((float)getR(x1,y2)*itx)+((float)getR(x2,y2)*tx);
float ty=y-(float)y1;

return ((unsigned char)(tmp1*(1.0f-ty)+tmp2*ty));
}

unsigned char Image::getG(float x, float y)
{
unsigned int x1=(int)x,x2=x1+1; if (x2>=width) x2=width-1;
unsigned int y1=(int)y,y2=y1+1; if (y2>=height) y2=height-1;

float tx=x-(float)x1,itx=1.0f-tx;
float tmp1=((float)getG(x1,y1)*itx)+((float)getG(x2,y1)*tx);
float tmp2=((float)getG(x1,y2)*itx)+((float)getG(x2,y2)*tx);
float ty=y-(float)y1;

return ((unsigned char)(tmp1*(1.0f-ty)+tmp2*ty));
}

unsigned char Image::getB(float x, float y)
{
unsigned int x1=(int)x,x2=x1+1; if (x2>=width) x2=width-1;
unsigned int y1=(int)y,y2=y1+1; if (y2>=height) y2=height-1;

float tx=x-(float)x1,itx=1.0f-tx;
float tmp1=((float)getB(x1,y1)*itx)+((float)getB(x2,y1)*tx);
float tmp2=((float)getB(x1,y2)*itx)+((float)getB(x2,y2)*tx);
float ty=y-(float)y1;

return ((unsigned char)(tmp1*(1.0f-ty)+tmp2*ty));
}

#if 0
//
// rotate the image by a given angle with respect to the image 
// position centerX, centerY.
// the mode indicates whether the image dimensions are kept 
// or if they are recomputed to accomodate the whole image
//
void Image::rotate(double Angle, unsigned int centerX, unsigned int centerY, int mode)
{
  unsigned int tmpWidth;
  unsigned int tmpHeight;
  unsigned int newcenterX;
  unsigned int newcenterY;
  unsigned char ** tmpR;
  unsigned char ** tmpG;
  unsigned char ** tmpB;

  double sinAngle = sin(Angle);
  double cosAngle = cos(Angle);

  if(mode==0) { // image size unchanged during rotation
    tmpWidth = width;
    tmpHeight = height;
    newcenterX = centerX;
    newcenterY = centerY;
  }
  else {
     if(cosAngle==0) { // rotating by +/- 90 degrees
       tmpWidth = height;
       tmpHeight = width;
       newcenterX = tmpWidth/2;
       newcenterY = tmpHeight/2;
     }
     else if(Angle==M_PI) {
       tmpWidth = width;
       tmpHeight = height;
       newcenterX = centerX;
       newcenterY = centerY;      
     }
     else {
       
       double trslY = (-(centerX*sinAngle)-(centerY*cosAngle)+centerY);
       double trslX = (-(centerX*cosAngle)+(centerY*sinAngle)+centerX);
       int maxX,minX,maxY,minY;
       
       if((cosAngle>=0) && (sinAngle>=0)) {
	 maxX = (int)((cosAngle*width)-(sinAngle*0)+trslX);
	 minX = (int)((cosAngle*0)-(sinAngle*height)+trslX);
	 maxY = (int)((cosAngle*height)+(sinAngle*width)+trslY);
	 minY = (int)trslY;
       }
       else if(cosAngle<=0 && sinAngle>=0) {
	 maxX = (int)trslX;
	 minX = (int)((cosAngle*width)-(sinAngle*height)+trslX);
	 maxY = (int)((sinAngle*width)+trslY);
	 minY = (int)((cosAngle*height)+trslY);
       }
       else if(cosAngle<=0 && sinAngle<=0) {
	 maxX = (int)(-(sinAngle*height)+trslX);
	 minX = (int)((cosAngle*width)+trslX);
	 maxY = (int)trslY;
	 minY = (int)((sinAngle*width)+(cosAngle*height)+trslY);
       }
       else if(cosAngle>=0 && sinAngle<=0) {
	 maxX = (int)((cosAngle*width)-(sinAngle*height)+trslX);
	 minX = (int)trslX;
	 maxY = (int)((cosAngle*height)+trslY);
	 minY = (int)((sinAngle*width)+trslY);
       }
       
      tmpHeight = maxY - minY;
      tmpWidth = maxX - minX;
      
      newcenterX = (int)((cosAngle*centerX)-(sinAngle*centerY)+(-(centerX*cosAngle)+(centerY*sinAngle)+(tmpHeight/2)));
      newcenterY = (int)(+(cosAngle*centerY)+(sinAngle*centerX)+(-(centerX*sinAngle)-(centerY*cosAngle)+(tmpWidth/2)));
      
     }
  }

  if(type==color)
    {
      Alloc2D(tmpR,unsigned char,tmpHeight,tmpWidth);
      Alloc2D(tmpG,unsigned char,tmpHeight,tmpWidth);
      Alloc2D(tmpB,unsigned char,tmpHeight,tmpWidth);
    }
  else Alloc2D(tmpR,unsigned char,tmpHeight,tmpWidth);
  
  double invDET = 1.0/((cosAngle*cosAngle)+(sinAngle*sinAngle));
  double Xtranslation = (-sinAngle*(-(centerX*sinAngle)-(centerY*cosAngle)+newcenterY))-(cosAngle*(-(centerX*cosAngle)+(centerY*sinAngle)+newcenterX));
  double Ytranslation = (sinAngle*(-(centerX*cosAngle)+(centerY*sinAngle)+newcenterX))-(cosAngle*(-(centerX*sinAngle)-(centerY*cosAngle)+newcenterY));

  if(antiAliasing) {
    for (unsigned int ypos=0; ypos < tmpHeight; ypos++)
      for (unsigned int xpos=0; xpos < tmpWidth; xpos++)
        {
	float Rxpos = (float) invDET*((cosAngle*xpos)+(sinAngle*ypos)+Xtranslation);
	float Rypos = (float) invDET*((cosAngle*ypos)-(sinAngle*xpos)+Ytranslation);
	  if(Rxpos>=0 && Rxpos<width && Rypos>=0 && Rypos<height)
	    if(type == color)
	      {
		tmpR[ypos][xpos] = getR(Rxpos,Rypos);
		tmpG[ypos][xpos] = getG(Rxpos,Rypos);//(G[Rypos][Rxpos]+G[Rypos-1][Rxpos]+G[Rypos+1][Rxpos]+G[Rypos][Rxpos-1]+G[Rypos][Rxpos+1])/5;
		tmpB[ypos][xpos] = getB(Rxpos,Rypos);//(B[Rypos][Rxpos]+B[Rypos-1][Rxpos]+B[Rypos+1][Rxpos]+B[Rypos][Rxpos-1]+B[Rypos][Rxpos+1])/5;
	      }
	    else tmpR[ypos][xpos] = get(Rxpos,Rypos);
	    //if(tmpR[ypos][xpos]==215) cout << " -->" << xpos << " " << ypos << "---> " << Rxpos << "," << Rypos << endl;}//(unsigned char) ((int)((int)R[Rypos][Rxpos]+(int)R[Rypos-1][Rxpos]+(int)R[Rypos+1][Rxpos]+(int)R[Rypos][Rxpos-1]+(int)R[Rypos][Rxpos+1])/(int)5);
	  else
	    if(type == color)
	      {
		tmpR[ypos][xpos] = R[0][0];
		tmpG[ypos][xpos] = G[0][0];
		tmpB[ypos][xpos] = B[0][0];
	      }
	    else tmpR[ypos][xpos] = R[0][0];
	}
    }
  else {
    for (unsigned int ypos=0; ypos < tmpHeight; ypos++)
      for (unsigned int xpos=0; xpos < tmpWidth; xpos++)
        {
          int Rxpos = (int) invDET*((cosAngle*xpos)+(sinAngle*ypos)+Xtranslation);
	  int Rypos = (int) invDET*((cosAngle*ypos)-(sinAngle*xpos)+Ytranslation);
	  if(Rxpos>=0 && Rxpos<width && Rypos>=0 && Rypos<height)
	    if(type == color)
	      {
		tmpR[ypos][xpos] = R[Rypos][Rxpos];
		tmpG[ypos][xpos] = G[Rypos][Rxpos];
		tmpB[ypos][xpos] = B[Rypos][Rxpos];
	      }
	    else tmpR[ypos][xpos] = R[Rypos][Rxpos];
	  else 
	    if(type == color)
	      {
		tmpR[ypos][xpos] = R[0][0];
		tmpG[ypos][xpos] = G[0][0];
		tmpB[ypos][xpos] = B[0][0];
	      }
	    else tmpR[ypos][xpos] = R[0][0];
	}
	
      }

  width = tmpWidth;
  height = tmpHeight;
  unsigned char ** tmpChar = R;
  R = tmpR;
  Free2D(tmpChar);
  if(type == color) {
    tmpChar = G;
    G = tmpG;
    Free2D(tmpChar);
    tmpChar = B;
    B = tmpB;
    Free2D(tmpChar);
  }
  else {
    G = R;
    B = R;
  }
}
#endif

void Image::resize(unsigned int newWidth,unsigned int newHeight, int aspectRatio)
{

  unsigned char ** tmpR=NULL;
  unsigned char ** tmpG=NULL;
  unsigned char ** tmpB=NULL;

  double invRatioWidth = (double)width/(double)newWidth;
  double invRatioHeight = (double)height/(double)newHeight;


  if(aspectRatio==1 && invRatioWidth!=invRatioHeight) { // image aspect ratio unchanged during rotation
    if(invRatioWidth<invRatioHeight)
		invRatioWidth=invRatioHeight;
	else invRatioHeight=invRatioWidth;
  }

  if(type==color)
    {
      Alloc2D(tmpR,unsigned char,newHeight,newWidth);
      Alloc2D(tmpG,unsigned char,newHeight,newWidth);
      Alloc2D(tmpB,unsigned char,newHeight,newWidth);
    }
  else Alloc2D(tmpR,unsigned char,newHeight,newWidth);

  unsigned int ypos;
  unsigned int xpos;
  

// cout << "HERE" << endl;
 
   if(antiAliasing) {
// 	cout << "Doing AntiAliasing!" << endl;
// 	// first row
// 	Sypos = 0;
// 	tmpR[0][0]= (R[0][0]+R[1][0]+R[0][1])/3;
// 	if(type == color) {
// 		tmpG[0][0]= (G[0][0]+G[1][0]+G[0][1])/3;
// 		tmpB[0][0]= (B[0][0]+B[1][0]+B[0][1])/3;
// 	}	
// 	for (xpos=1; xpos < newWidth-1; xpos++) {
// 		Sxpos = (int) (invRatioWidth*xpos);
// 		tmpR[0][xpos] = (R[Sypos][Sxpos]+R[Sypos+1][Sxpos]+R[Sypos][Sxpos-1]+R[Sypos][Sxpos+1])/4;
// 		if(type == color) {
// 			tmpG[0][xpos] = (G[Sypos][Sxpos]+G[Sypos+1][Sxpos]+G[Sypos][Sxpos-1]+G[Sypos][Sxpos+1])/4;
// 			tmpB[0][xpos] = (B[Sypos][Sxpos]+B[Sypos+1][Sxpos]+B[Sypos][Sxpos-1]+B[Sypos][Sxpos+1])/4;
// 		}
// 	}
// 	Sxpos = (int) (invRatioWidth*(newWidth-1));
// 	tmpR[0][newWidth-1] = (R[Sypos][Sxpos]+R[Sypos+1][Sxpos]+R[Sypos][Sxpos-1])/3;
// 	if(type == color) {
// 		tmpG[0][newWidth-1] = (G[Sypos][Sxpos]+G[Sypos+1][Sxpos]+G[Sypos][Sxpos-1])/3;
// 		tmpB[0][newWidth-1] = (B[Sypos][Sxpos]+B[Sypos+1][Sxpos]+B[Sypos][Sxpos-1])/3;
// 	}
// cout << "first row... done" << endl;
// 	// last row
// 	Sypos = (int) invRatioHeight*(newHeight-1);
// 	tmpR[newHeight-1][0] = (R[Sypos][0]+R[Sypos-1][0]+R[Sypos][1])/3;
// 	if(type == color) {
// 		tmpG[newHeight-1][0] = (G[Sypos][0]+G[Sypos-1][0]+G[Sypos][1])/3;
// 		tmpB[newHeight-1][0] = (B[Sypos][0]+B[Sypos-1][0]+B[Sypos][1])/3;
// 	}
// 	tmpR[newHeight-1][newWidth-1] = (R[Sypos][Sxpos]+R[Sypos-1][Sxpos]+R[Sypos][Sxpos-1])/3;
// 	if(type == color) {
// 		tmpG[newHeight-1][newWidth-1] = (G[Sypos][Sxpos]+G[Sypos-1][Sxpos]+G[Sypos][Sxpos-1])/3;
// 		tmpB[newHeight-1][newWidth-1] = (B[Sypos][Sxpos]+B[Sypos-1][Sxpos]+B[Sypos][Sxpos-1])/3;
// 	}
// 	for (xpos=1; xpos < newWidth-1; xpos++) {
// 		Sxpos = (int) (invRatioWidth*xpos);
// 		tmpR[newHeight-1][xpos] = (R[Sypos][Sxpos]+R[Sypos-1][Sxpos]+R[Sypos][Sxpos-1]+R[Sypos][Sxpos+1])/4;
// 		if(type == color) {
// 			tmpG[newHeight-1][xpos] = (G[Sypos][Sxpos]+G[Sypos-1][Sxpos]+G[Sypos][Sxpos-1]+G[Sypos][Sxpos+1])/4;
// 			tmpB[newHeight-1][xpos] = (B[Sypos][Sxpos]+B[Sypos-1][Sxpos]+B[Sypos][Sxpos-1]+B[Sypos][Sxpos+1])/4;
// 		}
// 	}
// cout << "last row... done" << endl;	
// 	// first col
// 	Sxpos = 0;
// 	for (ypos=1; ypos < newHeight-1; ypos++) {
// 		Sypos = (int) (invRatioHeight*ypos);
// 		tmpR[ypos][0] = (R[Sypos][Sxpos]+R[Sypos+1][Sxpos]+R[Sypos-1][Sxpos]+R[Sypos][Sxpos+1])/4;
// 		if(type == color) {
// 			tmpG[ypos][0] = (G[Sypos][Sxpos]+G[Sypos+1][Sxpos]+G[Sypos-1][Sxpos]+G[Sypos][Sxpos+1])/4;
// 			tmpB[ypos][0] = (B[Sypos][Sxpos]+B[Sypos+1][Sxpos]+B[Sypos-1][Sxpos]+B[Sypos][Sxpos+1])/4;
// 		}
// 	}
// cout << "first col... done" << endl;
// 	// last col
// 	Sxpos = (int) invRatioWidth*(newWidth-1);
// 	for (ypos=1; ypos < newHeight-1; ypos++) {
// 		Sypos = (int) invRatioHeight*ypos;
// 		tmpR[ypos][newWidth-1] = (R[Sypos][Sxpos]+R[Sypos+1][Sxpos]+R[Sypos-1][Sxpos]+R[Sypos][Sxpos-1])/4;
// 		if(type == color) {
// 			tmpG[ypos][newWidth-1] = (G[Sypos][Sxpos]+G[Sypos+1][Sxpos]+G[Sypos-1][Sxpos]+G[Sypos][Sxpos-1])/4;
// 			tmpB[ypos][newWidth-1] = (B[Sypos][Sxpos]+B[Sypos+1][Sxpos]+B[Sypos-1][Sxpos]+B[Sypos][Sxpos-1])/4;
// 		}
// 	}
// cout << "last col... done" << endl;
// 	// rest of the image!
// 	for (ypos=1; ypos < newHeight-1; ypos++) {
// 		Sypos = (int) invRatioHeight*ypos;
// 		for (xpos=1; xpos < newWidth-1; xpos++) {
// 			Sxpos = (int) (invRatioWidth*xpos);
// 			tmpR[ypos][xpos] = (unsigned char) ((R[Sypos][Sxpos]+R[Sypos+1][Sxpos]+R[Sypos-1][Sxpos]+R[Sypos][Sxpos-1]+R[Sypos][Sxpos+1])/5);
// 			if(type == color) {
// 				tmpG[ypos][xpos] = (G[Sypos][Sxpos]+G[Sypos+1][Sxpos]+G[Sypos-1][Sxpos]+G[Sypos][Sxpos-1]+G[Sypos][Sxpos+1])/5;
// 				tmpB[ypos][xpos] = (B[Sypos][Sxpos]+B[Sypos+1][Sxpos]+B[Sypos-1][Sxpos]+B[Sypos][Sxpos-1]+B[Sypos][Sxpos+1])/5;
// 			}
// 		}
// 	}
// cout << "rest... done" << endl;
	float Sxpos,Sypos;
	for (ypos=0; ypos < newHeight; ypos++) {
 		Sypos = (float) invRatioHeight*ypos;
 		for (xpos=0; xpos < newWidth; xpos++) {
 			Sxpos = (float) (invRatioWidth*xpos);
 			tmpR[ypos][xpos] = get(Sxpos,Sypos);
 			if(type == color) {
				tmpG[ypos][xpos] = getG(Sxpos,Sypos);
				tmpB[ypos][xpos] = getB(Sxpos,Sypos);			}
 		}
 	}
  }
  else {
  int Sxpos;
  int Sypos;
  //cout << "No doing AntiAliasing!" << endl;
	for (ypos=0; ypos < newHeight; ypos++) {
		Sypos = (int) (invRatioHeight*ypos);
		//cout << Sypos << endl;
		for (xpos=0; xpos < newWidth; xpos++)
		{
			Sxpos = (int) (invRatioWidth*xpos);		
			//cout << " " << Sxpos << endl;
			tmpR[ypos][xpos] = R[Sypos][Sxpos];	
			if(type == color) {
				tmpG[ypos][xpos] = G[Sypos][Sxpos];
				tmpB[ypos][xpos] = B[Sypos][Sxpos];
			}
		}
	}
  }

   //cout << "THERE" << endl;

  width = newWidth;
  height = newHeight;
  unsigned char ** tmpChar = R;
  R = tmpR;
  Free2D(tmpChar);
  if(type == color) {
    tmpChar = G;
    G = tmpG;
    Free2D(tmpChar);
    tmpChar = B;
    B = tmpB;
    Free2D(tmpChar);
  }
  else {
    G = R;
    B = R;
  }
}


void Image::drawline(unsigned int x_start, unsigned int y_start, unsigned int x_end, unsigned int y_end, unsigned char gray)
{
	int dx = x_end-x_start;
	int dy = y_end-y_start;
	int steps,k;
	double xInc,yInc,x=x_start,y=y_start;

	if (abs(dx) > abs(dy)) steps = abs(dx);
	else steps = abs(dy);
	xInc= dx/ (double)steps;
	yInc= dy/ (double)steps;

	R[(unsigned int)y][(unsigned int)x]=gray;
	for(k=0; k<steps ; k++) {
		x+=xInc;
		y+=yInc;
		R[(unsigned int)y][(unsigned int)x]=gray;
	}
}

void Image::drawline(unsigned int x_start, unsigned int y_start, unsigned int x_end, unsigned int y_end, unsigned char red, unsigned char green, unsigned char blue)
{
	int dx = x_end-x_start;
	int dy = y_end-y_start;
	int steps,k;
	double xInc,yInc,x=x_start,y=y_start;

	if (abs(dx) > abs(dy)) steps = abs(dx);
	else steps = abs(dy);
	xInc= dx/ (double)steps;
	yInc= dy/ (double)steps;

	R[(unsigned int)y][(unsigned int)x]=red;
	G[(unsigned int)y][(unsigned int)x]=green;
	B[(unsigned int)y][(unsigned int)x]=blue;
	for(k=0; k<steps ; k++) {
		x=(xInc*k)+x_start;
		y=(yInc*k)+y_start;
		R[(unsigned int)y][(unsigned int)x]=red;
		G[(unsigned int)y][(unsigned int)x]=green;
		B[(unsigned int)y][(unsigned int)x]=blue;
	}
}



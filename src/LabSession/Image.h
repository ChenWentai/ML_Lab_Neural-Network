 // Image Class (Written by Benoit Huet 06/97)
// 
// 
//

#ifndef IMAGE_H
#define IMAGE_H

#include <cstdlib>
#include <iostream>

using namespace std;

class Image {

public:
  enum imageType { mono , gray, color };
  enum imageFormat { unknown, mono_ascii, gray_ascii, color_ascii, 
		     mono_raw, gray_raw, color_raw };

  Image(char * fname);
  Image(unsigned int width,unsigned int height,int type);
  ~Image();

  int read(char * fname);

  int write() {return write(format,filename);};
  int write(char * fname) {return write(format,fname);};
  int write(int saveformat, char * fname);

  void setAntiAliasing(int Val) {antiAliasing = Val ? 1 : 0 ; };

  void put(unsigned int x,unsigned int y,unsigned char gray) {R[y][x]=gray;};
  void put(unsigned int x,unsigned int y,unsigned char red,unsigned char green,unsigned  char blue) {R[y][x]=red;G[y][x]=green;B[y][x]=blue;};
  void putRed(unsigned int x,unsigned int y,unsigned char red) {R[y][x]=red;};
  void putGreen(unsigned int x,unsigned int y,unsigned char green){G[y][x]=green;};
  void putBlue(unsigned int x,unsigned int y,unsigned char blue) {B[y][x]=blue;};
  unsigned char get(unsigned int x,unsigned int y) {return R[y][x];};
  unsigned char getR(unsigned int x,unsigned int y) {return R[y][x];};
  unsigned char getG(unsigned int x,unsigned int y) {return G[y][x];};
  unsigned char getB(unsigned int x,unsigned int y) {return B[y][x];};
  
  unsigned char get(float x, float y);
  unsigned char getR(float x, float y);
  unsigned char getG(float x, float y);
  unsigned char getB(float x, float y);
  
  unsigned int getWidth() {return width;};
  unsigned int getHeight() {return height;};

  void drawline(unsigned int x_start,unsigned int y_start,unsigned int x_end,unsigned int y_end,unsigned char gray);
  void drawline(unsigned int x_start,unsigned int y_start,unsigned int x_end,unsigned int y_end,unsigned char red,unsigned char green,unsigned char blue);

  void resize(unsigned int newWidth,unsigned int newHeight, int aspectRatio = 0);
  //void resize(double ratioWidth, double ratioHeight, int aspectRatio = 0) {cout << "+" << (unsigned int)(width*ratioWidth)<< " - " << (unsigned int)(height*ratioHeight) << endl;resize((unsigned int)(width*ratioWidth),(unsigned int)(height*ratioHeight),aspectRatio);};
  void rotate(double angle,int mode) { rotate(angle,width/2,height/2,mode);};
  void rotate(double angle,unsigned int centerX,unsigned int centerY,int mode);

private:
  int type;
  int format;
  int antiAliasing; // 0=off 1=on
  char * filename;
  unsigned int width;
  unsigned int height;
  
  unsigned char ** R; // stores the red intensity 
  unsigned char ** G; // stores the green intensity
  unsigned char ** B; // stores the blue intensity
};

#endif






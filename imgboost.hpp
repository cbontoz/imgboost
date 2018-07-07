/**
 * \file imgboost.hpp
 * \brief Image processing library. It loads/exports tiff or raw data, in grb &
 * grayscale, provides numerous filters and eases calculations between
 * 2D-arrays.
 *
 *
 * MIT License
 *
 * Copyright (c) 2017-2018 Christos Bontozoglou
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * \author Christos Bontozoglou
 * \email cbontoz@gmail.com
 */

#ifndef _IMGBOOST_HPP_
#define _IMGBOOST_HPP_

#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <smartix.hpp>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef TIFF_EN_ // Defined in makefile
#include <tiffio.h>
#endif // TIFF_EN_
#include <vector>

class imgboost {
public:
  struct Pixel {
    int a, r, g, b;
    Pixel() : a(255), r(0), g(0), b(0) {}
    Pixel(int Gray) : a(255), r(Gray), g(Gray), b(Gray) {}
    Pixel(int R, int G, int B) : a(255), r(R), g(G), b(B) {}
    Pixel(int A, int R, int G, int B) : a(A), r(R), g(G), b(B) {}
  };
  size_t width() const { return _width; };
  size_t height() const { return _height; };
  imgboost();
  imgboost(imgboost const &src);
  template <typename T>
  imgboost(smartix<T> const &sc, size_t Width, size_t Height);
  imgboost(std::string path);
  imgboost(size_t Width, size_t Height);
  imgboost(size_t Width, size_t Height, int fill);

  void setPixel(size_t Column, size_t Row, Pixel p);
  void setPixel(size_t Column, size_t Row, int Gray);
  void setRed(size_t Column, size_t Row, int Red);
  void setGreen(size_t Column, size_t Row, int Green);
  void setBlue(size_t Column, size_t Row, int Blue);
  Pixel getPixel(size_t Column, size_t Row);
  double getScalar(size_t Column, size_t Row);
  int getGray(size_t Column, size_t Row);
  int getRed(size_t Column, size_t Row);
  int getGreen(size_t Column, size_t Row);
  int getBlue(size_t Column, size_t Row);

  void DrawLine(size_t xstr, size_t ystr, size_t xend, size_t yend,
                Pixel color);

  imgboost &operator=(imgboost const &src);
  template <typename T> imgboost &operator=(T const &src);
  imgboost &operator+=(imgboost const &src);
  imgboost &operator-=(imgboost const &src);
  imgboost &operator*=(imgboost const &src);
  imgboost &operator/=(imgboost const &src);
  template <typename T> imgboost &operator+=(T const rhand);
  template <typename T> imgboost &operator-=(T const rhand);
  template <typename T> imgboost &operator*=(T const rhand);
  template <typename T> imgboost &operator/=(T const rhand);
  template <typename T> imgboost operator+(T const &rhand);
  template <typename T> imgboost operator-(T const &rhand);
  template <typename T> imgboost operator*(T const &rhand);
  template <typename T> imgboost operator/(T const &rhand);
  template <typename T>
  friend imgboost operator+(T const &rhand, imgboost const &src);
  template <typename T>
  friend imgboost operator-(T const &rhand, imgboost const &src);
  template <typename T>
  friend imgboost operator*(T const &rhand, imgboost const &src);
  template <typename T>
  friend imgboost operator/(T const &rhand, imgboost const &src);

  void savecsv(std::string filename, int BytesPerPixel);
  void loadcsv(std::string filename);
#ifdef TIFF_EN_
  void loadtiff(std::string tiffpath);
  void loadtiff(const char *tiffpath);
  void savetiff(std::string tiffpath);
  void savetiff(const char *tiffpath);
#endif // TIFF_EN_

  double Mean();
  double Stdev();
  void Crop(size_t w_border, size_t h_border);
  imgboost RegionOfInterest(size_t Column, size_t Row, size_t s_size);
  imgboost RegionOfInterest(size_t Column, size_t Row, size_t w_size,
                            size_t h_size);
  void StitchRoIBack(imgboost RoI, size_t Column, size_t Row);
  void MirrorPadding(size_t ExtraBorder);
  void averageWith(imgboost const &me);

  size_t findcell(size_t Row, size_t Column) { return Row * _width + Column; }

  smartix<double> scalar() { return _scalar; }
  double scalar(size_t w, size_t h) { return getScalar(w, h); }

private:
  size_t _width;
  size_t _height;
  smartix<int> _r;
  smartix<int> _g;
  smartix<int> _b;
  smartix<int> _a;
  smartix<int> _gray;
  smartix<double> _scalar;

  void copy(imgboost const &copyme);
  std::vector<size_t> _dimVector(size_t Width, size_t Height);
  void _checkCoordinates(size_t Column, size_t Row);
  void _compareDimensions(size_t Width, size_t Height);
  void _resize(size_t newWidth, size_t newHeight);
};

// Constructores/Destructors
imgboost::imgboost()
    : _width(0), _height(0), _r(_dimVector(_width, _height)),
      _g(_dimVector(_width, _height)), _b(_dimVector(_width, _height)),
      _a(_dimVector(_width, _height)), _gray(_dimVector(_width, _height)),
      _scalar(_dimVector(_width, _height)) {}

imgboost::imgboost(imgboost const &src)
    : _width(src.width()), _height(src.height()),
      _r(_dimVector(_width, _height)), _g(_dimVector(_width, _height)),
      _b(_dimVector(_width, _height)), _a(_dimVector(_width, _height)),
      _gray(_dimVector(_width, _height)), _scalar(_dimVector(_width, _height)) {
  copy(src);
}

template <typename T>
imgboost::imgboost(smartix<T> const &sc, size_t Width, size_t Height)
    : _width(Width), _height(Height), _r(_dimVector(_width, _height)),
      _g(_dimVector(_width, _height)), _b(_dimVector(_width, _height)),
      _a(_dimVector(_width, _height)), _gray(_dimVector(_width, _height)),
      _scalar(_dimVector(_width, _height)) {
  for (size_t h = 0; h < Height; ++h)
    for (size_t w = 0; w < Width; ++w) {
      setPixel(w, h, Pixel(round(sc.uacc_const(h * Width + w))));
    }
}

imgboost::imgboost(std::string path)
    : _width(0), _height(0), _r(_dimVector(_width, _height)),
      _g(_dimVector(_width, _height)), _b(_dimVector(_width, _height)),
      _a(_dimVector(_width, _height)), _gray(_dimVector(_width, _height)),
      _scalar(_dimVector(_width, _height)) {

  unsigned int idx = path.rfind('.');
  std::string extension = path.substr(idx + 1);

  if (extension.compare("tiff") == 0 || extension.compare("tif") == 0) {
    loadtiff(path);
  } else if (extension.compare("csv") == 0) {
    loadcsv(path);
  }
  // else if ... add your extensions
  else {
    throw std::invalid_argument("Extension not supported.");
  }
}

imgboost::imgboost(size_t Width, size_t Height)
    : _width(Width), _height(Height), _r(_dimVector(_width, _height)),
      _g(_dimVector(_width, _height)), _b(_dimVector(_width, _height)),
      _a(_dimVector(_width, _height)), _gray(_dimVector(_width, _height)),
      _scalar(_dimVector(_width, _height)) {}

imgboost::imgboost(size_t Width, size_t Height, int fill)
    : _width(Width), _height(Height), _r(_dimVector(_width, _height)),
      _g(_dimVector(_width, _height)), _b(_dimVector(_width, _height)),
      _a(_dimVector(_width, _height)), _gray(_dimVector(_width, _height)),
      _scalar(_dimVector(_width, _height)) {
  for (size_t h = 0; h < _height; ++h)
    for (size_t w = 0; w < _width; ++w) {
      setPixel(w, h, Pixel(round(fill)));
    }
}

// Basic pixel manipulation
void imgboost::setPixel(size_t Column, size_t Row, Pixel p) {
  _checkCoordinates(Column, Row);
  size_t position = findcell(Row, Column);
  _a.uacc(position) = (p.a > 255) ? 255 : p.a;
  _r.uacc(position) = (p.r > 255) ? 255 : p.r;
  _g.uacc(position) = (p.g > 255) ? 255 : p.g;
  _b.uacc(position) = (p.b > 255) ? 255 : p.b;
  _scalar.uacc(position) = static_cast<double>(
      (_r.uacc(position) + _g.uacc(position) + _b.uacc(position)) / 3.0);
  _gray.uacc(position) = round(_scalar.uacc(position));
}

void imgboost::setPixel(size_t Column, size_t Row, int Gray) {
  setPixel(Column, Row, Pixel(Gray));
}

void imgboost::setRed(size_t Column, size_t Row, int Red) {
  size_t position = findcell(Row, Column);
  _r.uacc(position) = (Red > 255) ? 255 : Red;
  _scalar.uacc(position) = static_cast<double>(
      (_r.uacc(position) + _g.uacc(position) + _b.uacc(position)) / 3.0);
  _gray.uacc(position) = round(_scalar.uacc(position));
}

void imgboost::setGreen(size_t Column, size_t Row, int Green) {
  size_t position = findcell(Row, Column);
  _g.uacc(position) = (Green > 255) ? 255 : Green;
  _scalar.uacc(position) = static_cast<double>(
      (_r.uacc(position) + _g.uacc(position) + _b.uacc(position)) / 3.0);
  _gray.uacc(position) = round(_scalar.uacc(position));
}

void imgboost::setBlue(size_t Column, size_t Row, int Blue) {
  size_t position = findcell(Row, Column);
  _b.uacc(position) = (Blue > 255) ? 255 : Blue;
  _scalar.uacc(position) = static_cast<double>(
      (_r.uacc(position) + _g.uacc(position) + _b.uacc(position)) / 3.0);
  _gray.uacc(position) = round(_scalar.uacc(position));
}

imgboost::Pixel imgboost::getPixel(size_t Column, size_t Row) {
  _checkCoordinates(Column, Row);
  size_t position = findcell(Row, Column);
  return Pixel(_a.uacc(position), _r.uacc(position), _g.uacc(position),
               _b.uacc(position));
}

double imgboost::getScalar(size_t Column, size_t Row) {
  return _scalar.uacc(findcell(Row, Column));
}

int imgboost::getGray(size_t Column, size_t Row) {
  return _gray.uacc(findcell(Row, Column));
}

int imgboost::getRed(size_t Column, size_t Row) {
  return _r.uacc(findcell(Row, Column));
}

int imgboost::getGreen(size_t Column, size_t Row) {
  return _g.uacc(findcell(Row, Column));
}

int imgboost::getBlue(size_t Column, size_t Row) {
  return _b.uacc(findcell(Row, Column));
}

void imgboost::DrawLine(size_t xstr, size_t ystr, size_t xend, size_t yend,
                        Pixel color) {
  double x1 = xstr, x2 = xend, y1 = ystr, y2 = yend;
  double b, a;
  if (x1 != x2)
    b = (y1 - y2) / (x1 - x2);
  else
    b = 0;
  a = y1 - b * x1;

  for (double x = xstr; x < xend; ++x) {
    double y = round(b * x + a);
    setPixel(x, y, color);
  }
}

// Operators overload
imgboost &imgboost::operator=(imgboost const &src) {
  copy(src);
  return *this;
}

template <typename T> imgboost &imgboost::operator=(T const &rhand) {
  for (size_t h = 0; h < _height; ++h)
    for (size_t w = 0; w < _width; ++w) {
      setPixel(w, h, Pixel(round(rhand)));
    }
  return *this;
}

imgboost &imgboost::operator+=(imgboost const &rhand) {
  _compareDimensions(rhand.width(), rhand.height());
  _r += rhand._r;
  _g += rhand._g;
  _b += rhand._b;
  _gray += rhand._gray;
  _scalar += rhand._scalar;
  return *this;
}

template <typename T> imgboost &imgboost::operator+=(T const rhand) {
  _r += rhand;
  _g += rhand;
  _b += rhand;
  _gray += rhand;
  _scalar += rhand;
  return *this;
}

imgboost &imgboost::operator-=(imgboost const &rhand) {
  _compareDimensions(rhand.width(), rhand.height());
  _r -= rhand._r;
  _g -= rhand._g;
  _b -= rhand._b;
  _gray -= rhand._gray;
  _scalar -= rhand._scalar;
  return *this;
}

template <typename T> imgboost &imgboost::operator-=(T const rhand) {
  _r -= rhand;
  _g -= rhand;
  _b -= rhand;
  _gray -= rhand;
  _scalar -= rhand;
  return *this;
}

imgboost &imgboost::operator*=(imgboost const &rhand) {
  _compareDimensions(rhand.width(), rhand.height());
  _r *= rhand._r;
  _g *= rhand._g;
  _b *= rhand._b;
  _gray *= rhand._gray;
  _scalar *= rhand._scalar;
  return *this;
}

template <typename T> imgboost &imgboost::operator*=(T const rhand) {
  _r *= rhand;
  _g *= rhand;
  _b *= rhand;
  _gray *= rhand;
  _scalar *= rhand;
  return *this;
}

imgboost &imgboost::operator/=(imgboost const &rhand) {
  _compareDimensions(rhand.width(), rhand.height());
  _r /= rhand._r;
  _g /= rhand._g;
  _b /= rhand._b;
  _gray /= rhand._gray;
  _scalar /= rhand._scalar;
  return *this;
}

template <typename T> imgboost &imgboost::operator/=(T const rhand) {
  if (rhand == 0)
    throw std::overflow_error("You attempt to divide image by zero.");
  _r /= rhand;
  _g /= rhand;
  _b /= rhand;
  _gray /= rhand;
  _scalar /= rhand;
  return *this;
}

template <typename T> imgboost imgboost::operator+(T const &rhand) {
  imgboost cthis(*this);
  cthis += rhand;
  return cthis;
}

template <typename T> imgboost operator+(T const &rhand, imgboost const &src) {
  imgboost cthis(src.width(), src.height());
  cthis = rhand;
  cthis += src;
  return cthis;
}

template <typename T> imgboost imgboost::operator-(T const &rhand) {
  imgboost cthis(*this);
  cthis -= rhand;
  return cthis;
}

template <typename T> imgboost operator-(T const &rhand, imgboost const &src) {
  imgboost cthis(src.width(), src.height());
  cthis = rhand;
  cthis -= src;
  return cthis;
}

template <typename T> imgboost imgboost::operator*(T const &rhand) {
  imgboost cthis(*this);
  cthis *= rhand;
  return cthis;
}

template <typename T> imgboost operator*(T const &rhand, imgboost const &src) {
  imgboost cthis(src.width(), src.height());
  cthis = rhand;
  cthis *= src;
  return cthis;
}

template <typename T> imgboost imgboost::operator/(T const &rhand) {
  imgboost cthis(*this);
  cthis /= rhand;
  return cthis;
}

template <typename T> imgboost operator/(T const &rhand, imgboost const &src) {
  imgboost cthis(src.width(), src.height());
  cthis = rhand;
  cthis /= src;
  return cthis;
}

// Load & Save functions
void imgboost::savecsv(std::string filename, int BytesPerPixel) {
  std::ofstream myfile;
  std::string buffer;
  myfile.open(filename.c_str());
  if (BytesPerPixel < 0 || BytesPerPixel == 2 || BytesPerPixel > 4)
    throw std::invalid_argument(
        "Unkonwn format, invalid number of bytes per pixel.");
  myfile << "Width," << std::to_string(_width) << ",Height,"
         << std::to_string(_height) << ",Bytes/Pixel,"
         << std::to_string(BytesPerPixel) << "\n";
  if (BytesPerPixel == 4) {
    for (size_t h = 0; h < _height; ++h) {
      buffer = "";
      for (size_t w = 0; w < _width - 1; ++w) {
        buffer += std::to_string(_r.uacc(findcell(h, w))) + std::string(",") +
                  std::to_string(_g.uacc(findcell(h, w))) + std::string(",") +
                  std::to_string(_b.uacc(findcell(h, w))) + std::string(",") +
                  std::to_string(_a.uacc(findcell(h, w))) + std::string(",");
      }
      buffer +=
          std::to_string(_r.uacc(findcell(h, _width - 1))) + std::string(",") +
          std::to_string(_g.uacc(findcell(h, _width - 1))) + std::string(",") +
          std::to_string(_b.uacc(findcell(h, _width - 1))) + std::string(",") +
          std::to_string(_a.uacc(findcell(h, _width - 1))) + std::string("\n");
      myfile << buffer.c_str();
    }
  } else if (BytesPerPixel == 3) {
    for (size_t h = 0; h < _height; ++h) {
      buffer = "";
      for (size_t w = 0; w < _width - 1; ++w) {
        buffer += std::to_string(_r.uacc(findcell(h, w))) + std::string(",") +
                  std::to_string(_g.uacc(findcell(h, w))) + std::string(",") +
                  std::to_string(_b.uacc(findcell(h, w))) + std::string(",");
      }
      buffer +=
          std::to_string(_r.uacc(findcell(h, _width - 1))) + std::string(",") +
          std::to_string(_g.uacc(findcell(h, _width - 1))) + std::string(",") +
          std::to_string(_b.uacc(findcell(h, _width - 1))) + std::string("\n");
      myfile << buffer.c_str();
    }
  } else if (BytesPerPixel == 1) {
    for (size_t h = 0; h < _height; ++h) {
      buffer = "";
      for (size_t w = 0; w < _width - 1; ++w) {
        buffer += std::to_string(_gray.uacc(findcell(h, w))) + std::string(",");
      }
      buffer += std::to_string(_gray.uacc(findcell(h, _width - 1))) +
                std::string("\n");
      myfile << buffer.c_str();
    }
  }
  myfile.close();
}

void imgboost::loadcsv(std::string filename) {

  std::ifstream csvin(filename);
  std::string line;
  std::getline(csvin, line);
  std::stringstream lineStream(line);

  std::string cell;
  size_t Width, Height;
  int pxlbytes = -1;
  std::getline(lineStream, cell, ',');
  // Dump first cell, string.
  std::getline(lineStream, cell, ',');
  Width = stoul(cell);
  std::getline(lineStream, cell, ',');
  // Dump third cell, string.
  std::getline(lineStream, cell, ',');
  Height = stoul(cell);
  std::getline(lineStream, cell, ',');
  // Dump fifth cell, string.
  std::getline(lineStream, cell, ',');
  pxlbytes = stoi(cell);
  _resize(Width, Height);

  size_t h = 0, w = 0;
  if (pxlbytes == 4) {
    while (std::getline(csvin, line)) {
      std::stringstream lineStream(line);
      std::string rstr, gstr, bstr, astr;
      while (std::getline(lineStream, rstr, ',')) {
        std::getline(lineStream, gstr, ',');
        std::getline(lineStream, bstr, ',');
        std::getline(lineStream, astr, ',');
        setPixel(w, h, Pixel(stoi(astr), stoi(rstr), stoi(gstr), stoi(bstr)));
        ++w;
      }
      w = 0;
      ++h;
    }
  } else if (pxlbytes == 3) {
    while (std::getline(csvin, line)) {
      std::stringstream lineStream(line);
      std::string rstr, gstr, bstr;
      while (std::getline(lineStream, rstr, ',')) {
        std::getline(lineStream, gstr, ',');
        std::getline(lineStream, bstr, ',');
        setPixel(w, h, Pixel(stoi(rstr), stoi(gstr), stoi(bstr)));
        ++w;
      }
      w = 0;
      ++h;
    }
  } else if (pxlbytes == 1) {
    while (std::getline(csvin, line)) {
      std::stringstream lineStream(line);
      std::string graystr;
      while (std::getline(lineStream, graystr, ',')) {
        setPixel(w, h, Pixel(stoi(graystr)));
        ++w;
      }
      w = 0;
      ++h;
    }
  } else {
    throw std::invalid_argument("Wrong number of bytes per pixel in CSV.");
  }
}

#ifdef TIFF_EN_
void imgboost::loadtiff(const char *tiffpath) {
  loadtiff(std::string(tiffpath));
}

void imgboost::loadtiff(std::string tiffpath) {
  const char *mode = "r";
  size_t npixels;
  int width, height;
  TIFF *img = TIFFOpen(tiffpath.c_str(), mode);
  uint32_t *rgbraster;

  if (!img)
    throw std::invalid_argument("Failed to open the tiff image");

  TIFFGetField(img, TIFFTAG_IMAGEWIDTH, &width);
  TIFFGetField(img, TIFFTAG_IMAGELENGTH, &height);

  _resize(width, height);

  npixels = (width) * (height);
  rgbraster = (uint32_t *)_TIFFmalloc(npixels * sizeof(uint32_t));

  if (rgbraster == NULL)
    throw "Failed to load tiff format image";
  TIFFReadRGBAImage(img, width, height, rgbraster, 0);
  TIFFClose(img);

  for (int y = 0; y < height; ++y)
    for (int x = 0; x < width; ++x) {
      size_t ind = (height - y - 1) * width + x;
      setPixel(x, y, Pixel(TIFFGetA(rgbraster[ind]), TIFFGetR(rgbraster[ind]),
                           TIFFGetG(rgbraster[ind]), TIFFGetB(rgbraster[ind])));
    }

  _TIFFfree(rgbraster);
}

void imgboost::savetiff(const char *tiffpath) {
  savetiff(std::string(tiffpath));
}

void imgboost::savetiff(std::string tiffpath) {
  TIFF *outimage;
  uint8_t *data = (uint8_t *)malloc(_height * _width * sizeof(uint32_t));
  uint8_t *datacp = data;

  size_t ind, rgb_ind;
  for (size_t y = 0; y < _height; ++y) {
    for (size_t x = 0; x < _width; ++x) {
      ind = ((_height - y - 1) * _width + x) * 4;
      rgb_ind = findcell(y, x);
      data[ind] = (uint8_t)_r.uacc(rgb_ind);
      data[ind + 1] = (uint8_t)_g.uacc(rgb_ind);
      data[ind + 2] = (uint8_t)_b.uacc(rgb_ind);
      data[ind + 3] = (uint8_t)_a.uacc(rgb_ind);
    }
  }

  if ((outimage = TIFFOpen(tiffpath.c_str(), "w")) == NULL)
    throw "Failed to create tiff output file";

  TIFFSetField(outimage, TIFFTAG_IMAGEWIDTH, _width);
  TIFFSetField(outimage, TIFFTAG_IMAGELENGTH, _height);
  TIFFSetField(outimage, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(outimage, TIFFTAG_SAMPLESPERPIXEL, 4);
  TIFFSetField(outimage, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(outimage, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(outimage, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
  TIFFSetField(outimage, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

  tsize_t linebytes = 4 * _width;

  unsigned char *buf = NULL;
  buf = (unsigned char *)_TIFFmalloc(linebytes);

  TIFFSetField(outimage, TIFFTAG_ROWSPERSTRIP,
               TIFFDefaultStripSize(outimage, _width * 4));

  for (uint32_t row = 0; row < _height; row++) {
    memcpy(buf, &data[(row)*linebytes], linebytes);
    if (TIFFWriteScanline(outimage, buf, row, 0) < 0)
      break;
  }

  TIFFClose(outimage);
  data = datacp;
  free(data);
  data = NULL;
  datacp = NULL;
}
#endif // TIFF_EN_

// Image Processing
double imgboost::Mean() {
  double mean = 0;
  for (size_t h = 0; h < _height; ++h)
    for (size_t w = 0; w < _width; ++w) {
      mean += _scalar.uacc(findcell(h, w));
    }
  mean /= static_cast<double>(_height * _width);
  return mean;
}

double imgboost::Stdev() {
  double stdev = 0;
  double mean = Mean();
  for (size_t h = 0; h < _height; ++h)
    for (size_t w = 0; w < _width; ++w) {
      stdev += pow(_scalar.uacc(findcell(h, w)) - mean, 2);
    }
  stdev = sqrt(stdev / static_cast<double>(_width * _height));
  return stdev;
}

void imgboost::Crop(size_t w_border, size_t h_border) {
  imgboost cropped = this->RegionOfInterest(
      w_border, h_border, _width - 2 * w_border, _height - 2 * h_border);

  _resize(cropped.width(), cropped.height());
  copy(cropped);
}

imgboost imgboost::RegionOfInterest(size_t Column, size_t Row, size_t s_size) {
  return RegionOfInterest(Column, Row, s_size, s_size);
}

imgboost imgboost::RegionOfInterest(size_t Column, size_t Row, size_t w_size,
                                    size_t h_size) {
  if (Column >= _width || (Column + w_size) >= _width || Row >= _height ||
      (Row + h_size) >= _height)
    throw std::length_error("Region requested out of image dimensions.");

  imgboost roi(w_size, h_size);
  for (size_t y = Row; y < Row + h_size; ++y)
    for (size_t x = Column; x < Column + w_size; ++x) {
      roi.setPixel(x - Column, y - Row, getPixel(x, y));
    }
  return roi;
}

void imgboost::StitchRoIBack(imgboost RoI, size_t Column, size_t Row) {
  _checkCoordinates(Column, Row);
  for (size_t h = 0; h < RoI.height(); ++h)
    for (size_t w = 0; w < RoI.width(); ++w) {
      setPixel(Column + w, Row + h, RoI.getPixel(w, h));
    }
}

void imgboost::MirrorPadding(size_t ExtraBorder) {
  size_t eb = ExtraBorder;
  imgboost old = *this;
  _resize(old.width() + 2 * eb, old.height() + 2 * eb);
  StitchRoIBack(old, eb, eb);

  // Pad corners
  for (size_t h = 0; h < eb; ++h)
    for (size_t w = 0; w < eb; ++w) {
      setPixel(w, h, old.getPixel(eb - w - 1, eb - h - 1));
      setPixel(_width - w - 1, _height - h - 1,
               old.getPixel(old.width() - (eb - w) - 1,
                            old.height() - (eb - h) - 1));
      setPixel(_width - w - 1, h,
               old.getPixel(old.width() - (eb - w) - 1, eb - h));
      setPixel(w, _height - h - 1,
               old.getPixel(eb - w, old.height() - (eb - h) - 1));
    }
  // Pad horizontal
  for (size_t h = 0; h < eb; ++h)
    for (size_t w = eb; w < _width - eb; ++w) {
      setPixel(w, h, old.getPixel(w - eb, eb - h));
      setPixel(_width - w - 1, _height - h - 1,
               old.getPixel(old.width() - (w - eb) - 1,
                            old.height() - (eb - h) - 1));
    }
  // Pad vertical
  for (size_t h = eb; h < _height - eb; ++h)
    for (size_t w = 0; w < eb; ++w) {
      setPixel(w, h, old.getPixel(eb - w, h - eb));
      setPixel(_width - w - 1, _height - h - 1,
               old.getPixel(old.width() - (eb - w) - 1,
                            old.height() - (h - eb) - 1));
    }
}

void imgboost::averageWith(imgboost const &me) {
  _compareDimensions(me.width(), me.height());
  *this = (*this + me) / 2.0;
}

// Private & Support functions
void imgboost::copy(imgboost const &copyme) {
  _compareDimensions(copyme.width(), copyme.height());
  _a = copyme._a;
  _r = copyme._r;
  _g = copyme._g;
  _b = copyme._b;
  _gray = copyme._gray;
  _scalar = copyme._scalar;
}

std::vector<size_t> imgboost::_dimVector(size_t Width, size_t Height) {
  std::vector<size_t> dimensions;
  dimensions.push_back(Height);
  dimensions.push_back(Width);
  return dimensions;
}

void imgboost::_checkCoordinates(size_t Column, size_t Row) {
  (Column >= _width)
      ? throw std::length_error("Column coordinate should be < Width")
      : 0;
  (Row >= _height)
      ? throw std::length_error("Row coordinate should be < Height")
      : 0;
}

void imgboost::_compareDimensions(size_t Width, size_t Height) {
  (Width != _width) ? throw std::length_error("Widths do not match") : 0;
  (Height != _height) ? throw std::length_error("Heights do not match") : 0;
}

void imgboost::_resize(size_t newWidth, size_t newHeight) {
  _a.~smartix<int>();
  _r.~smartix<int>();
  _g.~smartix<int>();
  _b.~smartix<int>();
  _gray.~smartix<int>();
  _scalar.~smartix<double>();

  _width = newWidth;
  _height = newHeight;

  new (&_a) smartix<int>(_dimVector(_width, _height));
  new (&_r) smartix<int>(_dimVector(_width, _height));
  new (&_g) smartix<int>(_dimVector(_width, _height));
  new (&_b) smartix<int>(_dimVector(_width, _height));
  new (&_gray) smartix<int>(_dimVector(_width, _height));
  new (&_scalar) smartix<double>(_dimVector(_width, _height));
}

#endif //_IMGBOOST_HPP_

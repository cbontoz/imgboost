/**
 * \file ib_example.cpp
 * \brief Image processing library. It loads/exports tiff or raw data, in grb &
 * grayscale, provides numerous filters and eases calculations between
 * 2D-arrays.
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

#include <ccmatrix.hpp>
#include <chrono>
#include <gaussian.hpp>
#include <imgboost.hpp>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
  std::string tfilename;
  if (argc > 1)
    tfilename = argv[1];
  else
    tfilename = "./sample-data/butterfly.tif";
  // Load a tiff
  imgboost img1(tfilename);
  // Load the grayscale from .csv
  // imgboost img1("./sample-data/grayscale.csv");
  // Save a rectangular RoI from the loaded image
  img1.RegionOfInterest(50, 50, 30, 50).savetiff("./RoI_1.tiff");
  // Save image as raw data in csv file. Where 1 grayscale, 3 rgb and 4 rgba
  img1.savecsv("./sample-data/grayscale.csv", 1);

  // Load image from CSV file
  imgboost dummy("./sample-data/grayscale.csv");
  dummy.savetiff("./Grayscale.tiff");

  // Crop 10 pixels from each side and save
  img1.Crop(10, 10);
  img1.savetiff("./Cropped.tiff");

  gaussian blur;
  gaussian_fast blur_fast;
  // Make a blur copy using 2D gaussian filter kernel
  std::cout << "Blur with 2D Kernel in.... ";
  auto start = std::chrono::steady_clock::now();
  imgboost img2 = blur.perform(img1, 10, 0.5, 5, 0.5);
  auto end = std::chrono::steady_clock::now();
  auto diff = end - start;
  std::cout << std::chrono::duration<double, std::milli>(diff).count() << "ms "
            << std::endl;

  // Make a blur copy using two vectors
  std::cout << "Blur with 2 Vectors.... ";
  start = std::chrono::steady_clock::now();
  imgboost img3 = blur_fast.perform(img1, 10, 0.5, 5, 0.5);
  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << std::chrono::duration<double, std::milli>(diff).count() << "ms "
            << std::endl;

  // Calculate StDev of error between two gaussian algorithms
  imgboost img4 = img3 - img2;
  std::cout << "StDev of error between two approaches: " << img4.Stdev()
            << std::endl;

  img2.savetiff("./2DK blur.tiff");
  img3.savetiff("./2V blur.tiff");

  // Gray-Level Co-occurence matrix with 7 levels, 1 degree resolution and 5
  // pixels displacement
  std::cout << "GLCM.... ";
  start = std::chrono::steady_clock::now();
  ccmatrix ccm(5, 180, 20);
  ccm.perform(img1, 0, 255);
  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << std::chrono::duration<double, std::milli>(diff).count() << "ms "
            << std::endl;
  // Output in csv files, the whole matrix, features and averaged features per
  // displacement
  ccm.printlast();
  ccm.printfeatures();
  ccm.printavgfeatures();
  return 0;
}

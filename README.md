## Synopsis

An image processing library. It eases operations between 2D-arrays and provides
standard image processing algoritms i.e. Gaussian/blur filtering, Co-occurrence matrix
features (energy, entropy, contrast, homogeneity & correlation) ...more to come.
It supports raw data from csv and, optionally, tiff image format.

## Dependences

* <a href="https://github.com/cbontoz/smartix.git" target="_blank">smartix</a>
 * You should include ```smartix.hpp``` in your project.
* <a href="http://www.libtiff.org/" target="_blank">tiffio</a> is optional
 * Install ```sudo apt-get install libtiff5-dev```
 * To remove tiff support remove ```TIFF_EN_`` from preprocessor & -ltiff.

## Using the Library

### Save/Load CSV
The first row of csv provides format information:

width | XXX | height | XXX | bytes/pixel | X
--- | --- | --- | --- | --- | ---


To save: ```yourimgboost.savecsv(path.csv, BytesPerPixel)```, where
```BytesPerPixel``` :
1. grayscale
3. rgb
4. rgba

To load: ```imgboost yourimgboost(path.csv)```, or load in existing instance
```yourimgboost.loadcsv(path.csv)``` .

Save & Load TIFF follows the same concept, but you need to use savetiff and
loadtiff instead.

### Example
```cpp
// Load a tiff
imgboost img1("./sample-data/butterfly.tif");
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

// Make a blur copy using 2D gaussian filter kernel
std::cout << "Blur with 2D Kernel in.... ";
gaussian blur;
auto start = std::chrono::steady_clock::now();
imgboost img2 = blur.perform(img1, 50, 0.5, 5, 0.5, true, false);
auto end = std::chrono::steady_clock::now();
auto diff = end - start;
std::cout << std::chrono::duration<double, std::milli>(diff).count() << "ms "
          << std::endl;

// Make a blur copy using two vectors
std::cout << "Blur with 2 Vector.... ";
start = std::chrono::steady_clock::now();
imgboost img3 = blur.perform(img1, 50, 0.5, 5, 0.5, true, true);
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
ccmatrix ccm(7, 180, 5);
ccm.perform(img3, 0, 256);
end = std::chrono::steady_clock::now();
diff = end - start;
std::cout << std::chrono::duration<double, std::milli>(diff).count() << "ms "
          << std::endl;
// Output in csv files, the whole matrix, features and averaged features per
// displacement
ccm.printlast();
ccm.printfeatures();
ccm.printavgfeatures();
```
Expected output:

Blur with 2D Kernel in.... 2236.33ms

Blur with 2 Vectors.... 215.8ms

StDev of error between two approaches: 1.03603

GLCM.... 1410.83ms

## License

MIT License

Copyright (c) 2017-2018 Christos Bontozoglou

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

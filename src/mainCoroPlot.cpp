/* 
 * CoroEval
 *
 * An evaluation tool for coronary artery reconstructions.
 *
 * Copyright © 2014:
 *
 * Christoph Forman, Universität Erlangen-Nürnberg
 * christoph.forman@cs.fau.de
 *
 * Chris Schwemmer, Universität Erlangen-Nürnberg
 * chris.schwemmer@cs.fau.de
 *
 * Jens Wetzl, Universität Erlangen-Nürnberg
 * jens.wetzl@cs.fau.de
 *
 * CoroEval is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CoroEval is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CoroEval.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "bsplineinterpolation.h"
#include "Data.h"

#include "measurementpoint.h"
#include "vesselsharpness.h"

#include "normalplot.h"
#include "measurementpoint.h"

#include "evaldialog.h"

#include "Wm3Vector3.h"

#include "CImg.h"

#include <QObject>
#include <QApplication>
#include <QWidget>
#include <QFile>

#include <iostream>
#include <fstream>


Data* loadRawFile(char* filename, size_t sizeX, size_t sizeY, size_t sizeZ, double sizeVx)
{
	if( sizeX < 1 || sizeY < 1 || sizeZ < 1 || sizeVx < 0.0)
        return 0;

	Data* data = new Data();
    data->create(Data::COL, sizeX,
                 Data::LIN, sizeY,
                 Data::PAR, sizeZ);
    data->setData(new float[data->size()]);
	data->setVoxelSize(sizeVx);
	std::ifstream in;
	in.open( filename, std::ios::binary );
    if( !in )
    {
		std::cerr << "File does not exists:" << filename << std::endl;
        delete [] data->data0();
        delete data;
        return 0;
    }

    size_t size = data->size();
    unsigned short* d = new unsigned short[size];
    in.read((char*)d, size*sizeof(unsigned short));
	in.close();

	if (!in.good() || in.fail() || (in.gcount() != size * sizeof(unsigned short))) {
		std::cerr << "Error reading file: Could not read from data file" << std::endl;
		delete [] data->data0();
		delete data;
		return 0;
	}

    float* data_data = data->data0();
    for(size_t i=0; i<size; i++)
		data_data[i] = d[i];
    delete [] d;

	return data;
}

std::vector< Wm3::Vector3<double> > loadControlPoints(char* filename)
{
	std::vector< Wm3::Vector3<double> > result;

	// Parse spline file
    std::ifstream file( filename );
    if (!file.good()) {
        std::cerr << "Could not open spline file" << std::endl;
        return result;
    }

    Wm3::Vector3<double> pt;
    std::stringstream line;
    std::string tmp;

    do {
        if (!std::getline(file, tmp))
            break;

        line.clear();
        line << tmp.c_str();

        line >> pt[0] >> pt[1] >> pt[2];
        if (line.bad())
            break;

        result.push_back(pt);

    } while (!file.eof() && file.good());

    file.close();

	return result;
}

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	if( argc < 7 )
	{
		std::cerr << "Usage: " << argv[0] << 
			" [rawFile] [sizeX] [sizeY] [sizeZ] [sizeVx] [splineFile]" << std::endl;
		return EXIT_FAILURE;
	}

	char*   rawFile		= argv[1];
	size_t	sizeX		= atoi(argv[2]);
    size_t	sizeY		= atoi(argv[3]);
    size_t	sizeZ		= atoi(argv[4]);
	double	sizeVx		= atof(argv[5]);
	char*   splineFile	= argv[6];
	char*   imageFile	= argv[7];

	Data* data = loadRawFile( rawFile, sizeX, sizeY, sizeZ, sizeVx );
	if( data == 0 )
		return EXIT_FAILURE;

	//std::cout << "Loaded " << argv[1] << "." << std::endl;
/*
	std::cout << "Loaded RawFile" << argv[1] << "." << std::endl;
	std::cout << "SizeX: " << data->getLen(Data::COL) << std::endl;
	std::cout << "SizeY: " << data->getLen(Data::LIN) << std::endl;
	std::cout << "SizeZ: " << data->getLen(Data::PAR) << std::endl;
	std::cout << "Voxel: " << data->getVoxelSize() << std::endl;
*/
	std::vector< Wm3::Vector3<double> > ctr = loadControlPoints( splineFile );
	if( ctr.empty() )
		return EXIT_FAILURE;

//	std::cout << "Loaded " << ctr.size() << " ControlPoints " << std::endl;

	BSplineInterpolation bSpline;
	bSpline.setSamples(ctr);

	double length = bSpline.getLength() * data->getVoxelSize();
	double delta = 0.5 / length;

	std::cout << "Total length: " << length << std::endl;
	
	// Compute the CenterOfMass of all ControlPoints
	Wm3::Vector3<double> centerOfmass = Wm3::Vector3<double>::ZERO;
	for(size_t i=0; i<ctr.size(); i++)
		centerOfmass += ctr[i];
	centerOfmass /= ctr.size();

	int sizeNormal = 42;
    double hwidth  = double( sizeNormal ) / 2.0;
	int width  = sizeNormal*2 + 1;
	int height = int(1.0 / delta + 1.5);

	//std::cout << width << " x " << height << std::endl;

	cimg_library::CImg<float> cImage( width, height );
	float* image = cImage.data();

    Wm3::Vector3<double> x   = centerOfmass - bSpline.pointAt( 0.5 );
	Wm3::Vector3<double> y   = bSpline.pointAt( 0.0 ) - bSpline.pointAt( 1.0 );
    Wm3::Vector3<double> z   = x.Cross(y);
    z.Normalize();

/*
	std::cout	<< "X: " << x << std::endl 
				<< "Y: " << y << std::endl 
				<< "Z: " << z << std::endl << std::endl;
*/
	std::vector<MeasurementPoint*> mpts;
	int pos = 0;
	for(double i=0.0; i<=1.0; i += delta)
    {
		MeasurementPoint::Normal n;
		//std::cout << i << " ";
        Wm3::Vector3<double> xi  = bSpline.pointAt( i );
		//std::cout << xi << std::endl;

        // Sample the normal (2x oversampled)
		for (double j = -hwidth; j <= hwidth; j += 0.5)
        {
            Wm3::Vector3<double> currPt = xi + z * j;
			float value = data->getValue( currPt.X(), currPt.Y(), currPt.Z() );
            image[pos++] = value;
			n.data.push_back( value );
        }
		VesselSharpness eval;
		eval.setNormal( n.data, data->getVoxelSize() );
		n.sharpness = eval.getSharpness(Interval::Both);
		n.diffMaxToCenter = eval.getCenterDistance()*0.5;
		n.centerInterval = eval.getCenterInterval();
		MeasurementPoint* mp = new MeasurementPoint(n, i*length, data);
		mpts.push_back( mp );
    }

	EvalDialog* d = new EvalDialog(imageFile, mpts, cImage );
	
	QFile file(":/res/style.qss");
    file.open(QFile::ReadOnly);
    QString styleSheet = QLatin1String(file.readAll());
    d->setStyleSheet(styleSheet);
	
	d->exec();

    // Draw result as bmp for debug purposes:
    //cImage.normalize(0, 255);
    //cImage.save( imageFile );

	// Clean up
	delete [] data->data0();
	delete data;

    return EXIT_SUCCESS;
}

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

#include "loaddialog.h"
#include "ui_loaddialog.h"

#include "ZipFile.h"
#include "Data.h"
#include "settings.h"

#include <QDoubleValidator>
#include <QIntValidator>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>

#include <gdcmImageReader.h>
#include <gdcmImage.h>
#include <gdcmAttribute.h>
#include <gdcmDirectory.h>

#include <limits>
#include <strstream>


#ifdef max
#undef max
#endif


LoadDialog::LoadDialog(QWidget* parent) :
	m_settings(Settings::getSettings()),
    QDialog(parent),
    ui(new Ui::LoadDialog)
{
    ui->setupUi(this);

    QIntValidator* val = new QIntValidator(1, std::numeric_limits<int>::max(), 0);
    ui->le_sizeX->setValidator(val);
    ui->le_sizeY->setValidator(val);
    ui->le_sizeZ->setValidator(val);

    QDoubleValidator* val2 = new QDoubleValidator(0.0, std::numeric_limits<double>::max(), 3, 0);
    ui->le_pxSize->setValidator(val2);

	ui->tabWidget->setCurrentIndex(1);

	ui->le_sizeX->setText(m_settings.getRawDefaultSizeX());
    ui->le_sizeY->setText(m_settings.getRawDefaultSizeY());
    ui->le_sizeZ->setText(m_settings.getRawDefaultSizeZ());
    ui->le_pxSize->setText(m_settings.getRawDefaultSizeVox());

    connect(ui->b_browseRaw, SIGNAL(clicked()), this, SLOT(browse_clicked()) );
    connect(ui->buttonBox, SIGNAL(accepted()), this, SLOT(loadFile()));
	connect(ui->b_browseCavarev, SIGNAL(clicked()), this, SLOT(browse_clicked()));
	connect(ui->b_browseDicom, SIGNAL(clicked()), this, SLOT(browse_clicked()));
	connect(ui->cb_loadFromZIP, SIGNAL(stateChanged(int)), this, SLOT(cb_changed(int)));

	m_lastLoadDir = m_settings.getDefaultDir();
}

LoadDialog::~LoadDialog()
{
    delete ui;
}

void LoadDialog::setLastLoadDir(const QString& dir)
{
	m_lastLoadDir = dir;
}

QString LoadDialog::getLastLoadDir() const
{
	return m_lastLoadDir;
}

void LoadDialog::browse_clicked()
{
	if (ui->tabWidget->currentWidget() == ui->w_raw) {
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Load Raw File"), m_lastLoadDir,
			tr("All Files (*)"));

		ui->le_filenameRaw->setText(fileName);

		if (m_settings.getRememberLastDir() && !fileName.isEmpty()) {
			QFileInfo info(fileName);
			m_lastLoadDir = info.path();
		}
	} else if (ui->tabWidget->currentWidget() == ui->w_dicom) {
		if (ui->cb_loadFromZIP->isChecked()) {
			QString fileName = QFileDialog::getOpenFileName(this,
				tr("Load DICOM Archive"), m_lastLoadDir,
				tr("ZIP Files (*.zip)"));

			ui->le_filenameDicom->setText(fileName);

			if (m_settings.getRememberLastDir() && !fileName.isEmpty()) {
				QFileInfo info(fileName);
				m_lastLoadDir = info.path();
			}
		} else {
			QString dirName = QFileDialog::getExistingDirectory(this,
				tr("Load DICOM files from"), m_lastLoadDir);
			
			ui->le_filenameDicom->setText(dirName);

			if (m_settings.getRememberLastDir() && !dirName.isEmpty())
				m_lastLoadDir = dirName;
		}
	} else if (ui->tabWidget->currentWidget() == ui->w_cavarev) {
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Load Cavarev File"), m_lastLoadDir,
			tr("Binary Files (*.bin)"));

		ui->le_filenameCavarev->setText(fileName);

		if (m_settings.getRememberLastDir() && !fileName.isEmpty()) {
			QFileInfo info(fileName);
			m_lastLoadDir = info.path();
		}
	} else {
		QMessageBox::critical(NULL, "Error", "Something is very wrong here");
		return;
	}
}

void LoadDialog::loadFile()
{
	if (ui->tabWidget->currentWidget() == ui->w_raw)
		loadRaw();
	else if (ui->tabWidget->currentWidget() == ui->w_dicom) {
		if (ui->cb_loadFromZIP->isChecked())
			loadDicomZIP();
		else
			loadDicomDir();
	} else if (ui->tabWidget->currentWidget() == ui->w_cavarev)
		loadCavarev();
	else
		QMessageBox::critical(NULL, "Error", "Something is very wrong here");  
}

void LoadDialog::loadRaw() {
	QString filename = ui->le_filenameRaw->text();
	size_t sizeX = ui->le_sizeX->text().toUInt();
    size_t sizeY = ui->le_sizeY->text().toUInt();
    size_t sizeZ = ui->le_sizeZ->text().toUInt();
    double sizePx = ui->le_pxSize->text().toDouble();

    if(filename == "")
        return;

    if( sizeX < 1 || sizeY < 1 || sizeZ < 1 || sizePx < 0.0)
        return;

	emit clearData();

    Data* data = new Data();
    data->create(Data::COL, sizeX,
                 Data::LIN, sizeY,
                 Data::PAR, sizeZ);
    data->setData(new float[data->size()]);

    std::ifstream in;
	QByteArray fn = filename.toLocal8Bit();
	in.open(fn.constData(), std::ios::binary );
    if( !in )
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Error");
        msgBox.setText( QString("File does not exists: %1").arg(filename) );
        msgBox.exec();

        delete [] data->data0();
        delete data;

        return;
    }


    size_t size = data->size();
    if( ui->comboBox->currentIndex() == 0 )
    {
        // is Short
        unsigned short* d = new unsigned short[size];
        in.read((char*)d, size*sizeof(unsigned short));
				
		if (!in.good() || in.fail() || (in.gcount() != size * sizeof(unsigned short))) {
			QMessageBox::critical(NULL, "Error reading file", "Could not read from data file");

			in.close();
			delete [] data->data0();
			delete data;

			return;
		}

        float* data_data = data->data0();
        for(size_t i=0; i<size; i++)
            data_data[i] = d[i];
        delete [] d;

    } else if(ui->comboBox->currentIndex() == 1) {
        // is Float
        in.read((char*)data->data0(), size*sizeof(float));
				
		if (!in.good() || in.fail() || (in.gcount() != size * sizeof(float))) {
			QMessageBox::critical(NULL, "Error reading file", "Could not read from data file");

			in.close();
			delete [] data->data0();
			delete data;

			return;
		}
    }
	
	in.close();

    data->setVoxelSize(sizePx);

    emit sendData(data);
    
	QStringList filenameList = filename.split("/");
    QString plainFilename = filenameList.back();
    QString proband = filenameList.at( filenameList.size()-2 );
	QString title("CoroEval: ");
	title.append(proband);
	title.append("/");
	title.append(plainFilename);
    
	emit sendFileName(title);
}

void LoadDialog::loadCavarev() {
	QString filename = ui->le_filenameCavarev->text();

	if (filename == "")
		return;

    std::ifstream in;
	QByteArray fn = filename.toLocal8Bit();
	in.open(fn.constData(), std::ios::binary);
    if (!in) {
		QMessageBox::critical(NULL, "Error", QString("File does not exists: %1").arg(filename));
        
        return;
    }

	// Read volume information
	float o[3]; // Origin (unused)
	in.read(reinterpret_cast<char*>(o), 3 * sizeof(float));

	unsigned int size[3];
	in.read(reinterpret_cast<char*>(size), 3 * sizeof(unsigned int));

	float sizePx;
	in.read(reinterpret_cast<char*>(&sizePx), sizeof(float));

	if (!in.good() || in.fail()) {
		QMessageBox::critical(NULL, "Error reading file", "Could not read from data file");

		return;
	}

	if (size[0] < 1 || size[0] < 1 || size[0] < 1 || sizePx < 0.0) {
		QMessageBox::critical(NULL, "Error reading file", "Invalid volume information");

		return;
	}

	emit clearData();

	Data* data = new Data();
	data->create(Data::COL, size[0],
				 Data::LIN, size[1],
				 Data::PAR, size[2]);
	data->setData(new float[data->size()]);
	
	size_t dataSize = data->size();
    
	// cavarev is unsigned char data
	unsigned char* d = new unsigned char[dataSize];
	in.read(reinterpret_cast<char*>(d), dataSize);

	if (!in.good() || in.fail() || (in.gcount() != dataSize)) {
		QMessageBox::critical(NULL, "Error reading file", "Could not read from data file");

		delete [] data->data0();
		delete data;

		return;
	}

	float* data_data = data->data0();
	for (size_t i = 0; i < dataSize; i++)
		data_data[i] = d[i];
	delete[] d;

	data->setVoxelSize(sizePx);

	emit sendData(data);
	QStringList filenameList = filename.split("/");
	QString plainFilename = filenameList.back();
	QString proband = filenameList.at( filenameList.size()-2 );
	emit sendFileName(QString("CoroEval: ").append(proband).append("/").append(plainFilename) );
}

void LoadDialog::cb_changed(int state) {
	if (state == Qt::Unchecked) {
		ui->l_filenameDicom->setText(tr("Directory name:"));
	} else {
		ui->l_filenameDicom->setText(tr("Filename:"));
	}
}

void LoadDialog::loadDicomDir() {
	QString dirname = ui->le_filenameDicom->text();

	if (dirname == "")
		return;

	Data* data = new Data();

	gdcm::Directory d;
	d.Load(dirname.toStdString());
	gdcm::Directory::FilenamesType dirList = d.GetFilenames();
	size_t nSlices = dirList.size();

	unsigned int ndim;
	unsigned int dims[2];
	double voxSize;
	size_t sliceOffset;

	// Load first slice
	gdcm::ImageReader reader;
	reader.SetFileName(dirList[0].c_str());
	if (!reader.Read()) {
		QString msg;
		msg.append("Loading from DICOM dir failed: ");
		msg.append("Cannot read slice");
		QMessageBox::critical(NULL, "Error", msg);
		delete data;
		return;
	}

	// Get and check slice information
	const gdcm::Image& image = reader.GetImage();
	ndim = image.GetNumberOfDimensions();
			
	const unsigned int* dimsP = image.GetDimensions();
	dims[0] = dimsP[0];
	dims[1] = dimsP[1];

	const double* voxSizeP = image.GetSpacing();
	voxSize = voxSizeP[0];

	gdcm::PixelFormat pixeltype = image.GetPixelFormat();
			 
	const gdcm::DataSet& ds = reader.GetFile().GetDataSet();
	gdcm::Attribute<0x0018, 0x0050> voxSizeSlicesAttr;
	//voxSizeSlicesAttr.SetFromDataSet(ds);
	voxSizeSlicesAttr.SetFromDataElement(ds.GetDataElement(voxSizeSlicesAttr.GetTag()));
	const gdcm::Attribute<0x0018, 0x0050>::ArrayType& voxSizeSlices = voxSizeSlicesAttr.GetValue();

	if (voxSizeP[0] != voxSizeSlices)
	{
		QString msg("Volume resolution is not isotropic, measurements reported by CoroEval may be wrong! In-plane resolution (%1 mm) will override slice resolution (%2 mm).");
		QMessageBox::warning(NULL, "Warning", msg.arg(QString::number(voxSizeP[0], 'g', 2)).arg(QString::number(voxSizeSlices, 'g', 2)));
	}

	if ((ndim != 2) || (voxSizeP[0] != voxSizeP[1]) || /*(voxSizeP[0] != voxSizeSlices) ||*/ (pixeltype.GetBitsAllocated() > 32)) {
		QString msg;
		msg.append("Loading from DICOM dir failed: ");
		msg.append("Unsupported DICOM image");
		QMessageBox::critical(NULL, "Error", msg);
		delete data;
		return;
	}

	emit clearData();

	// Allocate data
	data->create(Data::COL, dims[0],
					Data::LIN, dims[1],
					Data::PAR, nSlices);
	data->setData(new float[data->size()]);
	data->setVoxelSize(voxSize);

	size_t dataSize = data->size();
	sliceOffset = dims[0] * dims[1];
			
	// Store the first slice
	if (!loadDicomSlice(pixeltype, image, sliceOffset, 0, data->data0())) {
		QString msg;
		msg.append("Loading from DICOM dir failed: ");
		msg.append("Cannot store slice");
		QMessageBox::critical(NULL, "Error", msg);
		delete[] data->data0();
		delete data;
	}

	// Load remaining slices
	for (size_t slice = 1; slice < nSlices; ++slice) {
		gdcm::ImageReader reader;
		reader.SetFileName(dirList[slice].c_str());
		if (!reader.Read()) {
			QString msg;
			msg.append("Loading from DICOM dir failed: ");
			msg.append("Cannot read slice");
			QMessageBox::critical(NULL, "Error", msg);
			delete[] data->data0();
			delete data;
			return;
		}

		const gdcm::Image& image = reader.GetImage();
		unsigned int ndim2 = image.GetNumberOfDimensions();
		const unsigned int* dims2 = image.GetDimensions();
		const double* voxSize2 = image.GetSpacing();
		gdcm::PixelFormat pixeltype = image.GetPixelFormat();

		if ((ndim != 2) || (voxSizeP[0] != voxSizeP[1]) || /*(voxSizeP[0] != voxSizeSlices) ||*/ (pixeltype.GetBitsAllocated() > 32)) {
			QString msg;
			msg.append("Loading from DICOM dir failed: ");
			msg.append("Unsupported DICOM image");
			QMessageBox::critical(NULL, "Error", msg);
			delete[] data->data0();
			delete data;
			return;
		}

		if (!loadDicomSlice(pixeltype, image, sliceOffset, slice, data->data0())) {
			QString msg;
			msg.append("Loading from DICOM dir failed: ");
			msg.append("Cannot store slice");
			QMessageBox::critical(NULL, "Error", msg);
			delete[] data->data0();
			delete data;
			return;
		}
	}

	emit sendData(data);
	emit sendFileName(QString("CoroEval: ").append(dirname));
}

void LoadDialog::loadDicomZIP() {
	QString filename = ui->le_filenameDicom->text();

	if (filename == "")
		return;

	Data* data = new Data();

	try {
		// Open ZIP
		ZipFile zip(filename.toStdString());
		if (zip.nEntries() < 1)
			throw std::runtime_error("ZIP file does not contain any slices");

		char* tmpBuf;
		zip_uint64_t bufLen;

		unsigned int ndim;
		unsigned int dims[2];
		double voxSize;
		size_t sliceOffset;
		
		// Load first slice
		zip.loadFile(0, tmpBuf, bufLen);
		
		{
			std::istrstream tmpStream(const_cast<const char*>(tmpBuf), static_cast<std::streamsize>(bufLen));
			gdcm::ImageReader reader;
			reader.SetStream(tmpStream);
			if (!reader.Read()) {
				delete[] tmpBuf;
				throw std::runtime_error("Error reading DICOM slice");
			}
			
			// Get and check slice information
			const gdcm::Image& image = reader.GetImage();
			ndim = image.GetNumberOfDimensions();
			
			const unsigned int* dimsP = image.GetDimensions();
			dims[0] = dimsP[0];
			dims[1] = dimsP[1];

			const double* voxSizeP = image.GetSpacing();
			voxSize = voxSizeP[0];

			gdcm::PixelFormat pixeltype = image.GetPixelFormat();
			 
			const gdcm::DataSet& ds = reader.GetFile().GetDataSet();
			gdcm::Attribute<0x0018, 0x0050> voxSizeSlicesAttr;
			//voxSizeSlicesAttr.SetFromDataSet(ds);
			voxSizeSlicesAttr.SetFromDataElement(ds.GetDataElement(voxSizeSlicesAttr.GetTag()));
			const gdcm::Attribute<0x0018, 0x0050>::ArrayType& voxSizeSlices = voxSizeSlicesAttr.GetValue();

			if ((ndim != 2) || (voxSizeP[0] != voxSizeP[1]) || (voxSizeP[0] != voxSizeSlices) || (pixeltype.GetBitsAllocated() > 32)) {
				delete[] tmpBuf;
				throw std::runtime_error("Unsupported DICOM image");
			}

			emit clearData();

			// Allocate data
			data->create(Data::COL, dims[0],
						 Data::LIN, dims[1],
						 Data::PAR, zip.nEntries());
			data->setData(new float[data->size()]);
			data->setVoxelSize(voxSize);

			size_t dataSize = data->size();
			sliceOffset = dims[0] * dims[1];
			
			// Store the first slice
			if (!loadDicomSlice(pixeltype, image, sliceOffset, 0, data->data0())) {
				delete[] tmpBuf;

				throw std::runtime_error("Cannot read slice");
			}
		}

		// Load remaining slices
		for (zip_uint64_t i = 1; i < zip.nEntries(); ++i) {
			// Free memory from previous slice
			delete[] tmpBuf;

			zip.loadFile(i, tmpBuf, bufLen);

			std::istrstream tmpStream(const_cast<const char*>(tmpBuf), static_cast<std::streamsize>(bufLen));
			gdcm::ImageReader reader;
			reader.SetStream(tmpStream);
			if (!reader.Read()) {
				delete[] tmpBuf;
				throw std::runtime_error("Error reading DICOM slice");
			}

			const gdcm::Image& image = reader.GetImage();
			unsigned int ndim2 = image.GetNumberOfDimensions();
			const unsigned int* dims2 = image.GetDimensions();
			const double* voxSize2 = image.GetSpacing();
			gdcm::PixelFormat pixeltype = image.GetPixelFormat();

			if ((ndim2 != 2) || (voxSize2[0] != voxSize) || (voxSize2[1] != voxSize) || (pixeltype.GetBitsAllocated() > 32)) {
				delete[] tmpBuf;
				throw std::runtime_error("Unsupported DICOM image");
			}

			if (!loadDicomSlice(pixeltype, image, sliceOffset, i, data->data0())) {
				delete[] tmpBuf;

				throw std::runtime_error("Cannot read slice");
			}
		}

		// Free memory
		delete[] tmpBuf;
	} catch (std::exception& e) {
		QString msg;
		msg.append("Loading from ZIP file failed: ");
		msg.append(e.what());
		QMessageBox::critical(NULL, "Error", msg);
		if (data->data0())
			delete[] data->data0();
		delete data;
		return;
	}

	emit sendData(data);
	
	QStringList filenameList = filename.split("/");
	QString plainFilename = filenameList.back();
	QString proband = filenameList.at(filenameList.size() - 2);
	QString name("CoroEval: ");
	name.append(proband);
	name.append("/");
	name.append(plainFilename);
	
	emit sendFileName(name);
}

bool LoadDialog::loadDicomSlice(const gdcm::PixelFormat& pixeltype, const gdcm::Image& image, size_t sliceOffset, size_t index, float* data) {
	if (pixeltype == gdcm::PixelFormat::FLOAT32) {
		// Easy case
		if ((image.GetBufferLength() != sizeof(float) * sliceOffset) ||
			!image.GetBuffer(reinterpret_cast<char*>(data + index * sliceOffset))) {
			return false;
		}
	} else {
		// Conversion necessary
		float* slicePtr = data + index * sliceOffset;
		char* tmpImg = new char[image.GetBufferLength()];
		if (!image.GetBuffer(tmpImg)) {
			delete[] tmpImg;
			return false;
		}

		switch (pixeltype) {
		case gdcm::PixelFormat::INT8:
		case gdcm::PixelFormat::UINT8: {
			for (size_t i = 0; i < sliceOffset; ++i)
				slicePtr[i] = tmpImg[i];

			break;
		}
		case gdcm::PixelFormat::INT16: {
			short int* tmpImgP = reinterpret_cast<short int*>(tmpImg);
			for (size_t i = 0; i < sliceOffset; ++i)
				slicePtr[i] = tmpImgP[i];

			break;
		}
		case gdcm::PixelFormat::UINT16: {
			unsigned short int* tmpImgP = reinterpret_cast<unsigned short int*>(tmpImg);
			for (size_t i = 0; i < sliceOffset; ++i)
				slicePtr[i] = tmpImgP[i];

			break;
		}
		case gdcm::PixelFormat::INT32: {
			int* tmpImgP = reinterpret_cast<int*>(tmpImg);
			for (size_t i = 0; i < sliceOffset; ++i)
				slicePtr[i] = tmpImgP[i];

			break;
		}
		case gdcm::PixelFormat::UINT32: {
			unsigned int* tmpImgP = reinterpret_cast<unsigned int*>(tmpImg);
			for (size_t i = 0; i < sliceOffset; ++i)
				slicePtr[i] = tmpImgP[i];

			break;
		}
		default: {
			delete[] tmpImg;
			return false;

			break;
		}
		}

		delete[] tmpImg;
	}

	return true;
}

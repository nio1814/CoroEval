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

#include "settings.h"

#include <QDateTime>
#include <QMessageBox>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


Settings::Settings()
{
    m_numNormalsAtPoint     = 10;
    m_sizeNormal            = 21;
	m_minRelMagnitude       = 0.5;
	m_plateauThreshold      = 15.0;
	m_minVesselDiameter     = 1.0;
    m_sharpnessEvalInterval = 1.0;
	m_doPointCheck          = false;
	m_doNormalCheck			= true;
	m_useEllipse			= true;
	m_showMIP				= false;
	m_showNormals			= false;
	m_defaultDir			= "D:/data/";
	m_rememberDir			= true;
	m_rawDefaultSizeX		= "";
	m_rawDefaultSizeY		= "";
	m_rawDefaultSizeZ		= "";
	m_rawDefaultSizeVox		= "";
}

Settings& Settings::getSettings()
{
    static Settings instance;
    return instance;
}

void Settings::setNumNormalsAtPoint(size_t val)    { m_numNormalsAtPoint = val; }
size_t Settings::getNumNormalsAtPoint() const      { return m_numNormalsAtPoint; }

void Settings::setSizeNormal(size_t val)   { m_sizeNormal = val; }
size_t  Settings::getSizeNormal() const    { return m_sizeNormal; }

void Settings::setMinRelMagnitude(double val) { m_minRelMagnitude = val; }
double Settings::getMinRelMagniude() const    { return m_minRelMagnitude; }

void Settings::setPlateauThreshold(double val) { m_plateauThreshold = val; }
double Settings::getPlateauThreshold() const   { return m_plateauThreshold; }

void Settings::setMinVesselDiameter(double val) { m_minVesselDiameter = val; }
double Settings::getMinVesselDiameter() const   { return m_minVesselDiameter; }

void Settings::setSharpnessEvalInterval(double val) { m_sharpnessEvalInterval = val; }
double  Settings::getSharpnessEvalInterval() const  { return m_sharpnessEvalInterval; }

void Settings::setDoPointCheck(bool val) { m_doPointCheck = val; }
bool Settings::getDoPointCheck() const { return m_doPointCheck; }

void Settings::setDoNormalCheck(bool val) { m_doNormalCheck = val; }
bool Settings::getDoNormalCheck() const { return m_doNormalCheck; }

void Settings::setUseEllipse(bool val) { m_useEllipse = val; }
bool Settings::getUseEllipse() const { return m_useEllipse; }

void Settings::setShowMIP(bool val) { m_showMIP = val; }
bool Settings::getShowMIP() const { return m_showMIP; }

void Settings::setShowNormals(bool val) { m_showNormals = val; }
bool Settings::getShowNormals() const { return m_showNormals; }

void Settings::setDefaultDir(const QString& val) { m_defaultDir = val; }
const QString& Settings::getDefaultDir() const { return m_defaultDir; }

void Settings::setRememberLastDir(bool val) { m_rememberDir = val; }
bool Settings::getRememberLastDir() const { return m_rememberDir; }

const QString& Settings::getRawDefaultSizeX() const { return m_rawDefaultSizeX; }
const QString& Settings::getRawDefaultSizeY() const { return m_rawDefaultSizeY; }
const QString& Settings::getRawDefaultSizeZ() const { return m_rawDefaultSizeZ; }
const QString& Settings::getRawDefaultSizeVox() const { return m_rawDefaultSizeVox; }


static const QString numNormalsAtPointName("NumProfiles");
static const QString sizeNormalName("LengthProfile");
static const QString minRelMagnitudeName("MinRelMag");
static const QString plateauThresholdName("PlateauThreshold");
static const QString minVesselDiameterName("MinVesselDiameter");
static const QString doPointCheckName("Check20Within");
static const QString doNormalCheckName("CheckOutliers");
static const QString useEllipseName("EllipseFit");
static const QString sharpnessEvalIntervalName("SplineSampling");
static const QString showMIPName("ShowMIP");
static const QString showNormalsName("ShowNormals");
static const QString defaultDirName("DefaultDir");
static const QString rememberDirName("RememberLastDir");
static const QString rawDefaultSizeX("RawDefaultSizeX");
static const QString rawDefaultSizeY("RawDefaultSizeY");
static const QString rawDefaultSizeZ("RawDefaultSizeZ");
static const QString rawDefaultSizeVox("RawDefaultSizeVox");


bool Settings::load(const QString& filename) {
	QByteArray fn = filename.toLocal8Bit();
	std::ifstream config(fn.constData());
	if (!config.is_open() || !config.good())
		return false;

	while (!config.eof() && config.good()) {
			std::string line;
			if (!std::getline(config, line))
				break;

			if (line.empty())
				continue;

			if ((line[0] == '#') || (line[0] == '/'))
				continue;

			QString qLine = QString::fromStdString(line.c_str());
			qLine = qLine.trimmed();

			if (qLine.startsWith(numNormalsAtPointName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - numNormalsAtPointName.size()).trimmed();
				bool ok; 
				unsigned int i = param.toUInt(&ok);
				if (ok)
					m_numNormalsAtPoint = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(numNormalsAtPointName));
			} else if (qLine.startsWith(sizeNormalName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - sizeNormalName.size()).trimmed();
				bool ok; 
				unsigned int i = param.toUInt(&ok);
				if (ok)
					m_sizeNormal = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(sizeNormalName));
			} else if (qLine.startsWith(minRelMagnitudeName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - minRelMagnitudeName.size()).trimmed();
				bool ok; 
				double i = param.toDouble(&ok);
				if (ok)
					m_minRelMagnitude = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(minRelMagnitudeName));
			} else if (qLine.startsWith(plateauThresholdName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - plateauThresholdName.size()).trimmed();
				bool ok; 
				double i = param.toDouble(&ok);
				if (ok)
					m_plateauThreshold = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(plateauThresholdName));
			} else if (qLine.startsWith(minVesselDiameterName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - minVesselDiameterName.size()).trimmed();
				bool ok; 
				double i = param.toDouble(&ok);
				if (ok)
					m_minVesselDiameter = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(minVesselDiameterName));
			} else if (qLine.startsWith(doPointCheckName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - doPointCheckName.size()).trimmed();
        bool ok = false;
				bool i;

				if ((param.compare("true", Qt::CaseInsensitive) == 0) || (param.compare("1", Qt::CaseInsensitive) == 0) || (param.compare("yes", Qt::CaseInsensitive) == 0)) {
					i = true;
					ok = true;
				} else if ((param.compare("false", Qt::CaseInsensitive) == 0) || (param.compare("0", Qt::CaseInsensitive) == 0) || (param.compare("no", Qt::CaseInsensitive) == 0)) {
					i = false;
					ok = true;
				}

				if (ok)
					m_doPointCheck = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(doPointCheckName));
			} else if (qLine.startsWith(doNormalCheckName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - doNormalCheckName.size()).trimmed();
        bool ok = false;
				bool i;

				if ((param.compare("true", Qt::CaseInsensitive) == 0) || (param.compare("1", Qt::CaseInsensitive) == 0) || (param.compare("yes", Qt::CaseInsensitive) == 0)) {
					i = true;
					ok = true;
				} else if ((param.compare("false", Qt::CaseInsensitive) == 0) || (param.compare("0", Qt::CaseInsensitive) == 0) || (param.compare("no", Qt::CaseInsensitive) == 0)) {
					i = false;
					ok = true;
				}

				if (ok)
					m_doNormalCheck = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(doNormalCheckName));
			} else if (qLine.startsWith(useEllipseName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - useEllipseName.size()).trimmed();
        bool ok = false;
				bool i;

				if ((param.compare("true", Qt::CaseInsensitive) == 0) || (param.compare("1", Qt::CaseInsensitive) == 0) || (param.compare("yes", Qt::CaseInsensitive) == 0)) {
					i = true;
					ok = true;
				} else if ((param.compare("false", Qt::CaseInsensitive) == 0) || (param.compare("0", Qt::CaseInsensitive) == 0) || (param.compare("no", Qt::CaseInsensitive) == 0)) {
					i = false;
					ok = true;
				}

				if (ok)
					m_useEllipse = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(useEllipseName));
			} else if (qLine.startsWith(sharpnessEvalIntervalName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - sharpnessEvalIntervalName.size()).trimmed();
				bool ok; 
				double i = param.toDouble(&ok);
				if (ok)
					m_sharpnessEvalInterval = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(sharpnessEvalIntervalName));
			} else if (qLine.startsWith(showMIPName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - showMIPName.size()).trimmed();
        bool ok = false;
				bool i;

				if ((param.compare("true", Qt::CaseInsensitive) == 0) || (param.compare("1", Qt::CaseInsensitive) == 0) || (param.compare("yes", Qt::CaseInsensitive) == 0)) {
					i = true;
					ok = true;
				} else if ((param.compare("false", Qt::CaseInsensitive) == 0) || (param.compare("0", Qt::CaseInsensitive) == 0) || (param.compare("no", Qt::CaseInsensitive) == 0)) {
					i = false;
					ok = true;
				}

				if (ok)
					m_showMIP = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(showMIPName));
			} else if (qLine.startsWith(showNormalsName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - showNormalsName.size()).trimmed();
        bool ok = false;
				bool i;

				if ((param.compare("true", Qt::CaseInsensitive) == 0) || (param.compare("1", Qt::CaseInsensitive) == 0) || (param.compare("yes", Qt::CaseInsensitive) == 0)) {
					i = true;
					ok = true;
				} else if ((param.compare("false", Qt::CaseInsensitive) == 0) || (param.compare("0", Qt::CaseInsensitive) == 0) || (param.compare("no", Qt::CaseInsensitive) == 0)) {
					i = false;
					ok = true;
        }

				if (ok)
					m_showNormals = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(showNormalsName));
			} else if (qLine.startsWith(defaultDirName, Qt::CaseInsensitive)) {

			} else if (qLine.startsWith(rememberDirName, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - rememberDirName.size()).trimmed();
        bool ok = false;
				bool i;

				if ((param.compare("true", Qt::CaseInsensitive) == 0) || (param.compare("1", Qt::CaseInsensitive) == 0) || (param.compare("yes", Qt::CaseInsensitive) == 0)) {
					i = true;
					ok = true;
				} else if ((param.compare("false", Qt::CaseInsensitive) == 0) || (param.compare("0", Qt::CaseInsensitive) == 0) || (param.compare("no", Qt::CaseInsensitive) == 0)) {
					i = false;
					ok = true;
				}

				if (ok)
					m_rememberDir = i;
				else
					QMessageBox::warning(NULL, "Warning", param.append(" is not a valid value for ").append(rememberDirName));
			} else if (qLine.startsWith(rawDefaultSizeX, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - rawDefaultSizeX.size()).trimmed();
				if (!param.isEmpty())
					m_rawDefaultSizeX = param;
			} else if (qLine.startsWith(rawDefaultSizeY, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - rawDefaultSizeY.size()).trimmed();
				if (!param.isEmpty())
					m_rawDefaultSizeY = param;
			} else if (qLine.startsWith(rawDefaultSizeZ, Qt::CaseInsensitive)) {
				QString param = qLine.right(qLine.size() - rawDefaultSizeZ.size()).trimmed();
				if (!param.isEmpty())
					m_rawDefaultSizeZ = param;
			} else if (qLine.startsWith(rawDefaultSizeVox, Qt::CaseInsensitive)) {
			QString param = qLine.right(qLine.size() - rawDefaultSizeVox.size()).trimmed();
				if (!param.isEmpty())
					m_rawDefaultSizeVox = param;
			} else QMessageBox::warning(NULL, "Warning", QString("Encountered unknown option ").append(qLine));
	}

	config.close();

	return true;
}

bool Settings::save(const QString& filename) {
	QByteArray fn = filename.toLocal8Bit();
	std::ofstream config(fn.constData());
	if (!config.is_open() || !config.good())
		return false;

	QString dt(QDateTime::currentDateTime().toString("ddd MMMM d yyyy, hh:mm:ss"));
	config << "# CoroEval config file written " << dt.toStdString() << std::endl;

	config << numNormalsAtPointName.toStdString() << "\t\t\t" << m_numNormalsAtPoint << std::endl;
	config << sizeNormalName.toStdString() << "\t\t" << m_sizeNormal << std::endl;
	config << minRelMagnitudeName.toStdString() << "\t\t\t"   << m_minRelMagnitude << std::endl;
	config << plateauThresholdName.toStdString() << "\t"<< m_plateauThreshold << std::endl;
	config << minVesselDiameterName.toStdString() << "\t" << m_minVesselDiameter << std::endl;
	config << doPointCheckName.toStdString() << "\t\t"   << ((m_doPointCheck) ? "true" : "false") << std::endl;
	config << doNormalCheckName.toStdString() << "\t\t"   << ((m_doNormalCheck) ? "true" : "false") << std::endl;
	config << useEllipseName.toStdString() << "\t\t\t" << ((m_useEllipse) ? "true" : "false") << std::endl;
	config << sharpnessEvalIntervalName.toStdString() << "\t\t"  << m_sharpnessEvalInterval << std::endl;
	config << showMIPName.toStdString() << "\t\t\t\t"     << ((m_showMIP) ? "true" : "false") << std::endl;
	config << showNormalsName.toStdString() << "\t\t\t"   << ((m_showNormals) ? "true" : "false") << std::endl;
	config << defaultDirName.toStdString() << "\t\t\t"    << m_defaultDir.toStdString() << std::endl;
	config << rememberDirName.toStdString() << "\t\t" << ((m_rememberDir) ? "true" : "false") << std::endl;

	if (!m_rawDefaultSizeX.isEmpty())
		config << rawDefaultSizeX.toStdString() << "\t\t" << m_rawDefaultSizeX.toStdString() << std::endl;

	if (!m_rawDefaultSizeY.isEmpty())
		config << rawDefaultSizeY.toStdString() << "\t\t" << m_rawDefaultSizeY.toStdString() << std::endl;

	if (!m_rawDefaultSizeZ.isEmpty())
		config << rawDefaultSizeZ.toStdString() << "\t\t" << m_rawDefaultSizeZ.toStdString() << std::endl;

	if (!m_rawDefaultSizeVox.isEmpty())
		config << rawDefaultSizeVox.toStdString() << "\t" << m_rawDefaultSizeVox.toStdString() << std::endl;

	config.close();

	return true;
}


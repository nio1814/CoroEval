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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <QString>

#include <cstddef>


class Settings
{
private:
    Settings();
    
	size_t m_numNormalsAtPoint;
    size_t m_sizeNormal;
	double m_minRelMagnitude;
	double m_plateauThreshold;
	double m_minVesselDiameter;

    double m_sharpnessEvalInterval;

	bool   m_doPointCheck;
	bool   m_doNormalCheck;
	bool   m_useEllipse;

	bool   m_showMIP;
	bool   m_showNormals;

	QString m_defaultDir;
	bool   m_rememberDir;

	QString m_rawDefaultSizeX;
	QString m_rawDefaultSizeY;
	QString m_rawDefaultSizeZ;
	QString m_rawDefaultSizeVox;

public:
    static Settings& getSettings();

    void   setNumNormalsAtPoint(size_t val);
    size_t getNumNormalsAtPoint() const;

    void   setSizeNormal(size_t val);
    size_t getSizeNormal() const;

	void   setMinRelMagnitude(double val);
	double getMinRelMagniude() const;

	void   setPlateauThreshold(double val);
	double getPlateauThreshold() const;

	void   setMinVesselDiameter(double val);
	double getMinVesselDiameter() const;

    void   setSharpnessEvalInterval(double val);
    double getSharpnessEvalInterval() const;

	void   setDoPointCheck(bool val);
	bool   getDoPointCheck() const;

	void   setDoNormalCheck(bool val);
	bool   getDoNormalCheck() const;

	void   setUseEllipse(bool val);
	bool   getUseEllipse() const;

	void   setShowMIP(bool val);
	bool   getShowMIP() const;

	void   setShowNormals(bool val);
	bool   getShowNormals() const;

	void   setDefaultDir(const QString& val);
	const  QString& getDefaultDir() const;

	void   setRememberLastDir(bool val);
	bool   getRememberLastDir() const;

	const QString& getRawDefaultSizeX() const;
	const QString& getRawDefaultSizeY() const;
	const QString& getRawDefaultSizeZ() const;
	const QString& getRawDefaultSizeVox() const;

	bool   load(const QString& filename);
	bool   save(const QString& filename);
};

#endif // SETTINGS_H

